/*----------------------------------------------------------------------------- 
 * vlsSMLM Software
 * 
 * Copyright (C) 2014 Matthieu Palayret
 * Department of Chemistry
 * University of Cambridge, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

package smlm.vLS;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import smlm.util.SRutil;
import smlm.Params;
import smlm.fitting.FittingPeakFit;
import smlm.util.ResultsTableMt;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.Opener;
import ij.measure.CurveFitter;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.Zoom;
import ij.plugin.filter.MaximumFinder;
import ij.process.ByteProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;

public class Bead_Stack_Analysis implements PlugIn, DialogListener {

	Parameters param;
	private boolean cancel = false;

	CurveFitter reg;
	int guiNumber = 0;
	ImagePlus zAxis;
	int[] focalPlane = new int[2];
	ResultsTableMt rt;
	ImagePlus inFocus;
	ImagePlus outOfFocus;
	ImagePlus sigmaThreshPlot;
	ImagePlus iThreshPlot;
	int i0_ = ResultsTableMt.AMPLITUDE;// ResultsTableMt.I0;

	private double vLSWidth;
	private double sigma0, i0;

	public Bead_Stack_Analysis() {
	}

	class Parameters {
		double stepSize = 10; // nm
		int framesPerStep = 10;
		int fitting = -1;
		static final int NOFit = 0;
		static final int PeakFit = 1;
		static final int M2LE = 2;
		static final int QuickPALM = 3;
		static final int rapidSTORM = 4;
		static final int ThunderSTORM = 5;
		final String[] fitCodes = { "Previously fit", "PeakFit", "M2LE", "QuickPALM", "rapidSTORM", "ThunderSTORM" };
		int totalFrameIni;
		int optionsLength;

		String fileDirName;
		Params param;

		Parameters(String rootDirName, String fileName) {
			param = new Params(rootDirName, fileName);
			(new java.io.File(param.fileDirName)).mkdir();
			(new java.io.File(param.fileDirName + File.separator + "Temp")).mkdir();

			fileDirName = param.fileDirName;
		}

		void getParameters() {
			// Asking the parameters to the user
			GenericDialog gui = new GenericDialog("3D-PSF analysis");

			gui.addMessage("Parameters:");
			gui.addNumericField("Size of a pixel (nm):", 110, 0);
			gui.addNumericField("NA:", 1.49, 2);
			gui.addNumericField("Emission peak (nm):", 641, 0);

			gui.addMessage("");
			gui.addNumericField("z-step size (nm):", 10, 1);
			gui.addNumericField("Number of frames per step:", 10, 0);

			gui.addMessage("");

			gui.addNumericField("Gain (ADU/photon):", 250.0D / 11.5D, 2);
			gui.addNumericField("Camera offset (ADU):", 500, 0);

			gui.addMessage("");
			getAvailableFit();
			if (optionsLength > 0) {
				String[] FitCodesTemp = new String[2];
				for (int i = 0; i < 2; i++)
					FitCodesTemp[i] = fitCodes[i];
				gui.addChoice("Fitting method", FitCodesTemp, FitCodesTemp[0]);
			} else {
				String[] FitCodesTemp = new String[1];
				FitCodesTemp[0] = fitCodes[0];
				gui.addChoice("Fitting method", FitCodesTemp, FitCodesTemp[0]);
			}

			gui.showDialog();
			if(gui.wasCanceled()) {
				cancel = true;
				return;
			}

			double pixelSize = gui.getNextNumber();
			double na = gui.getNextNumber();
			double lambda = gui.getNextNumber();
			stepSize = gui.getNextNumber();
			framesPerStep = (int) gui.getNextNumber();
			double gain = gui.getNextNumber();
			double offset = gui.getNextNumber();
			if (optionsLength > 0)
				fitting = gui.getNextChoiceIndex();
			else
				fitting = gui.getNextChoiceIndex() + 1;

			param.hasBeenInitialised = true;
			param.pixelSize = pixelSize;
			param.na = na;
			param.lambda = lambda;
			param.fitting = Params.Fitting.FullGaussianFit;
			param.group = Params.Grouping.NoGrouping;
			param.rebuilding = Params.Rebuilding.NoRebuilt;
			param.rebuiltScale = 5;
			param.update();
			param.photonCountMode = false;// true;
			param.scaleBar = true;
			param.debug = false;
			param.checkMovie = false;
			param.preBleaching = true;
			param.gain = 1;
			param.emGain = gain * Params.ConversionGain[param.camera][param.gain];
			param.offset = offset;
		}
	}

	void newParameters(Params param, ResultsTableMt rt) {
		this.param = new Parameters(param.rootDirName, param.fileDirName + ".tif");
		this.param.param = param;
		this.rt = rt;
	}

	@Override
	public void run(String arg0) {

		// Ask for the .tif movie
		String[] temp = SRutil.getAFile("z-stack of the beads", "C://", "movie.tif");
		if(temp==null) return;

		// Ask for the parameters (step size and number of frames per step)
		param = new Parameters(temp[0], temp[1]);
		param.getParameters();
		if(cancel) return;

		// Combine all the beads in an average 3D-PSF
		Opener opener = new Opener();
		ImagePlus movie = opener.openImage(temp[0] + File.separator + temp[1]);
		param.totalFrameIni = movie.getStackSize();
		param.param.roi = new Rectangle(0, 0, movie.getWidth(), movie.getHeight());

		// Sum all the frames from each step together (one frame per step)
		ImageStack oneFramePerStep = new ImageStack(movie.getWidth(), movie.getHeight());
		{
			for (int step = 0; step < movie.getStackSize() / param.framesPerStep; step++) {
				ZProjector zProj = new ZProjector(movie);
				if (param.framesPerStep > 1) {
					zProj.setMethod(ZProjector.AVG_METHOD);
					zProj.setStartSlice(step * param.framesPerStep);
					zProj.setStopSlice((step + 1) * param.framesPerStep - 1);
					zProj.doProjection();
					oneFramePerStep.addSlice(zProj.getProjection().getProcessor());
				} else {
					oneFramePerStep.addSlice(movie.getImageStack().getProcessor(step + 1));
				}
			}
		}
		ImagePlus movieProj = new ImagePlus("One frame per step", oneFramePerStep);

		// Produce the movie and orthogonal projections of the 3D-PSF
		ResultsTableMt rtMax = new ResultsTableMt();
		// Get the bead
		{
			ZProjector zProj = new ZProjector(movieProj);
			zProj.setMethod(ZProjector.AVG_METHOD);
			zProj.doProjection();
			ImageProcessor ip = zProj.getProjection().getProcessor();

			// Gaussian filter Cf Matlab bpass.m
			double lowNoise = 1;
			double highNoise = 4.0D * param.param.psfSigmaInt;
			ImageProcessor ip2 = SRutil.gaussianFilter(ip, lowNoise, highNoise, param.param);

			// Get the maxima
			ip2.setRoi(param.param.roi);
			ImageProcessor ipCropped = ip2.crop().duplicate();

			ImageStatistics stat = ipCropped.getStatistics();
			double backgroundLevelFiltered = stat.mean;
			double stdBackground = stat.stdDev;
			MaximumFinder max = new MaximumFinder();
			double temp1 = 3.0D;// 1.0D;
			double temp2 = 7.0D;
			ByteProcessor ByteIp = max.findMaxima(ip2.crop(), temp1 * stdBackground, 1.2D
					* backgroundLevelFiltered + temp2 * stdBackground, MaximumFinder.SINGLE_POINTS, false,
					false);// 3*StdBackground,1.2*BackgroundLevelFiltered+4*StdBackground

			ip.setRoi(param.param.roi);
			ipCropped = ip.crop().duplicate();
			double BackgroundLevel = ipCropped.getStatistics().mean;

			// Store the maxima
			double intensity;
			int roiWidth = 3;
			for (int iii = 0; iii < ByteIp.getWidth(); iii++) {
				for (int jjj = 0; jjj < ByteIp.getHeight(); jjj++) {
					if (ByteIp.get(iii, jjj) == 255) {
						intensity = 0;

						for (int ii = 0; ii < param.param.psfSigmaInt * 2 * roiWidth + 1; ii++) {
							for (int jj = 0; jj < param.param.psfSigmaInt * 2 * roiWidth + 1; jj++) {
								intensity += Math
										.max(ip.get(
												Math.min(
														Math.max((int) iii + param.param.roi.x - roiWidth
																* param.param.psfSigmaInt + ii, 0),
														ip.getWidth() - 1),
												Math.min(
														Math.max((int) jjj + param.param.roi.y - roiWidth
																* param.param.psfSigmaInt + jj, 0),
														ip.getHeight() - 1))
												- BackgroundLevel, 0);
							}
						}
						rtMax.incrementCounter();
						rtMax.addValue(ResultsTableMt.Xmax, iii + param.param.roi.x);
						rtMax.addValue(ResultsTableMt.Ymax, jjj + param.param.roi.y);
						rtMax.addValue(ResultsTableMt.Imax, intensity);

					}
				}
			}

			// IJ.log(Param.rootDirName+File.separator+"TableMax.txt");

			try {
				rtMax.saveAs(param.param.fileDirName + File.separator + "TableMax.txt");
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		// Draw a box and save the xy-z movie
		int boxSize = 15; // in pixels
		ImageStack sumStack = new ImageStack(2 * boxSize + 1, 2 * boxSize + 1);
		{
			movieProj.show();

			for (int psf = 0; psf < rtMax.getCounter(); psf++) {

				if (SRutil.closestNeighbour(rtMax, psf)[0] > boxSize
						&& rtMax.getValueAsDouble(ResultsTableMt.Xmax, psf) > boxSize
						&& rtMax.getValueAsDouble(ResultsTableMt.Xmax, psf) < movie.getWidth() - boxSize
						&& rtMax.getValueAsDouble(ResultsTableMt.Ymax, psf) > boxSize
						&& rtMax.getValueAsDouble(ResultsTableMt.Ymax, psf) < movie.getHeight() - boxSize) {
					movieProj.setRoi(new Rectangle((int) rtMax.getValueAsDouble(ResultsTableMt.Xmax, psf)
							- boxSize, (int) rtMax.getValueAsDouble(ResultsTableMt.Ymax, psf) - boxSize,
							2 * boxSize + 1, 2 * boxSize + 1));

					ImageStack psfSquare = new ImageStack(2 * boxSize + 1, 2 * boxSize + 1);
					for (int frame = 0; frame < movieProj.getImageStackSize(); frame++) {
						ImageProcessor ipTemp = movieProj.getStack().getProcessor(frame + 1);
						ipTemp.setRoi(new Rectangle((int) rtMax.getValueAsDouble(ResultsTableMt.Xmax, psf)
								- boxSize, (int) rtMax.getValueAsDouble(ResultsTableMt.Ymax, psf) - boxSize,
								2 * boxSize + 1, 2 * boxSize + 1));
						ipTemp = ipTemp.crop().duplicate();
						psfSquare.addSlice("", ipTemp, frame);

						if (sumStack.getSize() < psfSquare.getSize())
							sumStack.addSlice("", ipTemp.duplicate(), frame);
						else {
							for (int x = 0; x < sumStack.getWidth(); x++) {
								for (int y = 0; y < sumStack.getHeight(); y++) {
									sumStack.setVoxel(x, y, frame,
											sumStack.getVoxel(x, y, frame) + ipTemp.get(x, y));
								}
							}
						}
					}

					SRutil.saveTiff(new ImagePlus("", psfSquare), param.param.fileDirName + File.separator
							+ "Temp" + File.separator + "PSF_" + (psf + 1) + ".tif", true);
				}

			}

			SRutil.saveTiff(new ImagePlus("", sumStack), param.param.fileDirName + File.separator
					+ "PSFSumXY.tif", false);

		}

		// Fit the 3D-PSF with a Gaussian
		{
			double[] zAxisProfile = new double[sumStack.getSize()];
			double[] x = new double[zAxisProfile.length];
			for (int frame = 0; frame < sumStack.getSize(); frame++) {
				ImageProcessor temp0 = sumStack.getProcessor(frame + 1);
				temp0.setRoi(new Rectangle(boxSize - 2, boxSize - 2, 5, 5));
				zAxisProfile[frame] = SRutil.averaging(temp0.crop().getIntArray());
				x[frame] = frame;
			}

			// Search for the local maxima and minima
			double[] zAxisProfileFiltered = SRutil.runningAverage(zAxisProfile, 10);

			zAxisProfileFiltered = SRutil.derivative(zAxisProfileFiltered);
			zAxisProfileFiltered = SRutil.runningAverage(zAxisProfileFiltered, 5);

			ResultsTableMt maxMin = new ResultsTableMt();
			maxMin.incrementCounter();
			maxMin.addValue("Max", 0);
			maxMin.addValue("Min", 0);
			int min = maxMin.getColumnIndex("Min"), max = maxMin.getColumnIndex("Max");
			for (int frame = 1; frame < zAxisProfile.length; frame++) {
				if (zAxisProfileFiltered[frame - 1] * zAxisProfileFiltered[frame] > 0) { // Same
					// sign,
					// no
					// extremum
				} else if (zAxisProfileFiltered[frame - 1] == 0) {
					if (zAxisProfileFiltered[frame] == 0) { // Do nothing, the
						// extremum has been
						// dealt with
						// previously
					} else if (zAxisProfileFiltered[frame] > 0) {
						if (maxMin.getValueAsDouble(min, maxMin.getCounter() - 1) > 0)
							maxMin.incrementCounter();
						maxMin.addValue(min, frame - 1);
					} else {
						if (maxMin.getValueAsDouble(max, maxMin.getCounter() - 1) > 0)
							maxMin.incrementCounter();
						maxMin.addValue(max, frame - 1);
					}
				} else if (zAxisProfileFiltered[frame] == 0) { // Do nothing,
					// the extremum
					// will be dealt
					// with in next
					// cycle
				} else if (zAxisProfileFiltered[frame - 1] > 0) {
					if (maxMin.getValueAsDouble(max, maxMin.getCounter() - 1) > 0)
						maxMin.incrementCounter();
					maxMin.addValue(max, frame - 1);
				} else if (zAxisProfileFiltered[frame - 1] < 0) {
					if (maxMin.getValueAsDouble(min, maxMin.getCounter() - 1) > 0)
						maxMin.incrementCounter();
					maxMin.addValue(min, frame - 1);
				}
			}

			// Only fit the central maxima up to its adjacent minima (avoid
			// mis-fitting due to out of focus light)
			double focalMax = 0;
			double distanceTo = 100000;
			for (int row = 0; row < maxMin.getCounter(); row++) {
				if (Math.abs(maxMin.getValue("Max", row) - (double) zAxisProfile.length / 2.0D) < distanceTo) {
					distanceTo = Math.abs(maxMin.getValue("Max", row) - (double) zAxisProfile.length / 2.0D);
					focalMax = maxMin.getValue("Max", row);
				}
			}
			double focalMinDown = 0;
			double focalMinUp = zAxisProfile.length - 1;
			double distanceToDown = 100000;
			double distanceToUp = 100000;
			for (int row = 0; row < maxMin.getCounter(); row++) {
				if (focalMax > maxMin.getValue("Min", row)
						&& focalMax - maxMin.getValue("Min", row) < distanceToDown) {
					distanceToDown = focalMax - maxMin.getValue("Min", row);
					focalMinDown = maxMin.getValue("Min", row);
				} else if (focalMax < maxMin.getValue("Min", row)
						&& maxMin.getValue("Min", row) - focalMax < distanceToUp) {
					distanceToUp = maxMin.getValue("Min", row) - focalMax;
					focalMinUp = maxMin.getValue("Min", row);
				}
			}

			double[] zAxisProfileForReg = new double[(int) (focalMinUp - focalMinDown + 1)];
			double[] xForReg = new double[zAxisProfileForReg.length];
			for (int frame = 0; frame < zAxisProfileForReg.length; frame++) {
				zAxisProfileForReg[frame] = zAxisProfile[frame + (int) focalMinDown];
				xForReg[frame] = frame + focalMinDown;
			}

			reg = new CurveFitter(x/*ForReg*/, zAxisProfile/*ForReg*/);
			reg.doFit(CurveFitter.GAUSSIAN);
			
			// QDot 655 - 65 -- 300
			// Beads 550 - 50 -- 300
			// Beads 660 - 75 -- 300
			// Beads 700 - 50 -- 300
			// Tom beads 505 - 180 -- 300
			//reg.getParams()[2] = (int)((double)(sumStack.getSize())/ 2.0 + 0.5); // TODO Comment this in general, except if not using beads but single fluorophores!!
			//reg.getParams()[3] = 283.0/param.stepSize;
			
		}

		// Interactively ask for a depth of the focal plane and define it
		{
			guiNumber = 1;
			NonBlockingGenericDialog gui = new NonBlockingGenericDialog("Depth of the focal plane");

			gui.addMessage("Fit of the z intensity profile of the PSF:");
			gui.addMessage("Gaussian fit: A0/(sigma*sqrt(2*Pi)) * exp(-1/2*((x-mu)/sigma)^2");
			gui.addMessage("Focal plane (nm): mu = " + IJ.d2s(reg.getParams()[2] * param.stepSize, 2));
			gui.addMessage("Standard deviation (nm): sigma = "
					+ IJ.d2s(reg.getParams()[3] * param.stepSize, 2));
			gui.addMessage("Full width at half maximum (nm): FWHM = "
					+ IJ.d2s(2.0 * Math.sqrt(2.0 * Math.log(2)) * reg.getParams()[3] * param.stepSize, 2));
			gui.addMessage("Half width at 1/e^2 (nm): w = "
					+ IJ.d2s(2.0 * reg.getParams()[3] * param.stepSize, 2));

			gui.addSlider("Width selector (nm)", param.stepSize, 5 * reg.getParams()[3] * param.stepSize,
					reg.getParams()[3] * param.stepSize);

			gui.addDialogListener(this);

			gui.showDialog();

			if (gui.wasCanceled())
				updateZAxis(reg.getParams()[3]);
			else if (gui.wasOKed())
				updateZAxis(selectedZWidth);
			SRutil.saveTiff(zAxis, param.param.fileDirName + File.separator + "FocalPlane_Definition.tif",
					false);
		}

		// Fit all PSFs
		if (param.fitting == Parameters.PeakFit) {
			movie.show();
			movie.setRoi(new Rectangle(0, 0, movie.getWidth(), movie.getHeight()));

			FittingPeakFit fitting = new FittingPeakFit();
			fitting.setConfig(param.param, FittingPeakFit.VLS);
			fitting.fitImage(movie);
			rt = fitting.getResults();

			try {
				rt.saveAs(param.param.fileDirName + File.separator + "TableFit_PeakFit.txt");
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		// Build the parameter diagrams
		{
			if (param.fitting == Parameters.NOFit) {
				boolean[] fitCode = getAvailableFit();
				String[] options = new String[param.optionsLength];

				int i = 0;
				for (int j = 0; j < fitCode.length; j++) {
					if (fitCode[j]) {
						options[i] = param.fitCodes[j];
						i++;
					}
				}

				if (options.length > 1) {
					guiNumber = 3;
					NonBlockingGenericDialog gui = new NonBlockingGenericDialog("Fit code to be used:");
					gui.addChoice("Fit results to be used", options, options[0]);
					gui.showDialog();
					
					int index = 0;
					if(gui.wasOKed()) index = gui.getNextChoiceIndex();

					rt = ResultsTableMt.open2(param.param.fileDirName + File.separator + "TableFit_"
							+ options[index] + ".txt");
				} else
					rt = ResultsTableMt.open2(param.param.fileDirName + File.separator + "TableFit_"
							+ options[0] + ".txt");
			}
			
			// Change all units from Counts and pixels to photons and nm
			rt = SRutil.translateCountsToPhotons(rt, i0_, param.param, false);
			rt = SRutil.translateCountsToPhotons(rt, ResultsTableMt.NOISE, param.param, false);
			param.param.photonCountMode = true;
			// Change amplitude units to photons/um^2
			if(i0_ == ResultsTableMt.AMPLITUDE) {
				for(int row=0; row<rt.getCounter(); row++) {
					rt.setValue(i0_, row, rt.getValueAsDouble(i0_, row)/Math.pow(param.param.pixelSize*0.001, 2));
				}
			}
			for(int row=0; row<rt.getCounter(); row++) {
				rt.setValue(ResultsTableMt.SIGMAX, row, rt.getValueAsDouble(ResultsTableMt.SIGMAX, row)*param.param.pixelSize);
				rt.setValue(ResultsTableMt.SIGMAY, row, rt.getValueAsDouble(ResultsTableMt.SIGMAY, row)*param.param.pixelSize);
			}

			analyseTenPerTenFrame(false);
		}

		// Calculate and plot the certainty metrics
		// Interactively ask for thresholds
		{
			guiNumber = 2;
			NonBlockingGenericDialog gui = new NonBlockingGenericDialog("Choose thresholds");

			gui.addMessage("Thresholds in width and amplitude:");
			gui.addMessage("(Use '-1' to ignore the corresponding threshold)");
			gui.addSlider("Width threshold (nm)", 1, (int) ((250.0 / scaleSigma)),
					223);
			gui.addSlider("Amplitude (photons/" + IJ.micronSymbol + "m^2)", (int) (1.0 / scaleI0),
					(int) (200.0 / scaleI0), 1502.0);

			gui.addDialogListener(this);
			gui.setLocation(405, 112);

			gui.showDialog();

			if (gui.wasCanceled()) {
				updateSigmaThresholdPlot(2.03, 1502, false);
				updateI0ThresholdPlot(2.03, 1502, false);
				updateOverlayThresholdedParamPlot(2.03, 1502, false);
			} else if (gui.wasOKed()) {
				updateSigmaThresholdPlot(previousSigma0, previousI0, false);
				updateI0ThresholdPlot(previousSigma0, previousI0, false);
				updateOverlayThresholdedParamPlot(previousSigma0, previousI0, false);
			}
			SRutil.saveTiff(sigmaThreshPlot, param.fileDirName + File.separator + "7_ThresholdingSigma.tif",
					false);
			SRutil.saveTiff(iThreshPlot, param.fileDirName + File.separator + "7_ThresholdingIntensity.tif",
					false);

			// Transfer both chosen thresholds to the QuickPSF_Filter plugin
			// through ImageJ preferences class
			Prefs.set("VLS.widthThreshold", previousSigma0);
			Prefs.set("VLS.amplitudeThreshold", previousI0);
		}

	}

	private boolean[] getAvailableFit() {
		boolean[] fitCode = new boolean[param.fitCodes.length];
		param.optionsLength = 0;
		for (int i = 0; i < fitCode.length; i++) {
			fitCode[i] = new File(param.param.fileDirName + File.separator + "TableFit_" + param.fitCodes[i]
					+ ".txt").isFile();
			param.optionsLength += (fitCode[i] ? 1 : 0);
		}
		return fitCode;
	}

	@Override
	public boolean dialogItemChanged(GenericDialog gui, AWTEvent e) {

		if (guiNumber == 1) {
			vLSWidth = gui.getNextNumber() / param.stepSize;
		} else {
			sigma0 = gui.getNextNumber();
			i0 = gui.getNextNumber();
		}
		updatePlot();
		updateImage();
		return true;
	}

	private boolean plotLock = false;

	private synchronized boolean aquirePlotLock() {
		if (plotLock)
			return false;
		return plotLock = true;
	}

	private void updatePlot() {
		if (aquirePlotLock()) {
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged) {
							// Store the parameters to be processed
							double myVlsWidth = vLSWidth;
							double mySigma0 = sigma0;
							double myI0 = i0;
							// Do something with parameters
							if (guiNumber == 1)
								updateZAxis(myVlsWidth);
							else if (guiNumber == 2) {
								updateSigmaThresholdPlot(mySigma0, myI0, true);
								updateOverlayThresholdedParamPlot(mySigma0, myI0, true);
							}
							// Check if the parameters have changed again
							parametersChanged = (myVlsWidth != vLSWidth) || (mySigma0 != sigma0)
									|| (myI0 != i0);
						}
					} finally {
						// Ensure the running flag is reset
						plotLock = false;
					}
				}
			}).start();
		}
	}

	private boolean imageLock = false;

	private synchronized boolean aquireImageLock() {
		if (imageLock)
			return false;
		return imageLock = true;
	}

	private void updateImage() {
		if (aquireImageLock()) {
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged) {
							// Store the parameters to be processed
							double myVlsWidth = vLSWidth;
							double mySigma0 = sigma0;
							double myI0 = i0;
							// Do something with parameters
							if (guiNumber == 2)
								updateI0ThresholdPlot(mySigma0, myI0, true);
							// Check if the parameters have changed again
							parametersChanged = (myVlsWidth != vLSWidth) || (mySigma0 != sigma0)
									|| (myI0 != i0);
						}
					} finally {
						// Ensure the running flag is reset
						imageLock = false;
					}
				}
			}).start();
		}
	}

	double selectedZWidth;

	private void updateZAxis(double width) {
		selectedZWidth = width;

		Plot plot = new Plot("z-axis profile of the PSF", "frame number", "intensity (A.U.)");
		double[] x = reg.getXPoints();
		double[] y = new double[x.length];
		int yInFocus0 = (int) reg.getParams()[2] - (int) (width);
		int yInFocusEnd = (int) reg.getParams()[2] + (int) (width);
		focalPlane[0] = yInFocus0 * param.framesPerStep;
		focalPlane[1] = yInFocusEnd * param.framesPerStep;
		//IJ.log(""+focalPlane[0]+"  "+focalPlane[1]);
		double[] yInFocus = new double[2 * (int) (width) + 1];
		double[] xInFocus = new double[yInFocus.length];
		double yMax = 0;
		for (int frame = 0; frame < x.length; frame++) {
			yMax = Math.max(yMax, reg.getYPoints()[frame]);
			y[frame] = CurveFitter.f(CurveFitter.GAUSSIAN, reg.getParams(), x[frame]);
			if (x[frame] >= yInFocus0 && x[frame] <= yInFocusEnd) {
				yInFocus[(int) x[frame] - yInFocus0] = CurveFitter.f(CurveFitter.GAUSSIAN, reg.getParams(),
						x[frame]);
				xInFocus[(int) x[frame] - yInFocus0] = x[frame];
			}
		}
		plot.setLimits(0, x.length - 1, 0, yMax);

		plot.setColor(Color.BLACK);
		plot.addPoints(x, reg.getYPoints(), Plot.CIRCLE);
		plot.draw();

		plot.setLineWidth(2);
		plot.setColor(Color.RED);
		plot.addPoints(x, y, Plot.LINE);
		plot.setColor(Color.GREEN);
		plot.addPoints(xInFocus, yInFocus, Plot.LINE);
		plot.draw();

		plot.setColor(SRutil.lightGrey);
		plot.drawLine((int) reg.getParams()[2] - (int) reg.getParams()[3], 0, (int) reg.getParams()[2]
				- (int) reg.getParams()[3], yMax - 1);
		plot.drawLine((int) reg.getParams()[2] + (int) reg.getParams()[3], 0, (int) reg.getParams()[2]
				+ (int) reg.getParams()[3], yMax - 1);
		plot.draw();

		if (zAxis == null)
			zAxis = new ImagePlus("", plot.getProcessor());
		else
			zAxis.setProcessor(plot.getProcessor());
		zAxis.updateAndDraw();
		zAxis.show();
	}

	double scaleSigma;
	double scaleI0;
	final static int widthPlot = 251;
	final static int heightPlot = 251;

	private void analyseTenPerTenFrame(boolean densityPlot) {
		int totalFrame = (int) SRutil.getMax(rt, ResultsTableMt.FRAME);
		Plot[] widthVsIntensityPlot = new Plot[(int) ((double) totalFrame / (double) param.framesPerStep + 0.5D) + 1];
		Plot[] ellipticityPlot = new Plot[widthVsIntensityPlot.length];

		ShortProcessor widthOverPlane = new ShortProcessor(
				(int) ((double) param.totalFrameIni / param.framesPerStep),
				(int) (10 * param.param.pixelSize));
		ShortProcessor intensityOverPlane = new ShortProcessor(
				(int) ((double) param.totalFrameIni / param.framesPerStep), 10000);

		// Initialise the plots
		double[] i = new double[1];
		double[] sigmaX = new double[1];
		double[] sigmaY = new double[1];
		double[] frame = new double[1];
		for (int z = 0; z < widthVsIntensityPlot.length; z++) {
			widthVsIntensityPlot[z] = new Plot("", "Amplitude (photons/" + IJ.micronSymbol + "m^2)",
					"Width (stdev in nm)");
			if (widthVsIntensityPlot.length > 250)
				widthVsIntensityPlot[z].setFrameSize(528, 264);
			else
				widthVsIntensityPlot[z].setFrameSize(1056, 528);
			widthVsIntensityPlot[z].setLimits(0, 50000, 0, 10* param.param.pixelSize);
			widthVsIntensityPlot[z].setLineWidth(2);

			ellipticityPlot[z] = new Plot("", "SigmaX (nm)", "SigmaY (nm)");
			if (widthVsIntensityPlot.length > 250)
				ellipticityPlot[z].setFrameSize(264, 264);
			else
				ellipticityPlot[z].setFrameSize(528, 528);
			ellipticityPlot[z].setLimits(0, 10* param.param.pixelSize, 0, 10* param.param.pixelSize);
			ellipticityPlot[z].setLineWidth(2);

			i[0] = 0;
			sigmaX[0] = 0;
			sigmaY[0] = 0;
			widthVsIntensityPlot[z].setColor(Color.GREEN);
			widthVsIntensityPlot[z].addPoints(i, sigmaX, Plot.DOT);

			ellipticityPlot[z].setColor(Color.WHITE);
			ellipticityPlot[z].addPoints(sigmaX, sigmaY, Plot.DOT);
		}

		Plot ellipticityPlotSum = new Plot("", "SigmaX (nm)", "SigmaY (nm)");
		ellipticityPlotSum.setFrameSize(528, 528);
		ellipticityPlotSum.setLimits(0, 10* param.param.pixelSize, 0, 10* param.param.pixelSize);
		ellipticityPlotSum.setLineWidth(2);

		// Fill in the plots
		for (int psf = 0; psf < rt.getCounter(); psf++) {
			if (rt.getValueAsDouble(ResultsTableMt.IS_FITTED, psf) != 0) {
				i[0] = rt.getValueAsDouble(i0_, psf);
				sigmaX[0] = rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf);
				sigmaY[0] = rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf);
				int z = (int) (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) / (double) param.framesPerStep);
				frame[0] = z;
				double ellipticity = Math.min(Math.abs(sigmaX[0] / sigmaY[0]),
						Math.abs(sigmaY[0] / sigmaX[0]));

				widthVsIntensityPlot[z].setColor(new Color((int) ((1.0D - ellipticity) * 255.0D),
						(int) ((1.0D - ellipticity) * 255.0D), (int) ((1.0D - ellipticity) * 255.0D)));
				widthVsIntensityPlot[z].addPoints(i, sigmaX, Plot.DOT);
				widthVsIntensityPlot[z].addPoints(i, sigmaY, Plot.DOT);

				if (sigmaX[0] < widthOverPlane.getHeight() - 1 && sigmaX[0] >= 0)
					widthOverPlane
							.set((int) frame[0], (int) (sigmaX[0]), widthOverPlane
									.get((int) frame[0], (int) (sigmaX[0])) + 1);
				if (sigmaY[0] < widthOverPlane.getHeight() - 1 && sigmaY[0] >= 0)
					widthOverPlane
							.set((int) frame[0], (int) (sigmaY[0]), widthOverPlane
									.get((int) frame[0], (int) (sigmaY[0])) + 1);

				if (i[0] < intensityOverPlane.getHeight() - 1 && i[0] >= 0)
					intensityOverPlane.set((int) frame[0], (int) (i[0]),
							intensityOverPlane.get((int) frame[0], (int) (i[0])) + 1);

				if (ellipticity < 1.0D / 1.7D)
					ellipticityPlotSum.setColor(Color.BLACK);
				else
					ellipticityPlotSum.setColor(SRutil.grey);
				ellipticityPlotSum.addPoints(sigmaX, sigmaY, Plot.DOT);

				if (ellipticity < 1.0D / 1.7D)
					ellipticityPlot[z].setColor(Color.BLACK);
				else
					ellipticityPlot[z].setColor(SRutil.grey);
				ellipticityPlot[z].addPoints(sigmaX, sigmaY, Plot.DOT);
			}
		}

		// Stack the plots and save them
		ImageStack widthVsIntensityStack = new ImageStack(widthVsIntensityPlot[0].getImagePlus().getWidth(),
				widthVsIntensityPlot[0].getImagePlus().getHeight());
		ImageStack ellipticityStack = new ImageStack(ellipticityPlot[0].getImagePlus().getWidth(),
				ellipticityPlot[0].getImagePlus().getHeight());
		for (int z = 0; z < widthVsIntensityPlot.length; z++) {
			widthVsIntensityPlot[z].draw();
			widthVsIntensityStack.addSlice("" + z, widthVsIntensityPlot[z].getImagePlus().getProcessor(), z);
			ellipticityPlot[z].draw();
			ellipticityStack.addSlice("" + z, ellipticityPlot[z].getImagePlus().getProcessor(), z);
		}
		ellipticityPlotSum.draw();
		intensityOverPlane.flipVertical();
		widthOverPlane.flipVertical();

		SRutil.saveTiff(new ImagePlus("", ellipticityStack), param.fileDirName + File.separator
				+ "5_Ellipticity_Stack.tif", true);
		SRutil.saveTiff(new ImagePlus("", ellipticityPlotSum.getProcessor()), param.fileDirName
				+ File.separator + "5_Ellipticity.tif", true);
		SRutil.saveTiff(new ImagePlus("", widthVsIntensityStack), param.fileDirName + File.separator
				+ "4_Amplitude-vs-Stdev.tif", true);
		SRutil.saveTiff(
				new ImagePlus("", intensityOverPlane.resize(10 * intensityOverPlane.getWidth(),
						intensityOverPlane.getHeight() / 10)), param.fileDirName + File.separator
						+ "4_AmplitudeOverZPlane.tif", true);
		SRutil.saveTiff(
				new ImagePlus("", widthOverPlane.resize(10 * widthOverPlane.getWidth(),
						1 * widthOverPlane.getHeight())), param.fileDirName + File.separator
						+ "4_WidthOverZPlane.tif", true);

		// Get the SigmaThreshold
		ByteProcessor sigmaThresh = new ByteProcessor(widthPlot, heightPlot);
		ByteProcessor sigmaThreshFocal = new ByteProcessor(widthPlot, heightPlot);
		ByteProcessor sigmaThreshOutOfFocus = new ByteProcessor(widthPlot, heightPlot);

		// TODO
		scaleSigma = Math.min(
				100,
				(double) sigmaThresh.getWidth()
						/ (SRutil.getMean(rt, ResultsTableMt.SIGMAX) + 2.0*SRutil.getIQR(rt, ResultsTableMt.SIGMAX)));// 100.0D;
		scaleI0 = Math.min(1,
				(double) sigmaThresh.getWidth() / (SRutil.getMean(rt, i0_) + 5.0*SRutil.getIQR(rt, i0_)));

		for (int psf = 0; psf < rt.getCounter(); psf++) {
			if (rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) < (double) (sigmaThresh.getHeight() - 1)
					/ scaleSigma
					&& rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) < (double) (sigmaThresh.getHeight() - 1)
							/ scaleSigma
					&& rt.getValueAsDouble(i0_, psf) < (double) (sigmaThresh.getWidth() - 1) / scaleI0
					&& Math.min(
							Math.min(rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
									rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf)),
							rt.getValueAsDouble(i0_, psf)) > 0) {
				sigmaThresh
						.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
								(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D),
								sigmaThresh.get(
										(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
										(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D)) + 1);
				sigmaThresh
						.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
								(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D),
								sigmaThresh.get(
										(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
										(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D)) + 1);

				if (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) >= focalPlane[0]
						&& rt.getValueAsDouble(ResultsTableMt.FRAME, psf) < focalPlane[1]) {
					sigmaThreshFocal
							.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
									(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D),
									sigmaThreshFocal.get(
											(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
											(int) (scaleSigma
													* rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D)) + 1);
					sigmaThreshFocal
							.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
									(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D),
									sigmaThreshFocal.get(
											(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
											(int) (scaleSigma
													* rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D)) + 1);
				} else {
					sigmaThreshOutOfFocus
							.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
									(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D),
									sigmaThreshOutOfFocus.get(
											(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
											(int) (scaleSigma
													* rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5D)) + 1);
					sigmaThreshOutOfFocus
							.set((int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
									(int) (scaleSigma * rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D),
									sigmaThreshOutOfFocus.get(
											(int) (rt.getValueAsDouble(i0_, psf) * scaleI0 + 0.5D),
											(int) (scaleSigma
													* rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5D)) + 1);
				}
			}
		}
		sigmaThresh.flipVertical();
		// SRutil.saveTiff(new ImagePlus("All", sigmaThresh),
		// Param.FileDirName+File.separator+"4_Intensity-Vs-Stdev_all.tif",
		// true);
		sigmaThreshFocal.flipVertical();
		// SRutil.saveTiff(new ImagePlus("In focus", sigmaThreshFocal),
		// Param.FileDirName+File.separator+"4_Intensity-Vs-Stdev_InFocus.tif",
		// true);
		sigmaThreshOutOfFocus.flipVertical();
		// SRutil.saveTiff(new ImagePlus("Out of focus", sigmaThreshOutOfFocus),
		// Param.FileDirName+File.separator+"4_Intensity-Vs-Stdev_OutOfFocus.tif",
		// true);

		ByteProcessor blank = new ByteProcessor(widthPlot, heightPlot);
		blank.setRoi(new Rectangle(0, 0, blank.getWidth(), blank.getHeight()));
		blank.setColor(Color.WHITE);
		blank.fill();
		blank.setMinAndMax(0, 255);

		ImageCalculator ic = new ImageCalculator();
		ImageProcessor sigmaThreshFocalG = ic.run("Substract create", new ImagePlus("", blank),
				new ImagePlus("", sigmaThreshFocal)).getProcessor();
		int min = SRutil.getMin(sigmaThreshFocalG.getIntArray());
		sigmaThreshFocalG.setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		ImagePlus[] temp = { new ImagePlus("", sigmaThreshFocalG), new ImagePlus("", blank),
				new ImagePlus("", sigmaThreshFocalG) };
		inFocus = RGBStackMerge.mergeChannels(temp, true);
		ImageConverter icon = new ImageConverter(inFocus);
		ImageConverter.setDoScaling(false);
		icon.convertRGBStackToRGB();
		inFocus.getProcessor().setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		inFocus.setOverlay(getScaleBarOverlay(sigmaThresh.getWidth(), sigmaThresh.getHeight()));
		inFocus.setHideOverlay(false);
		inFocus.updateAndDraw();
		inFocus.show();
		inFocus.getWindow().setLocation(536, 355);
		new Zoom().run("in");
		SRutil.saveTiff(inFocus, param.fileDirName + File.separator + "4_Amplitude-Vs-Stdev_InFocus.tif",
				false);

		ic = new ImageCalculator();
		ImageProcessor sigmaThreshOutOfFocusR = ic.run("Substract create", new ImagePlus("", blank),
				new ImagePlus("", sigmaThreshOutOfFocus)).getProcessor();
		min = SRutil.getMin(sigmaThreshOutOfFocusR.getIntArray());
		sigmaThreshOutOfFocusR.setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		ImagePlus[] Temp2 = { new ImagePlus("", blank), new ImagePlus("", sigmaThreshOutOfFocusR),
				new ImagePlus("", sigmaThreshOutOfFocusR) };
		outOfFocus = RGBStackMerge.mergeChannels(Temp2, true);
		icon = new ImageConverter(outOfFocus);
		icon.convertRGBStackToRGB();
		outOfFocus.getProcessor().setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		outOfFocus.setOverlay(getScaleBarOverlay(sigmaThresh.getWidth(), sigmaThresh.getHeight()));
		outOfFocus.setHideOverlay(false);
		outOfFocus.updateAndDraw();
		outOfFocus.show();
		outOfFocus.getWindow().setLocation(209, 355);
		new Zoom().run("in");
		SRutil.saveTiff(outOfFocus, param.fileDirName + File.separator
				+ "4_Amplitude-Vs-Stdev_OutOfFocus.tif", false);

	}

	Overlay getScaleBarOverlay(int widthSigmaThresh, int heightSigmaThresh) {
		// Scale bars
		int scaleBarWidth = 3;
		Rectangle scaleBarI = new Rectangle(widthSigmaThresh - 5 - scaleBarWidth - 3 - 30, heightSigmaThresh
				- 5 - scaleBarWidth - 3, 30, scaleBarWidth);
		Rectangle scaleBarW = new Rectangle(widthSigmaThresh - 5 - scaleBarWidth - 3, heightSigmaThresh - 5
				- scaleBarWidth - 3 - 30, scaleBarWidth, 30);

		Overlay ol = new Overlay();
		ImageRoi imgRoi2 = new ImageRoi(scaleBarI.x, scaleBarI.y, new ShortProcessor(scaleBarI.width,
				scaleBarI.height));
		imgRoi2.setOpacity(1);
		imgRoi2.setStrokeColor(Color.BLACK);
		imgRoi2.setStrokeWidth(2);
		imgRoi2.setFillColor(Color.BLACK);
		ol.add(imgRoi2);
		TextRoi label1 = new TextRoi(scaleBarI.x - 20, scaleBarI.y - 10, "" + (int) (30.0D / scaleI0)
				+ " photons/" + IJ.micronSymbol + "m^2", new Font(Font.SANS_SERIF, Font.PLAIN, 8));
		label1.setStrokeColor(Color.BLACK);
		label1.setNonScalable(false);
		ol.add(label1);
		ImageRoi imgRoi3 = new ImageRoi(scaleBarW.x, scaleBarW.y, new ShortProcessor(scaleBarW.width,
				scaleBarW.height));
		imgRoi3.setOpacity(1);
		imgRoi3.setStrokeColor(Color.BLACK);
		imgRoi3.setStrokeWidth(2);
		imgRoi3.setFillColor(Color.BLACK);
		ol.add(imgRoi3);
		TextRoi label2 = new TextRoi(scaleBarW.x - 25, scaleBarW.y + scaleBarW.height / 2 - 6, ""
				+ (int) (30.0D / scaleSigma) + " nm", new Font(Font.SANS_SERIF,
				Font.PLAIN, 8));
		label2.setStrokeColor(Color.BLACK);
		label2.setNonScalable(false);
		ol.add(label2);
		ol.setLabelColor(Color.BLACK);
		return ol;
	}

	double previousSigma0 = 0;
	double previousI0 = 0;

	private void updateSigmaThresholdPlot(double sigma0, double i0, boolean small) {
		if (sigma0 * scaleSigma < 1 || i0 * scaleI0 < 1)
			return;

		if (sigmaThresholdedPlot == null || !small)
			buildSigmaThresholdedPlot(small);

		if (sigmaThreshPlot == null) {
			sigmaThreshPlot = new ImagePlus("Width thresholding", sigmaThresholdedPlot.getProcessor());
			sigmaThreshPlot.show();
			sigmaThreshPlot.getWindow().setLocation(864, 35);
		}
		if (!small)
			sigmaThreshPlot.setProcessor(sigmaThresholdedPlot.getProcessor());
		Overlay ol = new Overlay();
		ImageRoi roi = new ImageRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (sigmaThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (250.0D / scaleSigma) * sigma0) - 1, Plot.TOP_MARGIN, new ShortProcessor(4,
				sigmaThreshPlot.getHeight() - Plot.BOTTOM_MARGIN - Plot.TOP_MARGIN));
		roi.setFillColor(Color.BLACK);
		roi.setStrokeColor(Color.BLACK);
		roi.setStrokeWidth(1);
		ol.add(roi);
		TextRoi label1 = new TextRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (sigmaThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (250.0 / scaleSigma) * sigma0) - 1 + 5, Plot.TOP_MARGIN + 20, IJ.d2s(
				sigmaBlue[(int) (sigma0 * scaleSigma)] * 100.0, 1) + "%");
		label1.setStrokeColor(SRutil.royalBlue);
		label1.setNonScalable(false);
		ol.add(label1);
		TextRoi label2 = new TextRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (sigmaThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (250.0D / scaleSigma) * sigma0) - 1 + 5, Plot.TOP_MARGIN + 5, IJ.d2s(
				sigmaOrange[(int) (sigma0 * scaleSigma)] * 100.0, 1) + "%");
		label2.setStrokeColor(SRutil.yellow);
		label2.setNonScalable(false);
		ol.add(label2);
		ol.setLabelColor(Color.BLACK);

		sigmaThreshPlot.setOverlay(ol);
		sigmaThreshPlot.setHideOverlay(false);
		sigmaThreshPlot.updateAndDraw();
		sigmaThreshPlot.show();
		if (!small)
			for (int i = 0; i < 2; i++)
				sigmaThreshPlot.getWindow().getCanvas()
						.zoomOut(sigmaThreshPlot.getWidth() / 4, sigmaThreshPlot.getHeight() / 4);
	}

	private void updateI0ThresholdPlot(double sigma0, double i0, boolean small) {
		if (sigma0 * scaleSigma < 1 || i0 * scaleI0 < 1)
			return;

		if (iThresholdedPlot == null || sigma0 != previousSigma0 || !small)
			buildIThresholdedPlot(sigma0, small);

		if (iThreshPlot == null) {
			iThreshPlot = new ImagePlus("Amplitude thresholding", iThresholdedPlot.getProcessor());
			iThreshPlot.show();
			iThreshPlot.getWindow().setLocation(864, 413);
		}
		if (sigma0 != previousSigma0 || !small)
			iThreshPlot.setProcessor(iThresholdedPlot.getProcessor());

		Overlay ol = new Overlay();
		ImageRoi roi = new ImageRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (iThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (200.0D / scaleI0) * i0) - 1, Plot.TOP_MARGIN, new ShortProcessor(4,
				iThreshPlot.getHeight() - Plot.BOTTOM_MARGIN - Plot.TOP_MARGIN));
		roi.setFillColor(Color.BLACK);
		roi.setStrokeColor(Color.BLACK);
		roi.setStrokeWidth(1);
		ol.add(roi);
		TextRoi label1 = new TextRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (iThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (200.0D / scaleI0) * i0) - 1 + 5, Plot.TOP_MARGIN + 20, IJ.d2s(
				iBlue[(int) (i0 * scaleI0)] * 100.0, 1) + "%");
		label1.setStrokeColor(SRutil.royalBlue);
		label1.setNonScalable(false);
		ol.add(label1);
		TextRoi label2 = new TextRoi(Plot.LEFT_MARGIN
				+ (int) ((double) (iThreshPlot.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (200.0D / scaleI0) * i0) - 1 + 5, Plot.TOP_MARGIN + 5, IJ.d2s(
				iOrange[(int) (i0 * scaleI0)] * 100.0, 1) + "%");
		label2.setStrokeColor(SRutil.yellow);
		label2.setNonScalable(false);
		Roi.setColor(SRutil.royalBlue);
		ol.add(label2);
		ol.setLabelColor(Color.BLACK);

		iThreshPlot.setOverlay(ol);
		// IThreshPlot.setOverlay(new Roi(new
		// Rectangle(Plot.LEFT_MARGIN+(int)((double)(IThreshPlot.getWidth()-Plot.LEFT_MARGIN-Plot.RIGHT_MARGIN)/2000.0D*I0)-1,
		// Plot.TOP_MARGIN, 4,
		// IThreshPlot.getHeight()-Plot.BOTTOM_MARGIN-Plot.TOP_MARGIN)), new
		// Color(0, 0, 0), 1, new Color(0, 0, 0));
		iThreshPlot.setHideOverlay(false);
		iThreshPlot.updateAndDraw();
		iThreshPlot.show();
		if (!small)
			for (int i = 0; i < 2; i++)
				iThreshPlot.getWindow().getCanvas()
						.zoomOut(iThreshPlot.getWidth() / 4, iThreshPlot.getHeight() / 4);

		previousSigma0 = sigma0;
		previousI0 = i0;
	}

	private void updateOverlayThresholdedParamPlot(double sigma0, double i0, boolean small) {
		if (sigma0 * scaleSigma < 1 || i0 * scaleI0 < 1)
			return;

		Overlay ol = getScaleBarOverlay(inFocus.getWidth(), inFocus.getHeight());
		ImageRoi imgRoi1 = new ImageRoi(0, 0, new ShortProcessor(inFocus.getWidth(), inFocus.getHeight()
				- (int) (sigma0 * scaleSigma)));
		imgRoi1.setOpacity(0.5);
		imgRoi1.setStrokeColor(Color.BLACK);
		imgRoi1.setStrokeWidth(2);
		imgRoi1.setFillColor(Color.BLACK);
		ol.add(imgRoi1);

		ImageRoi imgRoi2 = new ImageRoi(0, inFocus.getHeight() - (int) (sigma0 * scaleSigma),
				new ShortProcessor((int) (i0 * scaleI0), (int) (sigma0 * scaleSigma)));
		imgRoi2.setOpacity(0.5);
		imgRoi2.setStrokeColor(Color.BLACK);
		imgRoi2.setStrokeWidth(2);
		imgRoi2.setFillColor(Color.BLACK);
		ol.add(imgRoi2);

		inFocus.setOverlay(ol);
		inFocus.setHideOverlay(false);
		inFocus.updateAndDraw();
		inFocus.show();

		outOfFocus.setOverlay(ol);
		outOfFocus.setHideOverlay(false);
		outOfFocus.updateAndDraw();
		outOfFocus.show();
	}

	Plot sigmaThresholdedPlot;
	double[] sigmaBlue;
	double[] sigmaOrange;

	private void buildSigmaThresholdedPlot(boolean small) {
		double[] sigmaThresh = new double[heightPlot];

		double[] greenIn = new double[sigmaThresh.length];
		double[] redIn = new double[sigmaThresh.length];
		double[] greenTotal = new double[sigmaThresh.length];

		double[] confidence = new double[sigmaThresh.length];
		double[] recall = new double[sigmaThresh.length];

		for (int sigmaThr = 0; sigmaThr < sigmaThresh.length; sigmaThr++) {
			sigmaThresh[sigmaThr] = (double) sigmaThr / scaleSigma;

			for (int psf = 0; psf < rt.getCounter(); psf++) {
				if (rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) < sigmaThresh[sigmaThr]
						&& rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) < sigmaThresh[sigmaThr]
						&& Math.min(
								Math.min(rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
										rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf)),
								rt.getValueAsDouble(i0_, psf)) > 0) {
					if (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) >= focalPlane[0]
							&& rt.getValueAsDouble(ResultsTableMt.FRAME, psf) < focalPlane[1])
						greenIn[sigmaThr]++;
					else
						redIn[sigmaThr]++;
				}
				if (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) >= focalPlane[0]
						&& rt.getValueAsDouble(ResultsTableMt.FRAME, psf) < focalPlane[1])
					greenTotal[sigmaThr]++;
			}

			confidence[sigmaThr] = greenIn[sigmaThr] / (greenIn[sigmaThr] + redIn[sigmaThr]);
			recall[sigmaThr] = greenIn[sigmaThr] / greenTotal[sigmaThr];
		}

		int lessThan1pcPSFs = 0;
		while (lessThan1pcPSFs < sigmaThresh.length && recall[(int) lessThan1pcPSFs] < 0.01)
			lessThan1pcPSFs++;

		Plot plot = new Plot("", "Threshold in width (nm)",
				"Confidence (orange) and recall (navy blue) ratios");
		if (small)
			plot.setFrameSize(528, 264);
		else
			plot.setFrameSize(1056, 528);
		plot.setLimits(0, sigmaThresh[sigmaThresh.length - 1], 0, 1);
		plot.getProcessor().setRoi(
				new Rectangle(Plot.LEFT_MARGIN + 1, Plot.TOP_MARGIN + 1, (int) ((double) (plot.getProcessor()
						.getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
						/ (sigmaThresh[sigmaThresh.length - 1] - 0) * sigmaThresh[lessThan1pcPSFs]) - 1, plot
						.getProcessor().getHeight() - Plot.TOP_MARGIN - Plot.BOTTOM_MARGIN - 1));
		plot.getProcessor().setColor(SRutil.lightGrey);
		plot.getProcessor().fill();
		plot.setLineWidth(2);
		plot.setColor(SRutil.yellow);
		plot.addPoints(sigmaThresh, confidence, Plot.LINE);
		plot.setColor(SRutil.royalBlue);
		plot.addPoints(sigmaThresh, recall, Plot.LINE);
		plot.draw();

		sigmaThresholdedPlot = plot;
		sigmaBlue = recall;
		sigmaOrange = confidence;

		SRutil.saveVector(confidence, param.fileDirName + File.separator + "11_Width_Confidence.txt");
		SRutil.saveVector(recall, param.fileDirName + File.separator + "11_Width_Recall.txt");
	}

	Plot iThresholdedPlot;
	double[] iBlue;
	double[] iOrange;

	private void buildIThresholdedPlot(double sigma0, boolean small) {
		double sigmaThr = sigma0;
		double[] iThresh = new double[widthPlot];

		double[] greenIn = new double[iThresh.length];
		double[] redIn = new double[iThresh.length];
		double[] greenTotal = new double[iThresh.length];

		double[] confidence = new double[iThresh.length];
		double[] recall = new double[iThresh.length];
		for (int iThr = 0; iThr < iThresh.length; iThr++) {
			iThresh[iThr] = (double) iThr / scaleI0;

			for (int psf = 0; psf < rt.getCounter(); psf++) {
				if (rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) < sigmaThr
						&& rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) < sigmaThr
						&& rt.getValueAsDouble(i0_, psf) > iThr / scaleI0
						&& Math.min(
								Math.min(rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
										rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf)),
								rt.getValueAsDouble(i0_, psf)) > 0) {
					if (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) >= focalPlane[0]
							&& rt.getValueAsDouble(ResultsTableMt.FRAME, psf) < focalPlane[1])
						greenIn[iThr]++;
					else
						redIn[iThr]++;
				}
				if (rt.getValueAsDouble(ResultsTableMt.FRAME, psf) >= focalPlane[0]
						&& rt.getValueAsDouble(ResultsTableMt.FRAME, psf) < focalPlane[1])
					greenTotal[iThr]++;
			}

			confidence[iThr] = greenIn[iThr] / (greenIn[iThr] + redIn[iThr]);
			recall[iThr] = greenIn[iThr] / greenTotal[iThr];
		}

		int lessThan1pcPSFs = iThresh.length - 1;
		while (lessThan1pcPSFs > 0 && recall[lessThan1pcPSFs] < 0.01D)
			lessThan1pcPSFs--;

		Plot plot = new Plot("", "Threshold in amplitude (photons/" + IJ.micronSymbol + "m^2)",
				"Confidence (orange)  and recall (navy blue) ratios");
		if (small)
			plot.setFrameSize(528, 264);
		else
			plot.setFrameSize(1056, 528);
		plot.setLimits(0, iThresh[iThresh.length - 1], 0, 1);
		plot.getProcessor()
				.setRoi(new Rectangle(
						Plot.LEFT_MARGIN
								+ (int) ((double) (plot.getProcessor().getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
										/ (iThresh[iThresh.length - 1] - 0) * iThresh[lessThan1pcPSFs]),
						Plot.TOP_MARGIN + 1,
						plot.getProcessor().getWidth()
								- Plot.RIGHT_MARGIN
								- 1
								- Plot.LEFT_MARGIN
								- (int) ((double) (plot.getProcessor().getWidth() - Plot.LEFT_MARGIN - Plot.RIGHT_MARGIN)
										/ (iThresh[iThresh.length - 1] - 0) * iThresh[lessThan1pcPSFs]), plot
								.getProcessor().getHeight() - Plot.TOP_MARGIN - Plot.BOTTOM_MARGIN - 1));
		plot.getProcessor().setColor(SRutil.lightGrey);
		plot.getProcessor().fill();
		plot.setLineWidth(2);
		plot.setColor(SRutil.yellow);
		plot.addPoints(iThresh, confidence, Plot.LINE);
		plot.setColor(SRutil.royalBlue);
		plot.addPoints(iThresh, recall, Plot.LINE);
		plot.draw();

		iThresholdedPlot = plot;
		iBlue = recall;
		iOrange = confidence;

		SRutil.saveVector(confidence, param.fileDirName + File.separator + "11_Amplitude_Confidence.txt");
		SRutil.saveVector(recall, param.fileDirName + File.separator + "11_Amplitude_Recall.txt");
	}

}
