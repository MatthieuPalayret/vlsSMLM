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
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import smlm.util.SRutil;
import smlm.Params;
import smlm.fitting.FittingPeakFit;
import smlm.util.ResultsTableMt;


import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.settings.ResultsSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.Opener;
import ij.plugin.ImageCalculator;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

public class QuickPSF_Filter implements PlugIn, DialogListener {

	Parameters param;
	private boolean cancel = false;
	ResultsTableMt rtFit;
	ShortProcessor srRebuilt;
	int xmin;
	int ymin;
	ImagePlus impRebuilt;
	ImagePlus parameterPlot;
	private Bead_Stack_Analysis paramPlots;
	int i0_ = ResultsTableMt.AMPLITUDE;// ResultsTableMt.I0;

	double sigma0;
	double sigma0Min;
	double i0;
	double e0;
	double averagePrecision;

	Rebuilding rebuilding = Rebuilding.QuickMt;

	enum Rebuilding {
		QuickMt(ResultsImage.NONE, 0), Frame_Number(ResultsImage.FRAME_NUMBER, 1), Localisation_Dots(
				ResultsImage.LOCALISATIONS, 2), Localisation_PSFs(ResultsImage.LOCALISATIONS_PRECISION, 3), Localisation_Av_PSF(
				ResultsImage.LOCALISATIONS_AV_PRECISION, 4), Diffraction_Limited(ResultsImage.PSF, 5), Intensity_Dots(
				ResultsImage.SIGNAL_INTENSITY, 6), Intensity_PSFs(ResultsImage.SIGNAL_PRECISION, 7), Intensity_Av_PSF(
				ResultsImage.SIGNAL_AV_PRECISION, 8);

		private ResultsImage code;
		private int index;

		private Rebuilding(ResultsImage code, int index) {
			this.code = code;
			this.index = index;
		}

		ResultsImage getCode() {
			return code;
		}

		int getIndex() {
			return index;
		};
	}

	final String[] rebuildingMode = { "Quick fixed PSFs", "Frame number", "Localisation dots",
			"Localisation PSFs", "Localisation Av. PSF", "Diffraction limited", "Intensity dots",
			"Intensity PSFs", "Intensity Av. PSF" };

	ResultsSettings rebuildingSettings;

	public QuickPSF_Filter() {
	}

	private class Parameters {
		int fitting = -1;
		static final int NOFit = 0;
		static final int PeakFit = 1;
		final String[] fitCodes = { "Previously fit", "PeakFit" };
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
			GenericDialog gui = new GenericDialog("Quick PSF filter:");

			gui.addMessage("Parameters:");
			gui.addNumericField("Size of a pixel (nm):", 110, 0);
			gui.addNumericField("NA:", 1.49, 2);
			gui.addNumericField("Emission peak (nm):", 641, 0);

			gui.addMessage("");

			gui.addNumericField("Gain (ADU/photon):", 250.0 / 11.5, 2);
			gui.addNumericField("Camera offset (ADU):", 500, 0);

			gui.addMessage("");
			getAvailableFit();
			if (optionsLength > 0)
				gui.addChoice("Fitting method", fitCodes, fitCodes[0]);
			else {
				String[] FitCodesTemp = new String[fitCodes.length - 1];
				for (int i = 0; i < FitCodesTemp.length; i++)
					FitCodesTemp[i] = fitCodes[i + 1];
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
			param.rebuiltScale = 10;
			param.update();
			param.photonCountMode = false;// true;
			param.scaleBar = true;
			param.debug = false;
			param.checkMovie = false;
			param.preBleaching = true;
			param.camera = 0;
			param.gain = 0;
			param.emGain = gain * Params.ConversionGain[param.camera][param.gain];
			param.offset = offset;

			rebuildingSettings = new ResultsSettings();
			rebuildingSettings.imageScale = param.rebuiltScale;
			rebuildingSettings.weightedImage = true;
			rebuildingSettings.equalisedImage = true;
			rebuildingSettings.imageRollingWindow = 0;
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

		// Find a Peakfit that has not been PeakFitToMP-ed
		if (fitCode[Parameters.PeakFit] == false) {
			String filePeakFit = null;
			if (SRutil.getFileName(param.param.fileDirName, ".results.xls") != null)
				filePeakFit = param.param.fileDirName + File.separator
						+ SRutil.getFileName(param.param.fileDirName, "results.xls");
			else if (SRutil.getFileName(param.param.rootDirName, param.param.fileName + ".results.xls") != null)
				filePeakFit = param.param.rootDirName + File.separator
						+ SRutil.getFileName(param.param.rootDirName, param.param.fileName + ".results.xls");
			if (filePeakFit != null) {
				rtFit = SRutil.changeRtNamesPeakFitToMP(filePeakFit, param.param);
				try {
					rtFit.saveAs(param.param.fileDirName + File.separator + "TableFit_"
							+ param.fitCodes[Parameters.PeakFit]);
				} catch (IOException e) {
					e.printStackTrace();
				}

				fitCode[Parameters.PeakFit] = true;
				param.optionsLength += (fitCode[Parameters.PeakFit] ? 1 : 0);
			}
		}

		return fitCode;
	}

	@Override
	public void run(String arg0) {
		// Ask for the .tif movie
		String[] temp = SRutil.getAFile("Select movie to analyse...", "C://", "movie.tif");
		if(temp==null) return;

		// Ask for the parameters (step size and number of frames per step)
		param = new Parameters(temp[0], temp[1]);
		param.getParameters();
		if(cancel) return;

		Opener opener = new Opener();
		ImagePlus movie;
		if (param.fitting == Parameters.NOFit)
			movie = opener.openImage(temp[0] + File.separator + temp[1], 1);
		else
			movie = opener.openImage(temp[0] + File.separator + temp[1]);
		param.param.roi = new Rectangle(0, 0, movie.getWidth(), movie.getHeight());

		// Fit all PSFs
		if (param.fitting == Parameters.PeakFit) {
			movie.show();
			movie.setRoi(new Rectangle(0, 0, movie.getWidth(), movie.getHeight()));

			FittingPeakFit fitting = new FittingPeakFit();
			fitting.setConfig(param.param, FittingPeakFit.VLS);
			fitting.fitImage(movie);
			rtFit = fitting.getResults();

			try {
				rtFit.saveAs(param.param.fileDirName + File.separator + "TableFit_PeakFit.txt");
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		// Get the ResultsTableMt RtFit from a file
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
					NonBlockingGenericDialog gui = new NonBlockingGenericDialog("Fit code to be used:");
					gui.addChoice("Fit results to be used", options, options[0]);
					gui.showDialog();

					// IJ.write(Param.Param.FileDirName+File.separator+"TableFit_"+Options[Gui.getNextChoiceIndex()]+".txt");
					rtFit = ResultsTableMt.open2(param.param.fileDirName + File.separator + "TableFit_"
							+ options[gui.getNextChoiceIndex()] + ".txt");
				} else
					rtFit = ResultsTableMt.open2(param.param.fileDirName + File.separator + "TableFit_"
							+ options[0] + ".txt");
			}

		}
		averagePrecision = SRutil.getPrecisionHistogram(rtFit, param.param, false, false);
		
		// Change all units from Counts and pixels to photons and nm
		rtFit = SRutil.translateCountsToPhotons(rtFit, i0_, param.param, false);
		rtFit = SRutil.translateCountsToPhotons(rtFit, ResultsTableMt.NOISE, param.param, false);
		param.param.photonCountMode = true;
		// Change amplitude units to photons/um^2
		if(i0_ == ResultsTableMt.AMPLITUDE) {
			for(int row=0; row<rtFit.getCounter(); row++) {
				rtFit.setValue(i0_, row, rtFit.getValueAsDouble(i0_, row)/Math.pow(param.param.pixelSize*0.001, 2));
			}
		}
		for(int row=0; row<rtFit.getCounter(); row++) {
			rtFit.setValue(ResultsTableMt.SIGMAX, row, rtFit.getValueAsDouble(ResultsTableMt.SIGMAX, row)*param.param.pixelSize);
			rtFit.setValue(ResultsTableMt.SIGMAY, row, rtFit.getValueAsDouble(ResultsTableMt.SIGMAY, row)*param.param.pixelSize);
		}

		// Order the RtFit by max(SigmaX, SigmaY)
		{
			getParameterPlot();

			ResultsTableMt rtFitOrdered = new ResultsTableMt();

			for (int psf = 0; psf < rtFit.getCounter(); psf++) {
				rtFit.setValue(
						ResultsTableMt.MAX_SIGMA,
						psf,
						Math.max(rtFit.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
								rtFit.getValueAsDouble(ResultsTableMt.SIGMAY, psf)));
				rtFit.setValue(
						ResultsTableMt.MIN_SIGMA,
						psf,
						Math.min(rtFit.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
								rtFit.getValueAsDouble(ResultsTableMt.SIGMAY, psf)));
				rtFit.setValue(ResultsTableMt.ELLIPTICITY, psf, rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, psf)
						/ rtFit.getValueAsDouble(ResultsTableMt.MIN_SIGMA, psf));
			}

			double[] rows = SRutil.getSortedIndexes(rtFit.getColumnAsDoubles(ResultsTableMt.MAX_SIGMA));

			for (int row = 0; row < rows.length; row++) {
				SRutil.addRow(rtFit, rtFitOrdered, (int) rows[row]);
				rtFitOrdered.addValue(ResultsTableMt.ORIGINAL_FIT, rows[row]);
			}

			rtFit = rtFitOrdered;

		}

		// Ask the thresholds and rebuild accordingly and interactively
		{
			double widthThreshold = Prefs.get("VLS.widthThreshold", 1.0);
			double amplitudeThreshold = Prefs.get("VLS.amplitudeThreshold", 1);

			NonBlockingGenericDialog gui = new NonBlockingGenericDialog("Choose thresholds");

			gui.addMessage("Thresholds in width and amplitude:");
			gui.addMessage("(Use '-1' to ignore the corresponding threshold)");
			gui.addSlider("Max width threshold (nm)", 0.0, 1.5 * widthThreshold, widthThreshold);
			gui.addSlider("Min width threshold (nm)", 0.0, 1.5 * widthThreshold, 0.0);
			double[] minMax = SRutil.getMinMax(rtFit, i0_);
			gui.addSlider("Amplitude threshold (photons/" + IJ.micronSymbol + "m^2)", minMax[0], minMax[1],
					amplitudeThreshold);
			minMax = SRutil.getMinMax(rtFit, ResultsTableMt.ELLIPTICITY);
			gui.addSlider("Ellipticity threshold:", minMax[0], minMax[1], 1.7);

			gui.addMessage("Rebuild options:");
			gui.addChoice("Rebuilding mode", rebuildingMode, rebuildingMode[rebuilding.index]);
			gui.addNumericField("Width of the av. rebuilt PSF (nm):", averagePrecision, 1);
			gui.addMessage("Average precision of all fitted PSFs (nm):" + IJ.d2s(averagePrecision, 2));
			gui.addNumericField("Rebuilt scale :", 10, 0);
			gui.addCheckbox("Weighted image", true);
			gui.addCheckbox("Equalised image", true);
			// gui.addNumericField("Rolling window", 0, 0);

			gui.addDialogListener(this);

			gui.showDialog();

			if (gui.wasCanceled()) {
				rebuildingSettings.precision = (float) averagePrecision;
				rebuilding(2.03, -1, 1502, -1, Rebuilding.QuickMt, rebuildingSettings);
			} else if (gui.wasOKed()) {
				rebuilding(sigma0, sigma0Min, i0, e0, rebuilding, rebuildingSettings);
				Prefs.set("VLS.widthThreshold", sigma0);
				Prefs.set("VLS.amplitudeThreshold", i0);
			}
		}

		// Save threshold values as metadata and the rebuilt picture
		{
			String metaData = "Sigma threshold (nm): " + sigma0 + "\n" + "Sigma min threshold (nm): "
					+ sigma0Min + "\n" + "Amplitude threshold (photons/" + IJ.micronSymbol + "m^2): " + i0
					+ "\n" + "Ellipticity threshold: " + e0 + "\n"
					+ "Standard deviation of the rebuilt Gaussian blur: " + rebuildingSettings.precision
					+ "\n" + "Physical size of a pixel (nm): " + param.param.pixelSize;
			SRutil.addMetaData(impRebuilt, metaData);
			SRutil.saveTiff(impRebuilt, param.fileDirName + File.separator
					+ "Quick_Super-Resolved_Filtered.tif", false);
		}

		// Save RtFitFiltered
		{
			ResultsTableMt rtFinal = new ResultsTableMt();

			int row = 0;
			while (row < rtFit.getCounter()
					&& (sigma0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0)) {
				IJ.showProgress(row, rtFit.getCounter());

				if (sigma0Min != -1 && rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0Min)
					;
				else if ((i0 == -1 || rtFit.getValueAsDouble(i0_, row) > i0)
						&& (e0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.ELLIPTICITY, row) < e0)) {
					SRutil.addRow(rtFit, rtFinal, row);
				}
				row++;
			}

			try {
				rtFinal.saveAs(param.fileDirName + File.separator + "FilteredResults.txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	@Override
	public boolean dialogItemChanged(GenericDialog gui, AWTEvent e) {

		sigma0 = gui.getNextNumber();
		sigma0Min = gui.getNextNumber();
		if (sigma0Min > sigma0)
			sigma0Min = sigma0;
		i0 = gui.getNextNumber();
		e0 = gui.getNextNumber();
		rebuilding = Rebuilding.values()[gui.getNextChoiceIndex()];
		rebuildingSettings.setResultsImage(rebuilding.code);
		rebuildingSettings.precision = (float) gui.getNextNumber();
		param.param.rebuiltScale = (int) gui.getNextNumber();
		rebuildingSettings.imageScale = param.param.rebuiltScale;
		rebuildingSettings.weightedImage = gui.getNextBoolean();
		rebuildingSettings.equalisedImage = gui.getNextBoolean();
		// rebuildingSettings.imageRollingWindow = (int)gui.getNextNumber();

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
							double mySigma0 = sigma0;
							double mySigma0Min = sigma0Min;
							double myI0 = i0;
							// Do something with parameters
							updateParameterPlot(mySigma0, mySigma0Min, myI0);
							// Check if the parameters have changed again
							parametersChanged = (mySigma0 != sigma0) || (mySigma0Min != sigma0Min)
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
							double mySigma0 = sigma0;
							double mySigma0Min = sigma0Min;
							double myI0 = i0;
							double mye0 = e0;
							Rebuilding myRebuilding = rebuilding;
							int myRebuiltScale = param.param.rebuiltScale;
							ResultsSettings myRebuildingSettings = new ResultsSettings();
							myRebuildingSettings.setResultsImage(rebuilding.code);
							myRebuildingSettings.imageScale = myRebuiltScale;
							myRebuildingSettings.weightedImage = rebuildingSettings.weightedImage;
							myRebuildingSettings.equalisedImage = rebuildingSettings.equalisedImage;
							myRebuildingSettings.imageRollingWindow = rebuildingSettings.imageRollingWindow;
							myRebuildingSettings.precision = rebuildingSettings.precision;
							// Do something with parameters
							rebuilding(mySigma0, mySigma0Min, myI0, mye0, myRebuilding, myRebuildingSettings);
							// Check if the parameters have changed again
							parametersChanged = (mySigma0 != sigma0)
									|| (mySigma0Min != sigma0Min)
									|| (myI0 != i0)
									|| (mye0 != e0)
									|| (myRebuilding != rebuilding)
									|| (myRebuiltScale != param.param.rebuiltScale)
									|| (myRebuildingSettings.weightedImage != rebuildingSettings.weightedImage)
									|| (myRebuildingSettings.equalisedImage != rebuildingSettings.equalisedImage)
									|| (myRebuildingSettings.imageRollingWindow != rebuildingSettings.imageRollingWindow)
									|| (myRebuildingSettings.precision != rebuildingSettings.precision);
						}
					} finally {
						// Ensure the running flag is reset
						imageLock = false;
					}
				}
			}).start();
		}
	}

	private void updateParameterPlot(double sigma0, double sigma0Min, double i0) {
		if (parameterPlot == null)
			getParameterPlot();
		parameterPlot.show();
		
		Overlay ol = paramPlots.getScaleBarOverlay(parameterPlot.getWidth(), parameterPlot.getHeight());
		
		int yMinOverlayTemp = Math.max(1, parameterPlot.getHeight() - (int) (sigma0 * paramPlots.scaleSigma));
		if (sigma0 > 0) {
			ImageRoi imgRoi1 = new ImageRoi(0, 0, new ShortProcessor(parameterPlot.getWidth(), yMinOverlayTemp));
			imgRoi1.setOpacity(0.5);
			imgRoi1.setStrokeColor(Color.BLACK);
			imgRoi1.setStrokeWidth(2);
			imgRoi1.setFillColor(Color.BLACK);
			ol.add(imgRoi1);
		}

		int yMaxOverlayTemp = parameterPlot.getHeight() - Math.min(parameterPlot.getHeight(), (int) (sigma0Min * paramPlots.scaleSigma));
		if (sigma0Min > 0) {
			ImageRoi imgRoi3 = new ImageRoi(0, yMaxOverlayTemp,
					new ShortProcessor(parameterPlot.getWidth(), parameterPlot.getHeight() - yMaxOverlayTemp));
			imgRoi3.setOpacity(0.5);
			imgRoi3.setStrokeColor(Color.BLACK);
			imgRoi3.setStrokeWidth(2);
			imgRoi3.setFillColor(Color.BLACK);
			ol.add(imgRoi3);
		}

		if (i0 > 0) {
			ImageRoi imgRoi2 = new ImageRoi(0, yMinOverlayTemp, new ShortProcessor(Math.max(1,
					Math.min(parameterPlot.getWidth(), (int) (i0 * paramPlots.scaleI0))),
					yMaxOverlayTemp - yMinOverlayTemp));
			imgRoi2.setOpacity(0.5);
			imgRoi2.setStrokeColor(Color.BLACK);
			imgRoi2.setStrokeWidth(2);
			imgRoi2.setFillColor(Color.BLACK);
			ol.add(imgRoi2);
		}

		parameterPlot.setOverlay(ol);
		parameterPlot.setHideOverlay(false);
		parameterPlot.updateAndDraw();
		parameterPlot.show();
	}

	final static int widthPlot = 251;//401;
	final static int heightPlot = 251;//501;

	private ImagePlus getParameterPlot() {
		paramPlots = new Bead_Stack_Analysis();
		paramPlots.newParameters(param.param, rtFit);

		// Get the SigmaThreshold
		ByteProcessor sigmaThresh = new ByteProcessor(widthPlot, heightPlot);

		ResultsTableMt rt = rtFit;
		paramPlots.scaleSigma = Math.min(
				100,
				(double) sigmaThresh.getHeight()
						/ (SRutil.getMean(rt, ResultsTableMt.SIGMAX) + 2.0*SRutil.getIQR(rt,
								ResultsTableMt.SIGMAX)));// 100.0D;
		paramPlots.scaleI0 = Math.min(1,
				(double) sigmaThresh.getWidth() / (SRutil.getMean(rt, i0_) + 5.0*SRutil.getIQR(rt, i0_)));
		// IJ.log(""+scaleI0);

		for (int psf = 0; psf < rt.getCounter(); psf++) {
			if (rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) < (double) (sigmaThresh.getHeight() - 1)
					/ paramPlots.scaleSigma
					&& rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) < (double) (sigmaThresh.getHeight() - 1)
							/ paramPlots.scaleSigma
					&& rt.getValueAsDouble(i0_, psf) < (double) (sigmaThresh.getWidth() - 1)
							/ paramPlots.scaleI0
					&& Math.min(
							Math.min(rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf),
									rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf)),
							rt.getValueAsDouble(i0_, psf)) > 0) {
				sigmaThresh
						.set((int) (rt.getValueAsDouble(i0_, psf) * paramPlots.scaleI0 + 0.5),
								(int) (paramPlots.scaleSigma
										* rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5),
								sigmaThresh.get(
										(int) (rt.getValueAsDouble(i0_, psf) * paramPlots.scaleI0 + 0.5),
										(int) (paramPlots.scaleSigma
												* rt.getValueAsDouble(ResultsTableMt.SIGMAX, psf) + 0.5)) + 1);
				sigmaThresh
						.set((int) (rt.getValueAsDouble(i0_, psf) * paramPlots.scaleI0 + 0.5),
								(int) (paramPlots.scaleSigma
										* rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5),
								sigmaThresh.get(
										(int) (rt.getValueAsDouble(i0_, psf) * paramPlots.scaleI0 + 0.5),
										(int) (paramPlots.scaleSigma
												* rt.getValueAsDouble(ResultsTableMt.SIGMAY, psf) + 0.5)) + 1);
			}
		}
		sigmaThresh.flipVertical();

		ByteProcessor blank = new ByteProcessor(widthPlot, heightPlot);
		blank.setRoi(new Rectangle(0, 0, blank.getWidth(), blank.getHeight()));
		blank.setColor(Color.WHITE);
		blank.fill();
		blank.setMinAndMax(0, 255);

		ImageCalculator ic = new ImageCalculator();
		ImageProcessor sigmaThreshFocalG = ic.run("Substract create", new ImagePlus("", blank),
				new ImagePlus("", sigmaThresh)).getProcessor();
		int min = SRutil.getMin(sigmaThreshFocalG.getIntArray());
		sigmaThreshFocalG.setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		ImagePlus[] temp = { new ImagePlus("", sigmaThreshFocalG), new ImagePlus("", blank),
				new ImagePlus("", sigmaThreshFocalG) };
		ImagePlus parameterPlot = RGBStackMerge.mergeChannels(temp, true);
		ImageConverter icon = new ImageConverter(parameterPlot);
		ImageConverter.setDoScaling(false);
		icon.convertRGBStackToRGB();
		parameterPlot.getProcessor().setMinAndMax(min + (int) (0.33D * (255 - min)), 255);
		parameterPlot.setOverlay(paramPlots.getScaleBarOverlay(sigmaThresh.getWidth(),
				sigmaThresh.getHeight()));
		parameterPlot.setHideOverlay(false);
		parameterPlot.updateAndDraw();
		parameterPlot.show();
		parameterPlot.getWindow().setLocation(536, 355);
		// new Zoom().run("in");
		SRutil.saveTiff(parameterPlot, param.fileDirName + File.separator + "4_Amplitude-Vs-Stdev.tif", false);
		this.parameterPlot = parameterPlot;

		return parameterPlot;
	}

	private void rebuilding(double sigma0, double sigma0Min, double i0, double e0, Rebuilding rebuilding,
			ResultsSettings rebuildingSettings) {
		IJ.showStatus("Rebuilding the new picture: Sigma0 = " + IJ.d2s(sigma0, 0) + " nm / A0 = "
				+ IJ.d2s(i0, 0) + " photons/" + IJ.micronSymbol + "m^2 / e0 = " + IJ.d2s(e0, 2));

		ResultsTableMt rtSelected = new ResultsTableMt();

		if (rebuilding == Rebuilding.QuickMt) {

			// Build the ImageProcessor SRRebuilt
			double[] xMinMax = SRutil.getMinMax(rtFit, ResultsTableMt.X0);
			double[] yMinMax = SRutil.getMinMax(rtFit, ResultsTableMt.Y0);
			xmin = (int) xMinMax[0] - 5;
			ymin = (int) yMinMax[0] - 5;
			srRebuilt = new ShortProcessor((((int) xMinMax[1] + 5) - xmin)
					* (int) rebuildingSettings.imageScale, (((int) yMinMax[1] + 5) - ymin)
					* (int) rebuildingSettings.imageScale);
			if (WindowManager.getImage(" SuperRes") == null) {
				impRebuilt = new ImagePlus(" SuperRes", srRebuilt);
			} else {
				impRebuilt = WindowManager.getImage(" SuperRes");
				impRebuilt.setProcessor(srRebuilt);
			}

			WindowManager.setTempCurrentImage(impRebuilt);
			LutLoader lut = new LutLoader();
			lut.run("fire");
			WindowManager.setTempCurrentImage(null);

			impRebuilt.show();

			srRebuilt.setColor(Color.BLACK);
			srRebuilt.fill(new Roi(new Rectangle(0, 0, srRebuilt.getWidth(), srRebuilt.getHeight())));
			short[] pixels = (short[]) srRebuilt.getPixels();

			int row = 0;
			int x;
			int y;
			while (row < rtFit.getCounter()
					&& (sigma0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0)) {
				IJ.showProgress(row, rtFit.getCounter());

				if (sigma0Min != -1 && rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0Min)
					;
				else if ((i0 == -1 || rtFit.getValueAsDouble(i0_, row) > i0)
						&& (e0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.ELLIPTICITY, row) < e0)) {
					x = (int) ((rtFit.getValueAsDouble(ResultsTableMt.X0, row) - xmin) * (double) rebuildingSettings.imageScale);
					y = (int) ((rtFit.getValueAsDouble(ResultsTableMt.Y0, row) - ymin) * (double) rebuildingSettings.imageScale);
					int index = y * srRebuilt.getWidth() + x;
					pixels[index]++;
					srRebuilt.set(x, y, srRebuilt.get(x, y) + 1);
					rtSelected.incrementCounter();
					rtSelected.addValue(ResultsTableMt.SIGMArebuilt,
							rtFit.getValueAsDouble(ResultsTableMt.SIGMArebuilt, row));
					// No need to change back the units of SigmaX and SigmaY (from nm to pix)
				}
				row++;
			}

			if (rebuildingSettings.precision > 0) {
				srRebuilt.multiply(65535.0 / SRutil.getMax(srRebuilt.getIntArray()));
				new GaussianBlur().blurGaussian(srRebuilt, rebuildingSettings.precision
						/ param.param.pixelSize * rebuildingSettings.imageScale, rebuildingSettings.precision
						/ param.param.pixelSize * rebuildingSettings.imageScale, 0.1);
			}
			srRebuilt.setAutoThreshold(ImageProcessor.ISODATA2, ImageProcessor.NO_LUT_UPDATE);

			if (param.param.scaleBar) {
				int ScaleBarSize = (int) Math.round(1000.0 / param.param.pixelSize
						* rebuildingSettings.imageScale);
				if (ScaleBarSize + 4 * rebuildingSettings.imageScale < srRebuilt.getWidth()) {
					srRebuilt.fill(new Roi(new Rectangle(2 * (int) rebuildingSettings.imageScale,
							(param.param.roi.height - 2) * (int) rebuildingSettings.imageScale, ScaleBarSize,
							(int) rebuildingSettings.imageScale)));
				}
			}

		} else {

			int row = 0;
			while (row < rtFit.getCounter()
					&& (sigma0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0)) {
				IJ.showProgress(row, rtFit.getCounter());

				if (sigma0Min != -1 && rtFit.getValueAsDouble(ResultsTableMt.MAX_SIGMA, row) < sigma0Min)
					;
				else if ((i0 == -1 || rtFit.getValueAsDouble(i0_, row) > i0)
						&& (e0 == -1 || rtFit.getValueAsDouble(ResultsTableMt.ELLIPTICITY, row) < e0)) {
					SRutil.addRow(rtFit, rtSelected, row);
					
					// Change back from nm to pix for SigmaX and SigmaY for correct rebuilding (The amplitude can stay as is as it is just porportional)
					rtSelected.setValue(ResultsTableMt.SIGMAX, rtSelected.getCounter()-1, rtSelected.getValueAsDouble(ResultsTableMt.SIGMAX, rtSelected.getCounter()-1)/param.param.pixelSize);
					rtSelected.setValue(ResultsTableMt.SIGMAY, rtSelected.getCounter()-1, rtSelected.getValueAsDouble(ResultsTableMt.SIGMAY, rtSelected.getCounter()-1)/param.param.pixelSize);

				}
				
				row++;
			}

			impRebuilt = SRutil.rebuildingPeakFit(rtSelected, "", param.param, rebuildingSettings);
		}

		IJ.showStatus("Rebuilding picture done!");

		// SRutil.getPrecisionHistogram does not re-calculate the SigmaRebuilt, so no need to check the units of SigmaX and SigmaY
		double AverageSelectedPrecision = SRutil.getPrecisionHistogram(rtSelected, param.param, false, false);
		Overlay ov = new Overlay();
		TextRoi text = new TextRoi(10, 10, "Average precision of all selected PSFs (nm): "
				+ IJ.d2s(AverageSelectedPrecision, 2));
		text.setStrokeColor(Color.WHITE);
		text.setStrokeWidth(50);
		text.setNonScalable(true);
		ov.addElement(text);

		// impRebuilt.setProcessor(srRebuilt);
		impRebuilt.setOverlay(ov);
		impRebuilt.setHideOverlay(false);
		impRebuilt.updateAndDraw();
		impRebuilt.show();
	}

}
