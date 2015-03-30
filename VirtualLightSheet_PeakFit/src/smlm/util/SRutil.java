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

package smlm.util;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.settings.ResultsSettings;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.measure.ResultsTable;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.ImageCalculator;
import ij.plugin.LutLoader;
import ij.plugin.filter.GaussianBlur;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.Color;
import java.awt.Font;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import smlm.Params;
import smlm.fitting.FittingClassical;


public class SRutil {

	public final static Color lightGrey = new Color(200, 200, 200);
	public final static Color grey = new Color(150, 150, 150);
	public final static Color royalBlue = new Color(0, 0, 128);
	public final static Color yellow = new Color(255, 165, 0);

	// OPERATE FILES

	public static String getFileName(String localRootDirName, String ending) {

		File folder = new File(localRootDirName);
		File[] listOfFiles = folder.listFiles();

		// Find the image which name contains such an Ending
		int i = 0;
		for (i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].getName().indexOf(ending) >= 0)
				break;
		}

		if (i == listOfFiles.length)
			return null;
		else
			return listOfFiles[i].getName();

	}

	public static String[] getSubDirNames(String dir) {

		File folder = new File(dir);
		File[] listOfFiles = folder.listFiles();

		// Count the number of sub-directories
		int count = 0;
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isDirectory())
				count++;
		}

		// Fill SubDirNames
		String[] subDirNames = new String[count];
		count = 0;
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isDirectory()) {
				subDirNames[count] = listOfFiles[i].getAbsolutePath();
				count++;
			}
		}

		return subDirNames;
	}

	public static ImagePlus openImage(String localRootDirName, String ending, boolean virtual,
			boolean firstFrameOnly) {

		IJ.freeMemory();

		File folder = new File(localRootDirName);
		File[] listOfFiles = folder.listFiles();

		// Find the image which name contains such an Ending
		int i = 0;
		for (i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].getName().indexOf(ending) >= 0)
				break;
		}

		if (i == listOfFiles.length)
			return null;
		else {
			Opener opener = new Opener();
			ImagePlus temp = opener.openImage(listOfFiles[i].getAbsolutePath(), 1);
			FileInfo info = temp.getOriginalFileInfo();
			long memorySize = new File(listOfFiles[i].getAbsolutePath()).length();// info.getBytesPerPixel()*info.width*info.height*info.nImages;
			temp.close();
			if (firstFrameOnly)
				return temp;
			temp.flush();

			if (!virtual)
				virtual = (IJ.maxMemory() - IJ.currentMemory() < 1.3D * memorySize);
			// IJ.log("Max mem: "+IJ.maxMemory()+" Used mem: "+IJ.currentMemory());
			// IJ.log("Size: "+MemorySize);
			if (virtual) {
				IJ.log("Open virtual tiff (" + IJ.d2s((IJ.maxMemory() - IJ.currentMemory()) / 100000.0D, 2)
						+ " Gb): " + listOfFiles[i].getName());
				ImagePlus imp = new ImagePlus(listOfFiles[i].getName(), (new FileInfoVirtualStack(info)));
				imp.setFileInfo(info);
				imp.close();
				return imp;
			} else {
				IJ.log("Open tiff: " + listOfFiles[i].getAbsolutePath());
				return opener.openImage(listOfFiles[i].getAbsolutePath());
			}
		}
	}

	public static ImagePlus openImage(String localRootDirName, String ending, boolean virtual) {
		return openImage(localRootDirName, ending, virtual, false);
	}

	public static ImagePlus openImage(String localRootDirName, String ending) {
		return openImage(localRootDirName, ending, false);
	}

	public static String[] getAFile(String title, String initialRoot, String initialFile) {
		if (title == null || title == "")
			title = "Choose an analysed table:";
		if (initialRoot == null || initialRoot == "")
			initialRoot = "E:\\Data";
		OpenDialog dialog = new OpenDialog(title, initialRoot, initialFile);
		
		if(!new File(dialog.getDirectory()+File.separator+dialog.getFileName()).isFile())
			return null;
		
		String[] retour = { dialog.getDirectory(), dialog.getFileName() };
		return retour;
	}

	public static String[] getAFile() {
		return getAFile(null, null, "");
	}

	public static String[][] getFiles(String title, String initialRoot, String initialFile) {
		if (title == null || title == "")
			title = "Choose files";
		if (initialRoot == null || initialRoot == "")
			initialRoot = "E:\\Data";

		JFileChooser chooser = new JFileChooser();
		chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		FileNameExtensionFilter filter = new FileNameExtensionFilter("Tiff files", "tif", "tiff");
		chooser.setFileFilter(filter);
		chooser.setDialogTitle(title);
		chooser.setCurrentDirectory(new File(initialRoot));
		chooser.setMultiSelectionEnabled(true);

		int returnVal = chooser.showOpenDialog(null);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String[][] files = new String[chooser.getSelectedFiles().length][2];
			for (int i = 0; i < files.length; i++) {
				files[i][0] = chooser.getSelectedFiles()[i].getParent();
				files[i][1] = chooser.getSelectedFiles()[i].getName();
			}
			return files;
		} else
			return null;
	}

	public static void writeToFile(String h, String rootDirName) {
		try {
			FileWriter fw = new FileWriter(rootDirName, true);
			BufferedWriter bw = new BufferedWriter(fw);

			String stringFile = h;
			bw.write(stringFile);
			bw.newLine();
			bw.close();
		} catch (Exception e) {
			System.out.println("Exception: " + e);
		}
	}

	public static void saveTiff(ImagePlus imp, String fullName, boolean flush) {
		ij.io.FileSaver fs = new ij.io.FileSaver(imp);
		if (imp.getImageStackSize() > 1)
			fs.saveAsTiffStack(fullName);
		else
			fs.saveAsTiff(fullName);
		if (flush)
			imp.flush();
	}

	public static void addMetaData(ImagePlus imp, String metaData) {
		Object prop = imp.getProperty("Info");
		String oldInfo = (String) prop;

		String newInfo = metaData + "\n" + oldInfo;

		imp.setProperty("Info", newInfo);
	}

	public static Color getGradientColor(Color ini, Color fin, int numberOfSteps, int step) {
		if (numberOfSteps == 0)
			return ini;
		if (step % numberOfSteps == 0)
			return ini;
		if (step % numberOfSteps == numberOfSteps - 1)
			return fin;
		else {
			float[] hsbini = Color.RGBtoHSB(ini.getRed(), ini.getGreen(), ini.getBlue(), null);
			float[] hsbfin = Color.RGBtoHSB(fin.getRed(), fin.getGreen(), fin.getBlue(), null);
			float h = ((numberOfSteps - step % numberOfSteps) * hsbini[0] + (step % numberOfSteps)
					* hsbfin[0])
					/ numberOfSteps;
			float s = ((numberOfSteps - step % numberOfSteps) * hsbini[1] + (step % numberOfSteps)
					* hsbfin[1])
					/ numberOfSteps;
			float b = ((numberOfSteps - step % numberOfSteps) * hsbini[2] + (step % numberOfSteps)
					* hsbfin[2])
					/ numberOfSteps;
			return Color.getHSBColor(h, s, b);
		}
	}

	public static void saveVector(double[] vector, String path) {
		ResultsTableMt rt = new ResultsTableMt();
		for (int i = 0; i < vector.length; i++) {
			rt.incrementCounter();
			rt.addValue(0, vector[i]);
		}

		try {
			rt.saveAs(path);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static double[] openVector(String path) {
		ResultsTableMt vector = ResultsTableMt.open2(path);
		double[] retour = new double[vector.getCounter()];
		for (int row = 0; row < vector.getCounter(); row++) {
			retour[row] = vector.getValue(vector.getColumnHeading(0), row);
		}

		return retour;
	}

	// OPERATION ON VECTORS OR MATRICES (Average, Stdev...)

	public static double stdev(double[] vector, double mean, boolean donTCountZeros) {
		double retour = 0;
		int total = 0;
		if (vector == null)
			return 0;
		else {
			for (int j = 0; j < vector.length; j++) {
				if (donTCountZeros && vector[j] == 0) {
				} else {
					retour += java.lang.Math.pow((vector[j] - mean), 2);
					total++;
				}
			}
			if (total == 0)
				return 0;
			else
				return java.lang.Math.sqrt(retour / total);
		}
	}

	public static double stdev(double[][] matrix, double mean, int row, boolean donTCountZeros) {
		double retour = 0;
		int total = 0;
		if (matrix == null)
			return 0;
		else {
			for (int column = 0; column < matrix.length; column++) {
				if (matrix[column] != null && row < matrix[column].length) {
					if (donTCountZeros && matrix[column][row] == 0) {
					} else {
						retour += java.lang.Math.pow((matrix[column][row] - mean), 2);
						total++;
					}
				}
			}
			if (total == 0)
				return 0;
			else
				return java.lang.Math.sqrt(retour / total);
		}
	}

	public static double stdev(int[][] vector, double mean) {
		double retour = 0;
		for (int i = 0; i < vector.length; i++) {
			for (int j = 0; j < vector[0].length; j++) {
				retour += java.lang.Math.pow((vector[i][j] - mean), 2);
			}
		}
		return java.lang.Math.sqrt(retour / (vector.length * vector[0].length));
	}

	public static double averaging(int[][] vector) {
		double retour = 0;
		for (int i = 0; i < vector.length; i++) {
			for (int j = 0; j < vector[0].length; j++) {
				retour += vector[i][j];
			}
		}
		return retour / (vector.length * vector[0].length);
	}

	public static double averaging(double[] vector, boolean donTCountZeros) {
		double retour = 0;
		int total = 0;
		if (vector == null)
			return 0;
		else {
			for (int i = 0; i < vector.length; i++) {
				if (donTCountZeros && vector[i] == 0) {
				} else {
					retour += vector[i];
					total++;
				}
			}
			if (total == 0)
				return 0;
			else
				return retour / total;
		}
	}

	public static double averaging(double[][] matrix, int row, boolean donTCountZeros) {
		double retour = 0;
		int total = 0;
		if (matrix == null)
			return 0;
		else {
			for (int column = 0; column < matrix.length; column++) {
				if (matrix[column] != null && row < matrix[column].length) {
					if (donTCountZeros && matrix[column][row] == 0) {
					} else {
						retour += matrix[column][row];
						total++;
					}
				}
			}
			if (total == 0)
				return 0;
			else
				return retour / total;
		}
	}

	public static double[] runningAverage(double[] vector, int halfWindow) {
		double[] retour = new double[vector.length];
		double[] sum = new double[2 * halfWindow + 1];
		for (int i = 0; i < vector.length; i++) {
			for (int j = 0; j < sum.length; j++)
				sum[j] += vector[i];
			if (i < sum.length)
				retour[i] = sum[i % sum.length] / (double) (i + 1);
			else
				retour[i] = sum[i % sum.length] / (double) sum.length;
			sum[i % sum.length] = 0;
		}
		return retour;
	}

	public static double[] derivative(double[] vector) {
		double[] retour = new double[vector.length];
		for (int i = 1; i < vector.length - 1; i++) {
			retour[i] = (vector[i + 1] - vector[i - 1]) / 2.0D;
		}
		retour[0] = vector[1] - vector[0];
		retour[vector.length - 1] = vector[vector.length - 1] - vector[vector.length - 2];
		return retour;
	}

	public static double[] derivativeStack(ImageStack stack, int window) {
		double[] retour = new double[stack.getSize()];
		for (int i = window; i < retour.length - window; i++) {
			ImagePlus impTemp = new ImageCalculator().run("Subtract create 32-bit",
					new ImagePlus("", stack.getProcessor(i + window + 1)),
					new ImagePlus("", stack.getProcessor(i - window + 1)));
			impTemp.getProcessor().abs();
			impTemp.getProcessor().multiply(0.5D);
			retour[i] = impTemp.getStatistics().mean;
		}
		return retour;
	}

	public static double[][] derivativeStack(ImageStack stack) {
		double[][] retour = new double[stack.getSize() - 2][stack.getSize()];
		for (int window = 1; window < retour.length; window++)
			retour[window - 1] = derivativeStack(stack, window);

		double maxIntensityPerPixelFrame = 0;
		for (int i = 0; i < stack.getSize(); i++)
			maxIntensityPerPixelFrame = Math.max(maxIntensityPerPixelFrame, stack.getProcessor(i + 1)
					.getStatistics().mean);
		for (int i = 0; i < retour.length; i++) {
			for (int j = 0; j < retour[i].length; j++)
				retour[i][j] /= maxIntensityPerPixelFrame / 1000.0;
		}

		return retour;
	}

	// OPERATION ON Rt (TableFit)

	public static double translateCountsToPhotons(double i0, Params param, boolean removeOffset) {
		if (!param.photonCountMode) {
			return translateCountsToPhotons(i0, param.emGain
					/ Params.ConversionGain[param.camera][param.gain], removeOffset ? param.offset : 0);
		}
		return i0;
	}

	public static double translatePhotonsToCounts(double i0, Params param, boolean removeOffset) {
		return translatePhotonsToCounts(i0, param.emGain / Params.ConversionGain[param.camera][param.gain],
				removeOffset ? param.offset : 0);
	}

	public static double translatePhotonsToCounts(double i0, double aduPerPhoton, double offset) {
		return i0 * aduPerPhoton + offset;
	}

	public static double translateCountsToPhotons(double i0, double aduPerPhoton, double offset) {
		return Math.max((i0 - offset) / aduPerPhoton, 0);
	}

	public static ResultsTableMt translateCountsToPhotons(ResultsTableMt rt, String i0, Params param,
			boolean removeOffset) {
		if (rt == null)
			return null;

		double[] values = getColumnFromRt(rt, i0);
		if (values != null) {
			int index = rt.getColumnIndex(i0);
			for (int row = 0; row < rt.getCounter(); row++) {
				rt.setValue(index, row, translateCountsToPhotons(values[row], param, removeOffset));
			}
		}
		return rt;
	}

	public static ResultsTableMt translateCountsToPhotons(ResultsTableMt rt, int i0, Params param,
			boolean removeOffset) {
		if (rt == null)
			return null;

		double[] values = rt.getColumnAsDoubles(i0);
		for (int row = 0; row < rt.getCounter(); row++) {
			rt.setValue(i0, row, translateCountsToPhotons(values[row], param, removeOffset));
		}
		return rt;
	}

	public static ResultsTableMt translateCountsToPhotons(ResultsTableMt rt, String i0, double aduPerPhoton,
			double offset) {
		if (rt == null)
			return null;

		double[] values = getColumnFromRt(rt, i0);
		if (values != null) {
			int index = rt.getColumnIndex(i0);
			for (int row = 0; row < rt.getCounter(); row++) {
				rt.setValue(index, row, translateCountsToPhotons(values[row], aduPerPhoton, offset));
			}
		}
		return rt;
	}

	public static ResultsTableMt translateCountsToPhotons(ResultsTableMt rt, int i0, double aduPerPhoton,
			double offset) {
		if (rt == null)
			return null;

		double[] values = rt.getColumnAsDoubles(i0);
		for (int row = 0; row < rt.getCounter(); row++) {
			rt.setValue(i0, row, translateCountsToPhotons(values[row], aduPerPhoton, offset));
		}
		return rt;
	}

	public static double sigmaRebuilt(Params param, ResultsTableMt rt, int row) {
		return sigmaRebuilt(param, rt.getValueAsDouble(ResultsTableMt.SIGMAX, row),
				rt.getValueAsDouble(ResultsTableMt.SIGMAY, row), rt.getValueAsDouble(ResultsTableMt.I0, row),
				rt.getValueAsDouble(ResultsTableMt.NOISE, row));
	}

	public static double sigmaRebuilt(Params param, double sigmaX, double sigmaY, double i0, double noise) {
		return sigmaRebuilt(param, sigmaX, sigmaY, i0, noise, 2);
	}

	public static double sigmaRebuilt(Params param, double sigmaX, double sigmaY, double i0, double noise,
			int dimension) {

		if (!param.photonCountMode) {
			i0 = translateCountsToPhotons(i0, param, false);
			noise = translateCountsToPhotons(noise, param, true);
		}
		if (param.fitting == Params.Fitting.CountingFit
				|| param.rebuilding == Params.Rebuilding.IdenticalCirclePSF
				|| param.rebuilding == Params.Rebuilding.IdenticalGaussianPSF)
			return param.fixedPSFSizeRebuilt / param.pixelSize;
		if (dimension == 2) {
			double n = i0;
			if (param.theoreticalRebuilt)
				return Math.pow(1.83
						* (sigmaX * sigmaY + 1.0 / 12.0)
						/ n
						* (16.0 / 9.0 + 8.0 * Math.PI * (sigmaX * sigmaY + 1.0 / 12.0) * Math.pow(noise, 2)
								/ n), 0.5);
			else
				return 1.6142D * Math.pow(i0, -0.401);// Cf. MP069
			// Stability.xlsx
			// Math.pow(4200000D/(Math.pow(param.PixelSize,2)*I0),0.5);
		} else if (dimension == 1) {
			double n = i0;
			return Math.pow(
					1.83
							* (sigmaX * sigmaY + 1.0 / 12.0)
							/ n
							* (16.0 / 9.0 + 4.0 * Math.sqrt(Math.PI * (sigmaX * sigmaY + 1.0 / 12.0)) * noise
									/ n), 0.5);
		} else
			return 0;

	}

	public static double sigmaZRebuilt(Params param, ResultsTableMt rt, int row) {
		return sigmaRebuilt(param, rt.getValueAsDouble(ResultsTableMt.SIGMAZ, row),
				rt.getValueAsDouble(ResultsTableMt.SIGMAZ, row), rt.getValueAsDouble(ResultsTableMt.IZ, row),
				rt.getValueAsDouble(ResultsTableMt.NOISE, row), 1);
	}

	public static double getAmplitude(double i0, double sigmaX, double sigmaY) {
		return i0 / (2.0 * Math.PI * sigmaX * sigmaY);
	}

	public static double getAmplitude(ResultsTableMt rt, int row) {
		return getAmplitude(rt.getValueAsDouble(ResultsTableMt.I0, row),
				rt.getValueAsDouble(ResultsTableMt.SIGMAX, row),
				rt.getValueAsDouble(ResultsTableMt.SIGMAY, row));
	}

	public static ResultsTableMt refreshPostProcess(ResultsTableMt rt, ResultsTableMt tableMax,
			double[] backgroundPerFrame, Params param) {
		if (rt == null)
			return null;
		double[] vectorFrame = rt.getColumnAsDoubles(ResultsTableMt.FRAME), vectorX, vectorY;

		// Post-process
		if (tableMax != null && tableMax.getCounter() >= rt.getCounter()) {
			vectorX = tableMax.getColumnAsDoubles(ResultsTableMt.Xmax);
			vectorY = tableMax.getColumnAsDoubles(ResultsTableMt.Ymax);
		} else {
			vectorX = rt.getColumnAsDoubles(ResultsTableMt.X0);
			vectorY = rt.getColumnAsDoubles(ResultsTableMt.Y0);
		}
		if (vectorFrame != null && vectorX != null && vectorY != null) {
			for (int row = 0; row < rt.getCounter(); row++) {
				refreshPostProcessPerRow(rt, row, vectorX[row], vectorY[row],
						backgroundPerFrame[Math.max(0, (int) vectorFrame[row] - 1)], param);
			}
		}
		return rt;
	}

	public static ResultsTableMt refreshPostProcessPerRow(ResultsTableMt rt, int row, double Xmax,
			double Ymax, double BackgroundPerFrame, Params param) {

		// Post-process
		if (param.fitting != Params.Fitting.CountingFit) {
			if (param.i0Threshold > 0
					&& (rt.getValueAsDouble(ResultsTableMt.NOISE, row) < 100 && rt.getValueAsDouble(
							ResultsTableMt.I0, row) / rt.getValueAsDouble(ResultsTableMt.NOISE, row) < 10.0)
					|| (rt.getValueAsDouble(ResultsTableMt.NOISE, row) >= 100 && rt.getValueAsDouble(
							ResultsTableMt.I0, row) / rt.getValueAsDouble(ResultsTableMt.NOISE, row) < 2.0)) {
				rt.setValue(ResultsTableMt.IS_FITTED, row, -1);
			} else if ((rt.getValueAsDouble(ResultsTableMt.NOISE, row) < 100 && rt.getValueAsDouble(
					ResultsTableMt.I0, row)
					/ (2.0 * Math.PI * rt.getValueAsDouble(ResultsTableMt.SIGMAX, row) * rt.getValueAsDouble(
							ResultsTableMt.SIGMAY, row)) < 1.0 * BackgroundPerFrame)
					|| (rt.getValueAsDouble(ResultsTableMt.NOISE, row) >= 100 && rt.getValueAsDouble(
							ResultsTableMt.I0, row)
							/ (2.0 * Math.PI * rt.getValueAsDouble(ResultsTableMt.SIGMAX, row) * rt
									.getValueAsDouble(ResultsTableMt.SIGMAY, row)) < 1.0 * BackgroundPerFrame)) { // Centered
				// position higher than 1.0*Background
				rt.setValue(ResultsTableMt.IS_FITTED, row, -2);
			} else if (param.ellipticityThreshold > 0
					&& (rt.getValueAsDouble(ResultsTableMt.SIGMAX, row)
							/ rt.getValueAsDouble(ResultsTableMt.SIGMAY, row) > Math.max(
							param.ellipticityThreshold, 1.0 / param.ellipticityThreshold) || rt
							.getValueAsDouble(ResultsTableMt.SIGMAY, row)
							/ rt.getValueAsDouble(ResultsTableMt.SIGMAX, row) > Math.max(
							param.ellipticityThreshold, 1.0 / param.ellipticityThreshold))) {
				rt.setValue(ResultsTableMt.IS_FITTED, row, -3);
			} else if ((param.sigmaThreshold > 0 && rt.getValueAsDouble(ResultsTableMt.SIGMAX, row) > param.sigmaThreshold)
					|| rt.getValueAsDouble(ResultsTableMt.SIGMAX, row) <= 0
					|| (param.sigmaThreshold > 0 && rt.getValueAsDouble(ResultsTableMt.SIGMAY, row) > param.sigmaThreshold)
					|| rt.getValueAsDouble(ResultsTableMt.SIGMAY, row) <= 0) {
				rt.setValue(ResultsTableMt.IS_FITTED, row, -4);
			} else if (param.fitting != Params.Fitting.PeakFit
					&& Math.pow(rt.getValueAsDouble(ResultsTableMt.X0, row) - Xmax, 2)
							+ Math.pow(rt.getValueAsDouble(ResultsTableMt.Y0, row) - Ymax, 2) > Math.pow(2.0,
							2)) {
				rt.setValue(ResultsTableMt.IS_FITTED, row, -5);
			} else
				rt.setValue(ResultsTableMt.IS_FITTED, row, 1);
		}

		rt.setValue(ResultsTableMt.SIGMArebuilt, row, sigmaRebuilt(param, rt, row));
		if (rt.getValueAsDouble(ResultsTableMt.SIGMAZ, row) != 0)
			rt.setValue(ResultsTableMt.SIGMArebuiltZ, row, sigmaZRebuilt(param, rt, row));
		rt.setValue(ResultsTableMt.AMPLITUDE, row, getAmplitude(rt, row));

		return rt;
	}

	public static int getFrameNumberPerCycle(ResultsTableMt Rt) {
		int FrameNumberPerCycle = 0;
		for (int row = 0; row < Rt.getCounter(); row++) {
			if (!Rt.columnExists(ResultsTableMt.CYCLE) || Rt.getValueAsDouble(ResultsTableMt.CYCLE, row) != 1) {
				FrameNumberPerCycle = Math.max(FrameNumberPerCycle, (int) Rt.getValue("Frame", row));
			}
		}
		return FrameNumberPerCycle;
	}

	public static double[][] transpose(double[][] M) {
		return new Array2DRowRealMatrix(M).transpose().getData();
	}

	public static double[] getSortedIndexes(double[] vector) {
		double[] Rows = new double[vector.length];
		for (int row = 0; row < Rows.length; row++)
			Rows[row] = row;

		// Sort Rt following the "frame" column
		double[][] RtDouble = new double[2][vector.length];
		RtDouble[0] = Rows;
		RtDouble[1] = vector;

		RtDouble = transpose(RtDouble);
		class Comp implements Comparator<double[]> {
			@Override
			public int compare(double[] o1, double[] o2) {
				return ((Double) o1[1]).compareTo(o2[1]);
			}
		}
		Comp comp = new Comp();
		Arrays.sort(RtDouble, comp);
		return transpose(RtDouble)[0];
	}

	public static ResultsTableMt sortRt(ResultsTableMt rt, String Column) {
		if (rt.getColumnIndex(Column) == -1)
			return null;

		IJ.showStatus("Sorting results by " + Column);

		double[] rows = getSortedIndexes(rt.getColumnAsDoubles(rt.getColumnIndex(Column)));

		ResultsTableMt rt_Sorted = new ResultsTableMt();
		rt_Sorted.incrementCounter();
		rt_Sorted.addValue("OriginalFit", 0);
		rt_Sorted.deleteRow(0);
		int originalFit = rt_Sorted.getColumnIndex("OriginalFit");
		for (int row = 0; row < rows.length; row++) {
			addRow(rt, rt_Sorted, (int) rows[row]);
			rt_Sorted.addValue(originalFit, rows[row]);
		}

		return rt_Sorted;
	}

	public static double[] log10(double[] vector) {
		double[] retour = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			retour[i] = Math.log10(vector[i]);
		}
		return retour;
	}

	public static double getSum(double[] vector) {
		double sum = 0;
		for (int i = 0; i < vector.length; i++)
			sum += vector[i];
		return sum;
	}

	public static ResultsTable[] getGroups(ResultsTable[] Rts_Groups) {
		if (Rts_Groups.length <= 3)
			return null;

		ResultsTable[] Groups = new ResultsTable[Rts_Groups.length - 3];
		for (int groupi = 3; groupi < Rts_Groups.length; groupi++) {
			Groups[groupi - 3] = Rts_Groups[groupi];
		}
		return Groups;
	}

	public static int[] getIntensityTraceOverTime(ResultsTableMt Groups_group, int group, ImageStack vStack,
			Params param, boolean TraceAndSaveData) {
		return getIntensityTraceOverTime(Groups_group, group, vStack, param, TraceAndSaveData, null);
	}

	public static int[] getIntensityTraceOverTime(ResultsTableMt Groups_group, int group, ImageStack vStack,
			Params param, boolean TraceAndSaveData, ResultsTableMt tON_OFF) {
		int totalFrame = vStack.getSize();// (int)SRutil.getMax(Rt, "Frame");

		double[][] TraceFitted = new double[2][totalFrame];
		double[] x = new double[totalFrame];
		int totalIntensity = 0;
		for (int fit = 0; fit < Groups_group.getCounter(); fit++) {
			TraceFitted[0][(int) Groups_group.getValue("RealFrame", fit) - 1] = Groups_group
					.getValueAsDouble(ResultsTableMt.I0, fit);
			totalIntensity += Groups_group.getValueAsDouble(ResultsTableMt.I0, fit);
		}

		// Group everything within 50 frames
		ResultsTableMt TimeGroups = new ResultsTableMt();
		boolean Grouping = false;
		int lastFramePartOfGroup = 0;
		for (int frame = 0; frame < totalFrame; frame++) {
			if (frame - lastFramePartOfGroup >= 50)
				Grouping = false;
			if (!Grouping) {
				while (frame < totalFrame && TraceFitted[0][frame] == 0)
					frame++;
				if (frame < totalFrame) {
					Grouping = true;
					lastFramePartOfGroup = frame;
					TimeGroups.incrementCounter();
					TimeGroups.addValue("NumberOfPSFs", 1);
					TraceFitted[1][frame] = TimeGroups.getCounter();
				}
			} else {
				if (TraceFitted[0][frame] != 0) {
					lastFramePartOfGroup = frame;
					TimeGroups.addValue("NumberOfPSFs",
							TimeGroups.getValue("NumberOfPSFs", TimeGroups.getCounter() - 1) + 1);
					TraceFitted[1][frame] = TimeGroups.getCounter();
				}
			}
		}

		// Count groups having more than two PSFs (noise otherwise)
		int numberOfGroups = 0;
		for (int group_ = 0; group_ < TimeGroups.getCounter(); group_++) {
			if (TimeGroups.getValue("NumberOfPSFs", group_) > 2) {
				numberOfGroups++;
			} else {
				for (int frame = 0; frame < totalFrame; frame++) {
					if (TraceFitted[1][frame] == group_ + 1)
						TraceFitted[1][frame] = -group_ - 1;
				}
			}
		}

		// Fill tON_OFF
		if (tON_OFF != null) {
			tON_OFF.incrementCounter();
			tON_OFF.addValue("tON_Monomer", 0);
			tON_OFF.addValue("tOFF_Monomer", 0);
			tON_OFF.addValue("tON_Adjacent", 0);
			tON_OFF.addValue("tOFF_Adjacent", 0);
			tON_OFF.addValue("NbreOfPSFsPerMonomer", 0);

			for (int group_ = 0; group_ < TimeGroups.getCounter(); group_++) {
				if (TimeGroups.getValue("NumberOfPSFs", group_) > 2) {
					tON_OFF.setValue("NbreOfPSFsPerMonomer", tON_OFF.getCounter() - 1,
							TimeGroups.getValue("NumberOfPSFs", group_));
					tON_OFF.incrementCounter();
				}
			}

			boolean ON_Monomer = false;
			int previousFrameMonomer = 0;
			boolean ON_Adjacent = false;
			int previousFrameAdj = 0;
			int[] rowsON_OFF = new int[4];
			for (int frame = 0; frame < totalFrame; frame++) {

				if (ON_Monomer && TraceFitted[1][frame] <= 0) {
					if (tON_OFF.getCounter() <= rowsON_OFF[0])
						tON_OFF.incrementCounter();
					tON_OFF.setValue("tON_Monomer", rowsON_OFF[0], frame - previousFrameMonomer);
					previousFrameMonomer = frame;
					ON_Monomer = false;
					rowsON_OFF[0]++;
				} else if (!ON_Monomer && TraceFitted[1][frame] > 0) {
					if (tON_OFF.getCounter() <= rowsON_OFF[1])
						tON_OFF.incrementCounter();
					tON_OFF.setValue("tOFF_Monomer", rowsON_OFF[1], frame - previousFrameMonomer);
					previousFrameMonomer = frame;
					ON_Monomer = true;
					rowsON_OFF[1]++;
				}

				if (ON_Adjacent && TraceFitted[1][frame] == 0) {
					if (tON_OFF.getCounter() <= rowsON_OFF[2])
						tON_OFF.incrementCounter();
					tON_OFF.setValue("tON_Adjacent", rowsON_OFF[2], frame - previousFrameAdj);
					previousFrameAdj = frame;
					ON_Adjacent = false;
					rowsON_OFF[2]++;
				} else if (!ON_Adjacent && TraceFitted[1][frame] != 0) {
					if (tON_OFF.getCounter() <= rowsON_OFF[3])
						tON_OFF.incrementCounter();
					tON_OFF.setValue("tOFF_Adjacent", rowsON_OFF[3], frame - previousFrameAdj);
					previousFrameAdj = frame;
					ON_Adjacent = true;
					rowsON_OFF[3]++;
				}

			}

			// for(int row=0; row<tON_OFF.getCounter(); row++)
			// if(tON_OFF.getValue("NbreOfPSFsPerMonomer", row) != 0)
			// IJ.log("c " + tON_OFF.getValue("NbreOfPSFsPerMonomer", row));
		}

		if (TraceAndSaveData) {
			double[] TraceData = new double[totalFrame];
			ResultsTableMt OutputGroups = new ResultsTableMt();
			OutputGroups = ResultsTableMt.open2(param.fileDirName + File.separator + "OutPut_Groups.txt");

			int row = 0;
			while (row < OutputGroups.getCounter()
					&& OutputGroups.getValueAsDouble(ResultsTableMt.GROUP, row) < group + 1)
				row++;
			if (OutputGroups.getValueAsDouble(ResultsTableMt.GROUP, row) != group + 1)
				IJ.log("" + group + " // "
						+ (int) (OutputGroups.getValueAsDouble(ResultsTableMt.X0, row) + 0.5D) + " // "
						+ (int) (OutputGroups.getValueAsDouble(ResultsTableMt.Y0, row) + 0.5D));

			for (int frame = 0; frame < totalFrame; frame++) {
				if (!vStack.isVirtual()) {
					ImageProcessor ipTemp = vStack.getProcessor(frame + 1);
					ipTemp.setRoi(new Rectangle(
							(int) (OutputGroups.getValueAsDouble(ResultsTableMt.X0, row) + 0.5D) - 4,
							(int) (OutputGroups.getValueAsDouble(ResultsTableMt.Y0, row) + 0.5D) - 4, 9, 9));

					TraceData[frame] = ipTemp.crop().getStatistics().mean * 10.0D;
				}
				x[frame] = frame;
			}

			double[] minMaxIni;
			double MaxFin = SRutil.getMax(TraceFitted[0]);
			if (!vStack.isVirtual()) {
				minMaxIni = SRutil.getMinMax(TraceData);
				TraceData = new Array2DRowRealMatrix(TraceData).scalarAdd(-minMaxIni[0])
						.scalarMultiply(MaxFin / (minMaxIni[1] - minMaxIni[0])).getColumn(0);
			}

			Plot plot = new Plot("", "frame #  (number of groups: " + numberOfGroups + " )",
					"intensity of the fit (G) (photons) or of the data (R) (A.U.) (group #" + (group + 1)
							+ ")");
			plot.setSize(2048, 512);
			plot.setLimits(0, totalFrame - 1, 0, MaxFin);
			if (!vStack.isVirtual()) {
				plot.setLineWidth(1);
				plot.setColor(Color.RED);
				plot.addPoints(x, TraceData, Plot.LINE);
				plot.draw();
			}
			plot.setLineWidth(2);
			plot.setColor(Color.GREEN);
			plot.addPoints(x, TraceFitted[0], Plot.LINE);
			plot.draw();
			plot.setColor(Color.GRAY);
			for (int frame = 0; frame < totalFrame - 1; frame++) {
				if (TraceFitted[1][frame] > 0 || TraceFitted[1][frame + 1] > 0) {
					plot.drawLine(x[frame], TraceFitted[0][frame], x[frame + 1], TraceFitted[0][frame + 1]);
				}
			}
			plot.draw();
			SRutil.saveTiff(plot.getImagePlus(), param.fileDirName + File.separator + "Temp" + File.separator
					+ "Trace_Group-" + (group + 1) + ".tif", true);
		}

		return new int[] { numberOfGroups, totalIntensity };
	}

	public static void changeGroupFromTo(ResultsTableMt rt, int previousGroup, int newGroup) {
		double[] vectorGroup = rt.getColumnAsDoubles(ResultsTableMt.GROUP);
		if (vectorGroup == null)
			return;

		for (int fit = 0; fit < rt.getCounter(); fit++) {
			if (vectorGroup[fit] - 1 == previousGroup)
				rt.setValue(ResultsTableMt.GROUP, fit, newGroup + 1);
		}
	}

	public static double getPrecisionHistogram(ResultsTableMt rt, Params param, boolean show) {
		return getPrecisionHistogram(rt, param, false, show);
	}

	public static double getPrecisionHistogram(ResultsTableMt rt, Params param, boolean Zdirection,
			boolean show) {
		int SigmaRebuilt = ResultsTableMt.SIGMArebuilt;
		if (Zdirection)
			SigmaRebuilt = ResultsTableMt.SIGMArebuiltZ;

		if (rt == null || rt.getCounter() == 0)
			return 0;

		double[] MinMax = getMinMax(rt, SigmaRebuilt);
		double[] bin = new double[1001];
		double PixelSize = param.pixelSize;
		if (Zdirection)
			PixelSize = param.pixelSizeZ;
		for (int i = 0; i < bin.length; i++)
			bin[i] = (MinMax[0] + (MinMax[1] - MinMax[0]) * i / (bin.length - 1)) * PixelSize;

		double[] hist = new double[bin.length];
		for (int row = 0; row < rt.getCounter(); row++)
			hist[Math
					.min((int) ((rt.getValueAsDouble(SigmaRebuilt, row) - MinMax[0]) * (bin.length - 1) / (MinMax[1] - MinMax[0])),
							hist.length - 1)]++;
		int max = 0;
		for (int i = 1; i < bin.length; i++)
			if (hist[i] > hist[max])
				max = i;

		Plot histogram = new Plot("Precision histogram (max: " + (int) (bin[max]) + " nm)", "Bins (nm)",
				"Number of PSFs");
		histogram.setLimits(0, bin[bin.length - 1], 0, hist[max]);
		histogram.addPoints(bin, hist, Plot.LINE);
		histogram.draw();
		histogram.setColor(Color.GRAY);
		histogram.setLineWidth(2);
		double[] histRunningAverage = runningAverage(hist, 10);
		histogram.addPoints(bin, histRunningAverage, Plot.LINE);
		histogram.draw();
		if (show)
			histogram.show();
		String saveName;
		if (new File(param.fileDirName).isDirectory())
			saveName = param.fileDirName;
		else
			saveName = param.rootDirName;
		if (Zdirection)
			saveName += File.separator + "10_PrecisionZHistogram.tif";
		else
			saveName += File.separator + "10_PrecisionHistogram.tif";
		saveTiff(histogram.getImagePlus(), saveName, false);

		max = 0;
		for (int i = 1; i < bin.length; i++)
			if (histRunningAverage[i] > histRunningAverage[max])
				max = i;
		return bin[max];
	}

	public static void addRow(ResultsTableMt from, ResultsTableMt to, int row) {
		to.incrementCounter();
		for (int column = 0; column <= from.getLastColumn(); column++) {
			if(from.columnExists(column) && from.getColumnHeading(column) != null)
				to.addValue(from.getColumnHeading(column), from.getValueAsDouble(column, row));
		}
	}

	public static ResultsTableMt changeRtNamesSFLtoMP(ResultsTableMt steveRt, Params param) {
		ResultsTableMt mtRt = new ResultsTableMt();
		for (int row = 0; row < steveRt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, steveRt.getValue("Frame", row));
			mtRt.addValue(ResultsTableMt.X0, steveRt.getValue("X Center", row));
			mtRt.addValue(ResultsTableMt.Y0, steveRt.getValue("Y Center", row));
			mtRt.addValue(ResultsTableMt.SIGMAX, steveRt.getValue("Width", row));
			mtRt.addValue(ResultsTableMt.SIGMAY, steveRt.getValue("Width", row));
			mtRt.addValue(ResultsTableMt.I0, steveRt.getValue("Amplitude", row));
			mtRt.addValue(ResultsTableMt.NOISE, steveRt.getValue("Offset", row));
			mtRt.addValue(ResultsTableMt.IS_FITTED, steveRt.getValue("Good Fit?", row));
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, row));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
		}

		return mtRt;
	}

	public static MemoryPeakResults getPeakFitResultsFromRt(ResultsTableMt rt) {
		MemoryPeakResults results = new MemoryPeakResults(rt.getCounter());
		int[] rtIndexes = { ResultsTableMt.X0, ResultsTableMt.Y0, ResultsTableMt.SIGMAX,
				ResultsTableMt.SIGMAY, ResultsTableMt.NOISE, ResultsTableMt.FRAME, ResultsTableMt.I0,
				ResultsTableMt.IS_FITTED, ResultsTableMt.OFFSET };
		double[][] values = getColumnsFrom(rt, rtIndexes);

		for (int psf = 0; psf < rt.getCounter(); psf++) {
			float[] param = new float[7];
			param[Gaussian2DFunction.X_POSITION] = (float) values[0][psf];
			param[Gaussian2DFunction.Y_POSITION] = (float) values[1][psf];
			param[Gaussian2DFunction.X_WIDTH] = (float) values[2][psf];
			param[Gaussian2DFunction.Y_WIDTH] = (float) values[3][psf];
			param[Gaussian2DFunction.AMPLITUDE] = (float) getAmplitude(rt, psf);
			param[Gaussian2DFunction.BACKGROUND] = (float) values[4][psf];
			results.add((int) values[5][psf], (int) (values[0][psf] + 0.5), (int) (values[1][psf] + 0.5),
					(float) values[6][psf], (float) values[7][psf], (float) values[8][psf], param, null);
		}

		return results;
	}

	public static ResultsTableMt changeRtNamesPeakFitToMP(String file, Params param) {
		// Remove the 7 first lines from the PeakFit results
		File tempFile = null;
		try {
			File inputFile = new File(file);
			tempFile = new File(param.fileDirName + File.separator + "TableFit_PeakFit.Temp.txt");

			BufferedReader reader = new BufferedReader(new FileReader(inputFile));
			BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));

			String currentLine;
			int i = 1;

			while ((currentLine = reader.readLine()) != null) {
				if (i > 7)
					writer.write(currentLine + "\n");
				// if(i>7) IJ.log(""+currentLine);
				i++;
			}
			writer.close();
			reader.close();
		} catch (Exception e) {
		}

		// Open the new trimmed result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(param.fileDirName + File.separator
				+ "TableFit_PeakFit.Temp.txt");

		int[] mtRtIndexes = { ResultsTableMt.FRAME, ResultsTableMt.X0, ResultsTableMt.Y0,
				ResultsTableMt.SIGMAX, ResultsTableMt.SIGMAY, ResultsTableMt.I0, ResultsTableMt.NOISE,
				ResultsTableMt.IS_FITTED, ResultsTableMt.SIGMArebuilt, ResultsTableMt.AMPLITUDE };

		String[] rtHeadings = { "#Frame", "X", "Y", "X SD", "Y SD", "Signal", "Noise", "Error" };
		double[][] values = getColumnsFrom(rt, rtHeadings);

		for (int row = 0; row < rt.getCounter(); row++) {
			int column = 0;
			mtRt.incrementCounter();
			for (; column < rtHeadings.length; column++)
				mtRt.addValue(mtRtIndexes[column], values[column][row]);
			mtRt.addValue(mtRtIndexes[column], sigmaRebuilt(param, mtRt, mtRt.getCounter() - 1));
			mtRt.addValue(mtRtIndexes[column + 1], getAmplitude(mtRt, mtRt.getCounter() - 1));
		}

		tempFile.delete();
		return mtRt;
	}

	public static ResultsTableMt changeRtNamesQuickPALMToMP(String File, Params param) {

		// Open the new trimmed result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(File);

		for (int row = 0; row < rt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, rt.getValue("Frame Number", row));
			mtRt.addValue(ResultsTableMt.X0, rt.getValue("X (px)", row));
			mtRt.addValue(ResultsTableMt.Y0, rt.getValue("Y (px)", row));
			mtRt.addValue(ResultsTableMt.SIGMAX,
					0.5 * (rt.getValue("Left-Width(px)", row) + rt.getValue("Right-Width (px)", row)));
			mtRt.addValue(ResultsTableMt.SIGMAY,
					0.5 * (rt.getValue("Up-Height (px)", row) + rt.getValue("Down-Height (px)", row)));
			mtRt.addValue(ResultsTableMt.I0, rt.getValue("Intensity", row));
			mtRt.addValue(ResultsTableMt.NOISE, 0);
			mtRt.addValue(ResultsTableMt.IS_FITTED, 1);
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, row));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
		}

		return mtRt;
	}

	public static ResultsTableMt changeRtNamesrapidSTORMToMP(String file, Params param) {
		// Remove the first line from the rapidSTORM results
		File tempFile = null;
		try {
			File inputFile = new File(file);
			tempFile = new File(param.fileDirName + File.separator + "TableFit_rapidSTORM.Temp.txt");

			BufferedReader reader = new BufferedReader(new FileReader(inputFile));
			BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));

			String currentLine;
			int i = 1;

			while ((currentLine = reader.readLine()) != null) {
				if (i == 1)
					writer.write("X_nm" + "\t" + "Y_nm" + "\t" + "Frame" + "\t" + "Sigma_X_nm" + "\t"
							+ "Sigma_Y_nm" + "\t" + "Intensity" + "\t" + "Chi_Squared" + "\t" + "Noise"
							+ "\n");
				if (i > 1)
					writer.write(currentLine.replace(" ", "\t") + "\n");
				// if(i>7) IJ.log(""+currentLine);
				i++;
			}
			writer.close();
			reader.close();
		} catch (Exception e) {
		}

		// Open the new trimmed result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(param.fileDirName + File.separator
				+ "TableFit_rapidSTORM.Temp.txt");

		for (int row = 0; row < rt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, rt.getValue("Frame", row));
			mtRt.addValue(ResultsTableMt.X0, rt.getValue("X_nm", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.Y0, rt.getValue("Y_nm", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.SIGMAX,
					fwhmToStdev(rt.getValue("Sigma_Y_nm", row) / param.pixelSize));
			mtRt.addValue(ResultsTableMt.SIGMAY,
					fwhmToStdev(rt.getValue("Sigma_Y_nm", row) / param.pixelSize));
			mtRt.addValue(ResultsTableMt.I0, rt.getValue("Intensity", row));
			mtRt.addValue(ResultsTableMt.NOISE, rt.getValue("Chi_Squared", row));
			mtRt.addValue(ResultsTableMt.IS_FITTED, rt.getValue("Noise", row));
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, mtRt.getCounter() - 1));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
			// IJ.log(""+row);
		}

		tempFile.delete();

		return mtRt;
	}

	public static ResultsTableMt changeRtNamesThunderSTORMToMP(String file, Params param) {
		// Open the result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(file);

		for (int row = 0; row < rt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, rt.getValue("\"frame\"", row));
			mtRt.addValue(ResultsTableMt.X0, rt.getValue("\"x [nm]\"", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.Y0, rt.getValue("\"y [nm]\"", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.SIGMAX, rt.getValue("\"sigma [nm]\"", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.SIGMAY, rt.getValue("\"sigma [nm]\"", row) / param.pixelSize);
			mtRt.addValue(ResultsTableMt.I0,
					translatePhotonsToCounts(rt.getValue("\"intensity [photon]\"", row), param, false));
			mtRt.addValue(ResultsTableMt.NOISE,
					translatePhotonsToCounts(rt.getValue("\"bkgstd [photon]\"", row), param, false));
			mtRt.addValue(ResultsTableMt.IS_FITTED, rt.getValue("\"uncertainty [nm]\"", row));
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, mtRt.getCounter() - 1));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
			// IJ.log(""+row);
		}
		return mtRt;
	}

	public static ResultsTableMt changeRtNamesM2LEToMP(String file, Params param) {
		// Open the result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(file);

		for (int row = 0; row < rt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, rt.getValue("Frame", row));
			mtRt.addValue(ResultsTableMt.X0, rt.getValue("x (px)", row));
			mtRt.addValue(ResultsTableMt.Y0, rt.getValue("y (px)", row));
			mtRt.addValue(ResultsTableMt.SIGMAX, fwhmToStdev(rt.getValue("Width x", row)));
			mtRt.addValue(ResultsTableMt.SIGMAY, fwhmToStdev(rt.getValue("Width y", row)));
			mtRt.addValue(ResultsTableMt.I0,
					(rt.getValue("Intensity x", row) + rt.getValue("Intensity y", row)) * 5e5D);
			mtRt.addValue(ResultsTableMt.NOISE,
					rt.getValue("Background x", row) + rt.getValue("Background y", row));
			mtRt.addValue(ResultsTableMt.IS_FITTED, 1);
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, mtRt.getCounter() - 1));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
			// IJ.log(""+row);
		}
		return mtRt;
	}

	public static ResultsTableMt changeRtNamesGraspJToMP(String file, Params param) {
		// Open the result file as a ResultsTable
		ResultsTableMt mtRt = new ResultsTableMt();
		ResultsTableMt rt = ResultsTableMt.open2(file);

		for (int row = 0; row < rt.getCounter(); row++) {
			mtRt.incrementCounter();
			mtRt.addValue(ResultsTableMt.FRAME, rt.getValue("Frame", row));
			mtRt.addValue(ResultsTableMt.X0, rt.getValue("x (px)", row));
			mtRt.addValue(ResultsTableMt.Y0, rt.getValue("y (px)", row));
			mtRt.addValue(ResultsTableMt.SIGMAX, rt.getValue("Width x", row));
			mtRt.addValue(ResultsTableMt.SIGMAY, rt.getValue("Width y", row));
			mtRt.addValue(
					ResultsTableMt.I0,
					translatePhotonsToCounts(
							(rt.getValue("Intensity x", row) + rt.getValue("Intensity y", row)), param, false));
			mtRt.addValue(
					ResultsTableMt.NOISE,
					translatePhotonsToCounts(
							(rt.getValue("Background x", row) + rt.getValue("Background y", row)), param,
							false));
			mtRt.addValue(ResultsTableMt.IS_FITTED, 1);
			mtRt.addValue(ResultsTableMt.SIGMArebuilt, sigmaRebuilt(param, mtRt, mtRt.getCounter() - 1));
			mtRt.addValue(ResultsTableMt.AMPLITUDE, getAmplitude(mtRt, mtRt.getCounter() - 1));
			// IJ.log(""+row);
		}
		return mtRt;
	}

	public static double fwhmToStdev(double fwhm) {
		final double factor = Math.sqrt(8.0 * Math.log(2));
		return fwhm / factor;
	}

	public static double stdevtoFWHM(double sigma) {
		final double factor = Math.sqrt(8.0 * Math.log(2));
		return factor * sigma;
	}

	public static double[] getColumnFromRt(ResultsTableMt rt, String x0) {
		if (rt == null)
			return null;

		int index = rt.getColumnIndex(x0);
		if (index == ResultsTable.COLUMN_NOT_FOUND)
			return null;

		return rt.getColumnAsDoubles(index);
	}

	public static int[] getIndexOf(ResultsTableMt rt, String[] headings) {
		if (rt == null)
			return null;

		int[] indexes = new int[headings.length];
		for (int column = 0; column < headings.length; column++)
			indexes[column] = rt.getColumnIndex(headings[column]);

		return indexes;
	}

	public static int[] addHeadingsTo(ResultsTableMt rt, String[] headings) {
		if (rt == null)
			return null;

		int[] indexes = getIndexOf(rt, headings);
		boolean deleteRow0 = false;
		if (rt.getCounter() == 0) {
			rt.incrementCounter();
			deleteRow0 = true;
		}

		for (int column = 0; column < headings.length; column++) {
			if (indexes[column] == ResultsTable.COLUMN_NOT_FOUND) {
				rt.addValue(headings[column], 0);
				indexes[column] = rt.getColumnIndex(headings[column]);
			}
		}

		if (deleteRow0)
			rt.deleteRow(0);

		return indexes;
	}

	public static double[][] getColumnsFrom(ResultsTableMt rt, String[] headings) {
		if (rt == null)
			return null;

		double[][] values = new double[headings.length][];
		int[] indexes = getIndexOf(rt, headings);
		for (int column = 0; column < headings.length; column++) {
			if (indexes[column] != ResultsTable.COLUMN_NOT_FOUND)
				values[column] = rt.getColumnAsDoubles(indexes[column]);
		}

		return values;
	}

	public static double[][] getColumnsFrom(ResultsTableMt rt, int[] indexes) {
		if (rt == null)
			return null;

		double[][] values = new double[indexes.length][];
		for (int column = 0; column < indexes.length; column++)
			values[column] = rt.getColumnAsDoubles(indexes[column]);

		return values;
	}

	public static ResultsTableMt getTableMaxFromTableFit(ResultsTableMt tableFit) {
		ResultsTableMt tableMax = new ResultsTableMt();
		int[] tableMaxIndexes = { ResultsTableMt.FRAME, ResultsTableMt.Xmax, ResultsTableMt.Ymax,
				ResultsTableMt.Imax };// addHeadingsTo(tableMax,
		// tableMaxHeadings);

		int[] tableFitIndexes = { ResultsTableMt.FRAME, ResultsTableMt.X0, ResultsTableMt.Y0,
				ResultsTableMt.I0 };
		double[][] values = getColumnsFrom(tableFit, tableFitIndexes);

		for (int row = 0; row < tableFit.getCounter(); row++) {
			tableMax.incrementCounter();
			for (int column = 0; column < tableMaxIndexes.length; column++)
				tableMax.addValue(tableMaxIndexes[column], (int) values[column][row]);
		}

		return tableMax;
	}

	public static ResultsTableMt concatenate(ResultsTableMt rt1, ResultsTableMt rt2) {
		if (rt2 == null)
			return rt1;

		IJ.showStatus("Combining results...");
		String[] rt1Headings = rt1.getHeadings();
		double[][] rt2Values = getColumnsFrom(rt2, rt1Headings);

		for (int i = 0; i < rt2.getCounter(); i++) {
			if (i % 100 == 0)
				IJ.showProgress(i, rt2.getCounter());

			rt1.incrementCounter();
			for (int j = 0; j <= rt1.getLastColumn(); j++)
				if (rt2Values[j] != null)
					rt1.addValue(j, rt2Values[j][i]);
		}

		IJ.showProgress(1);

		return rt1;
	}

	public static double[] concatenate(double[] vector1, double[] vector2) {
		double[] result = new double[vector1.length + vector2.length];
		for (int i = 0; i < vector1.length; i++)
			result[i] = vector1[i];
		for (int i = 0; i < vector2.length; i++)
			result[vector1.length + i] = vector2[i];
		return result;
	}

	public static int getRealFrameNumber(ResultsTableMt rt, int rowInRt, int frameNumberPerCycle) {
		if (rt.columnExists(ResultsTableMt.CYCLE))
			return (int) (rt.getValueAsDouble(ResultsTableMt.FRAME, rowInRt) + (rt.getValueAsDouble(
					ResultsTableMt.CYCLE, rowInRt) - 2) * frameNumberPerCycle);
		else
			return (int) (rt.getValueAsDouble(ResultsTableMt.FRAME, rowInRt));
	}

	public static double getMax(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getMinMax(rt.getColumnAsDoubles(rt.getColumnIndex(x0)))[1];
		else
			return -1;
	}

	public static double getMax(ResultsTableMt rt, int x0) {
		if (x0 >= 0)
			return getMinMax(rt.getColumnAsDoubles(x0))[1];
		else
			return -1;
	}

	public static double getMin(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getMinMax(rt.getColumnAsDoubles(rt.getColumnIndex(x0)))[0];
		else
			return -1;
	}

	public static double getMin(ResultsTableMt rt, int x0) {
		if (x0 >= 0)
			return getMinMax(rt.getColumnAsDoubles(x0))[0];
		else
			return -1;
	}

	public static double[] getMinMax(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getMinMax(rt.getColumnAsDoubles(rt.getColumnIndex(x0)));
		else
			return null;
	}

	public static double[] getMinMax(ResultsTableMt rt, int x0) {
		return getMinMax(rt.getColumnAsDoubles(x0));
	}

	public static double getIQR(double[] vector, boolean discardZeros) {
		double[] ordered = vector.clone();
		Arrays.sort(ordered);

		if (discardZeros) {
			int row = 0;
			while (row < ordered.length && ordered[row] == 0)
				row++;
			if (row == ordered.length)
				return vector.length;
			else
				return ordered[(int) (0.75 * (double) (ordered.length - row)) + row + 1]
						- ordered[(int) (0.25 * (double) (ordered.length - row)) + row];
		} else
			return ordered[Math.min((int) (0.75 * (double) ordered.length) + 1, ordered.length - 1)]
					- ordered[(int) (0.25 * (double) ordered.length)];
	}

	public static double getIQR(double[] vector) {
		return getIQR(vector, false);
	}

	public static double getIQR(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getIQR(rt.getColumnAsDoubles(rt.getColumnIndex(x0)));
		else
			return -1;
	}

	public static double getIQR(ResultsTableMt rt, int x0) {
		return getIQR(rt.getColumnAsDoubles(x0));
	}

	public static double getMean(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getMean(rt.getColumnAsDoubles(rt.getColumnIndex(x0)));
		else
			return -1;
	}

	public static double getMean(ResultsTableMt rt, int x0) {
		return getMean(rt.getColumnAsDoubles(x0));
	}

	public static double[] getMeanAndStdev(ResultsTableMt rt, String x0) {
		if (rt.getColumnIndex(x0) != -1)
			return getMeanAndStdev(rt.getColumnAsDoubles(rt.getColumnIndex(x0)));
		else
			return null;
	}

	public static double getMean(double[] vector, boolean discardZeros) {
		double mean = 0;
		int tot = 0;
		for (int row = 0; row < vector.length; row++) {
			if (!discardZeros || vector[row] != 0) {
				mean += vector[row];
				tot++;
			}
		}
		return mean / (double) tot;
	}

	public static double getMean(double[] vector) {
		return getMean(vector, false);
	}

	public static double[] getMeanAndStdev(double[] vector) {
		double[] meanAndStdev = new double[2];
		meanAndStdev[0] = getMean(vector);
		for (int row = 0; row < vector.length; row++) {
			meanAndStdev[1] += Math.pow(vector[row] - meanAndStdev[0], 2);
		}
		meanAndStdev[1] = Math.sqrt(meanAndStdev[1] / (double) vector.length);
		return meanAndStdev;
	}

	public static double getMax(double[] vector) {
		return getMinMax(vector)[1];
	}

	public static double getMin(double[] vector) {
		return getMinMax(vector)[0];
	}

	public static ResultsTableMt getMinMaxS(double[] vector) {
		ResultsTableMt minMaxS = new ResultsTableMt();
		int minNumber = 0;
		int maxNumber = 0;

		double[] derivative = derivative(vector);
		minMaxS.incrementCounter();
		if (derivative[0] > 0) {
			minMaxS.addValue("xmin", 0);
			minMaxS.addValue("Min", vector[0]);
			minNumber++;
		} else {
			minMaxS.addValue("xmax", 0);
			minMaxS.addValue("Max", vector[0]);
			maxNumber++;
		}
		for (int i = 1; i < vector.length; i++) {
			if (derivative[i] * derivative[i - 1] < 0) {
				if (derivative[i - 1] > 0) {
					while (minMaxS.getCounter() <= maxNumber)
						minMaxS.incrementCounter();
					minMaxS.setValue("xmax", maxNumber, i);
					minMaxS.setValue("Max", maxNumber, vector[i]);
					maxNumber++;
				} else {
					while (minMaxS.getCounter() <= minNumber)
						minMaxS.incrementCounter();
					minMaxS.setValue("xmin", minNumber, i);
					minMaxS.setValue("Min", minNumber, vector[i]);
					minNumber++;
				}
			}
		}

		return minMaxS;
	}

	public static double[] getMinMax(double[] vector) {
		double[] minMax = { Math.pow(2, 128), 0 };
		for (int row = 0; row < vector.length; row++) {
			if (vector[row] > minMax[1])
				minMax[1] = vector[row];
			if (vector[row] < minMax[0])
				minMax[0] = vector[row];
		}
		return minMax;
	}

	public static int getMax(int[] vector) {
		int max = 0;
		for (int row = 0; row < vector.length; row++) {
			if (vector[row] > max)
				max = vector[row];
		}
		return max;
	}

	public static int getMax(int[][] matrix) {
		int max = 0;
		for (int column = 0; column < matrix.length; column++) {
			int temp = getMax(matrix[column]);
			if (temp > max)
				max = temp;
		}
		return max;
	}

	public static int getMin(int[] vector) {
		int min = (int) Math.pow(2, 128);
		for (int row = 0; row < vector.length; row++) {
			if (vector[row] < min)
				min = vector[row];
		}
		return min;
	}

	public static int getMin(int[][] matrix) {
		int min = (int) Math.pow(2, 128);
		for (int column = 0; column < matrix.length; column++) {
			int temp = getMin(matrix[column]);
			if (temp < min)
				min = temp;
		}
		return min;
	}

	public static int[] doubleToInt(double[] vector) {
		int[] retour = new int[vector.length];
		for (int i = 0; i < vector.length; i++) {
			retour[i] = (int) (vector[i] + 0.5D);
		}
		return retour;
	}

	public static double[] intToDouble(int[] vector) {
		double[] retour = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			retour[i] = vector[i];
		}
		return retour;
	}

	public static final int VoronoiDensity = 10000;
	public static final int CircleDensity = 10001;
	public static final int SquareDensity = 10002;

	public static double linearInterpolation(double xi, double yi, ResultsTableMt closeBy) {
		// http://en.wikipedia.org/wiki/Barycentric_coordinate_system

		double det = (closeBy.getValueAsDouble(ResultsTableMt.Y0, 1) - closeBy.getValueAsDouble(
				ResultsTableMt.Y0, 2))
				* (closeBy.getValueAsDouble(ResultsTableMt.X0, 0) - closeBy.getValueAsDouble(
						ResultsTableMt.X0, 2))
				+ (closeBy.getValueAsDouble(ResultsTableMt.X0, 2) - closeBy.getValueAsDouble(
						ResultsTableMt.X0, 1))
				* (closeBy.getValueAsDouble(ResultsTableMt.Y0, 0) - closeBy.getValueAsDouble(
						ResultsTableMt.Y0, 2));
		if (det == 0) {
			closeBy.incrementCounter();
			closeBy.addValue(ResultsTableMt.X0, xi);
			closeBy.addValue(ResultsTableMt.Y0, yi);
			double d14 = getDistanceBetweenPix(closeBy, closeBy.getCounter() - 1, 0);
			double d24 = getDistanceBetweenPix(closeBy, closeBy.getCounter() - 1, 1);
			double d34 = getDistanceBetweenPix(closeBy, closeBy.getCounter() - 1, 2);
			double max = Math.max(d14, Math.max(d24, d34));
			if (max == d14)
				return (d24 * closeBy.getValue("Density", 1) + d34 * closeBy.getValue("Density", 2))
						/ (d24 + d34);
			else if (max == d24)
				return (d14 * closeBy.getValue("Density", 0) + d34 * closeBy.getValue("Density", 2))
						/ (d14 + d34);
			else if (max == d34)
				return (d14 * closeBy.getValue("Density", 0) + d24 * closeBy.getValue("Density", 1))
						/ (d14 + d24);
			else
				return closeBy.getValue("Density", 0);
		} else {
			double lambda1 = ((closeBy.getValueAsDouble(ResultsTableMt.Y0, 1) - closeBy.getValueAsDouble(
					ResultsTableMt.Y0, 2)) * (xi - closeBy.getValueAsDouble(ResultsTableMt.X0, 2)) + (closeBy
					.getValueAsDouble(ResultsTableMt.X0, 2) - closeBy.getValueAsDouble(ResultsTableMt.X0, 1))
					* (yi - closeBy.getValueAsDouble(ResultsTableMt.Y0, 2))) / det;
			double lambda2 = ((closeBy.getValueAsDouble(ResultsTableMt.Y0, 2) - closeBy.getValueAsDouble(
					ResultsTableMt.Y0, 0)) * (xi - closeBy.getValueAsDouble(ResultsTableMt.X0, 2)) + (closeBy
					.getValueAsDouble(ResultsTableMt.X0, 0) - closeBy.getValueAsDouble(ResultsTableMt.X0, 2))
					* (yi - closeBy.getValueAsDouble(ResultsTableMt.Y0, 2))) / det;
			double lambda3 = 1.0D - lambda1 - lambda2;
			// IJ.log("lamda1: "+lambda1+" / lambda2: "+lambda2);
			return (lambda1 * closeBy.getValue("Density", 0) + lambda2 * closeBy.getValue("Density", 1) + lambda3
					* closeBy.getValue("Density", 2));
		}
	}

	public static double[][] closestNeighbours(double[] x, double[] y, double xi, double yi, int neighbours) {
		ResultsTableMt rt = new ResultsTableMt();

		for (int i = 0; i < Math.min(x.length, y.length); i++) {
			rt.incrementCounter();
			rt.addValue(ResultsTableMt.X0, x[i]);
			rt.addValue(ResultsTableMt.Y0, y[i]);
		}

		rt.incrementCounter();
		rt.addValue(ResultsTableMt.X0, xi);
		rt.addValue(ResultsTableMt.Y0, yi);

		return closestNeighbours(rt, rt.getCounter() - 1, neighbours);
	}

	public static double[] closestNeighbour(ResultsTableMt rt, int psf) {
		return closestNeighbours(rt, psf, 1)[0];
	}

	public static double[][] closestNeighbours(ResultsTableMt rt, int psf, int neighbours) {
		double[][] retour = new double[neighbours][2];
		for (int n = 0; n < neighbours; n++)
			retour[n][0] = 10000;
		double temp;

		for (int i = 0; i < rt.getCounter(); i++) {
			if (i != psf) {
				temp = getDistanceBetweenPix(rt, i, psf);

				for (int n = neighbours - 1; n >= 0; n--) {
					if (n == neighbours - 1 && temp < retour[n][0]) {
						retour[n][0] = temp;
						retour[n][1] = i;
					} else if (n < neighbours - 1 && retour[n + 1][0] < retour[n][0]) {
						temp = retour[n][0];
						retour[n][0] = retour[n + 1][0];
						retour[n + 1][0] = temp;

						temp = retour[n][1];
						retour[n][1] = retour[n + 1][1];
						retour[n + 1][1] = temp;
					} else
						break;
				}
			}
		}

		return retour;
	}

	public static double getDistanceBetweenPix(ResultsTableMt rt, int rowInRt1, int rowInRt2) {
		int x0 = ResultsTableMt.X0, y0 = ResultsTableMt.Y0;
		if(!rt.columnExists(x0)) {
			if(rt.columnExists(ResultsTableMt.Xmax)) {
				x0 = ResultsTableMt.Xmax;
				y0 = ResultsTableMt.Ymax;
			} else if(rt.columnExists(ResultsTableMt.X_CENTROID)) {
				x0 = ResultsTableMt.X_CENTROID;
				y0 = ResultsTableMt.Y_CENTROID;
			} else if(rt.columnExists(ResultsTableMt.X_CENTER_OF_MASS)) {
				x0 = ResultsTableMt.X_CENTER_OF_MASS;
				y0 = ResultsTableMt.Y_CENTER_OF_MASS;
			}
		}
		
		return Math.sqrt(Math.pow(
				rt.getValueAsDouble(x0, rowInRt1)
						- rt.getValueAsDouble(x0, rowInRt2), 2.0)
				+ Math.pow(
						rt.getValueAsDouble(y0, rowInRt1)
								- rt.getValueAsDouble(y0, rowInRt2), 2.0));
	}

	public static ImagePlus rebuildingPeakFit(ResultsTableMt rt, String nameOfTheFile, Params param,
			ResultsSettings resultsSettings) {
		if (resultsSettings.getResultsImage() == null)
			resultsSettings.setResultsImage(ResultsImage.LOCALISATIONS_PRECISION);
		resultsSettings.showDeviations = false;
		Rectangle rec = param.roi;

		// Display the configured output
		IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
				resultsSettings.getResultsImage(), resultsSettings.weightedImage,
				resultsSettings.equalisedImage, nameOfTheFile, rec, param.pixelSize, 1f, resultsSettings.imageScale,
				resultsSettings.precision, ResultsMode.ADD);
		image.setRollingWindowSize(resultsSettings.imageRollingWindow);

		MemoryPeakResults results = getPeakFitResultsFromRt(rt);

		image.begin();
		final int totalSteps = results.size();
		final int step = (totalSteps > 400) ? totalSteps / 200 : 2;
		int slice = 0;
		for (PeakResult r : results.getResults()) {
		    if (++slice % step == 0)
		        IJ.showProgress(slice, totalSteps);
		    image.add(r.peak, r.origX, r.origY, r.origValue, r.error, r.noise, r.params, r.paramsStdDev);
		}
		IJ.showProgress(1);
		image.end();

		if (param.scaleBar) {
			int ScaleBarSize = (int) Math.round(1000.0 / param.pixelSize * param.rebuiltScale);
			if (ScaleBarSize + 4 * param.rebuiltScale < image.getImagePlus().getWidth()) {
				image.getImagePlus().setColor(Color.WHITE);
				image.getImagePlus()
						.getProcessor()
						.fill(new Roi(new Rectangle(2 * param.rebuiltScale, (param.roi.height - 2)
								* param.rebuiltScale, ScaleBarSize, param.rebuiltScale)));
				image.getImagePlus().updateImage();
			}
		}

		String FileName = param.fileDirName + File.separator + "RebuiltSRPicture_PeakFit" + nameOfTheFile
				+ ".tif";
		saveTiff(image.getImagePlus(), FileName, false);

		return image.getImagePlus();

	}

	public static ImagePlus rebuilding(ResultsTableMt rtGrouped, String nameOfTheFile, Params param) {

		if (!param.hasBeenInitialised)
			param.getNewParameters();

		// Preparing the Rebuilt image and window
		int rebuiltScale = param.rebuiltScale;
		if (param.rebuilding == Params.Rebuilding.DiffractionLimited)
			rebuiltScale = 1;
		double pixelSize = param.pixelSize;
		boolean colouring = param.rebuilding == Params.Rebuilding.ColouringCirclePSF;
		boolean boolean_ = (param.rebuilding == Params.Rebuilding.BooleanAddPSF || param.rebuilding == Params.Rebuilding.BooleanReplacePSF);
		ImageProcessor ipRebuilt;

		if (!nameOfTheFile.equalsIgnoreCase("None")) {
			if (colouring)
				nameOfTheFile += "-TimeStamp";
			else if (boolean_)
				nameOfTheFile += "-Boolean";
			else if (param.rebuilding == Params.Rebuilding.ColorCodedPerGroup)
				nameOfTheFile += "-ShowGroups";
			else if (param.rebuilding == Params.Rebuilding.DiffractionLimited)
				nameOfTheFile += "-DiffractionLimited";
			else
				nameOfTheFile += "-Proba";
		}
		String FileName = param.fileDirName + File.separator + "RebuiltSRPicture" + nameOfTheFile;

		if (param.roi == null) {
			if (new File(param.rootDirName + File.separator + param.fileName).exists()) {
				ImagePlus Temp = new Opener().openImage(param.rootDirName + File.separator + param.fileName,
						1);
				param.roi = new Rectangle(0, 0, Temp.getWidth(), Temp.getHeight());
			} else {
				param.roi = new Rectangle(0, 0, (int) (getMax(rtGrouped, ResultsTableMt.X0) + 10.5),
						(int) (getMax(rtGrouped, ResultsTableMt.Y0) + 10.5));
			}
		}

		if (colouring || param.rebuilding == Params.Rebuilding.ColorCodedPerGroup) { // Showing
			// the legend
			ipRebuilt = new ShortProcessor((rebuiltScale) * (param.roi.width + 2), rebuiltScale
					* param.roi.height);

			for (int k = 1 * rebuiltScale; k < (param.roi.height - 2) * rebuiltScale; k++) {
				for (int l = 0; l < rebuiltScale; l++) {
					ipRebuilt.set(l + param.roi.width * rebuiltScale, k, (int) ((k - 1.0 * rebuiltScale)
							/ ((param.roi.height - 3.0) * rebuiltScale) * 65535.0));
				}
			}

			if (param.rebuilding == Params.Rebuilding.ColorCodedPerGroup)
				maxGroup = getMax(rtGrouped, ResultsTableMt.GROUP);
		} else {
			if (param.rebuilding == Params.Rebuilding.DiffractionLimited)
				ipRebuilt = /* Float */new ShortProcessor((rebuiltScale) * param.roi.width, rebuiltScale
						* param.roi.height);
			else
				ipRebuilt = new ShortProcessor((rebuiltScale) * param.roi.width, rebuiltScale
						* param.roi.height);
		}

		ImagePlus impRebuilt = new ImagePlus("Rebuilt super-resolved picture", ipRebuilt);
		WindowManager.setTempCurrentImage(impRebuilt);
		LutLoader lut = new LutLoader();
		if (colouring || param.rebuilding == Params.Rebuilding.ColorCodedPerGroup)
			lut.run("Red Hot"); // Or any other inbuilt LUT
		else
			lut.run("Blue Green Red");
		WindowManager.setTempCurrentImage(null);

		int frameLength = (int) getMax(rtGrouped.getColumnAsDoubles(ResultsTableMt.FRAME));

		scaleIntensity = 1.0D;
		llhFunct = new FittingClassical(param, ipRebuilt, 10, 10, 0, 0, 0, 0, 0, 0);
		llhFunct.setFast(true); // No need to have PixelPrecision>1
		for (int i = 0; i < rtGrouped.getCounter(); i++) {
			// Update the rebuilt picture
			if (i % 1000 == 0) {
				impRebuilt = addToRebuiltImp(rtGrouped, i, param, frameLength, impRebuilt, true);
			} else
				impRebuilt = addToRebuiltImp(rtGrouped, i, param, frameLength, impRebuilt, false);
		}

		if (param.scaleBar && !nameOfTheFile.equalsIgnoreCase("None")) {
			int ScaleBarSize = (int) Math.round(1000.0 / pixelSize * rebuiltScale);
			if (ScaleBarSize + 4 * param.rebuiltScale < ipRebuilt.getWidth()) {
				ipRebuilt.fill(new Roi(new Rectangle(2 * param.rebuiltScale, (param.roi.height - 2)
						* param.rebuiltScale, ScaleBarSize, param.rebuiltScale)));
			}
		}
		impRebuilt.setProcessor(ipRebuilt);
		impRebuilt.updateAndDraw();

		if (!nameOfTheFile.equalsIgnoreCase("None")) {
			saveTiff(impRebuilt, FileName + ".tif", false);

			if (boolean_) {
				double precisionForBlurring = 20.0; // in nm
				new GaussianBlur().blurGaussian(ipRebuilt, precisionForBlurring / pixelSize * rebuiltScale,
						precisionForBlurring / pixelSize * rebuiltScale, 0.01);
				ipRebuilt.setAutoThreshold(ImageProcessor.ISODATA2, ImageProcessor.NO_LUT_UPDATE);
				impRebuilt.updateAndDraw();
				saveTiff(impRebuilt, FileName + "_GaussianBlurred.tif", false);
			}

		}

		return impRebuilt;
	}

	private static double scaleIntensity = 1.0;
	private static double maxGroup = 0;
	private static FittingClassical llhFunct;
	static final double MAX = Math.pow(2, 15);

	public static ImagePlus addToRebuiltImp(ResultsTableMt rtGrouped, int i, Params param, int frameLength,
			ImagePlus impRebuilt, boolean refresh) {

		// if(RtGrouped.getValue(ResultsTableMt.IS_FITTED, i)>=1) {
		if (rtGrouped.getValueAsDouble(ResultsTableMt.IS_FITTED, i) > 0
				&& (param.sigmaThreshold < 0 || (rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAX, i) < param.sigmaThreshold && rtGrouped
						.getValueAsDouble(ResultsTableMt.SIGMAY, i) < param.sigmaThreshold))
				&& (param.i0Threshold < 0 || rtGrouped.getValueAsDouble(ResultsTableMt.I0, i) > param.i0Threshold)
				&& (param.ellipticityThreshold < 0 || Math.max(
						rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAX, i)
								/ rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAY, i),
						rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAY, i)
								/ rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAX, i)) < param.ellipticityThreshold)) {

			IJ.showProgress(i, rtGrouped.getCounter() - 1);
			IJ.showStatus("Rebuilding the final SR picture...");

			int rebuiltScale = param.rebuiltScale;
			if (param.rebuilding == Params.Rebuilding.DiffractionLimited)
				rebuiltScale = 1;

			ImageProcessor ipRebuilt = impRebuilt.getProcessor();
			short[] pixels = (short[]) ipRebuilt.getPixels();
			int index = ((int) (rtGrouped.getValueAsDouble(ResultsTableMt.Y0, i) * rebuiltScale) - param.roi.y
					* rebuiltScale)
					* ipRebuilt.getWidth()
					+ ((int) (rtGrouped.getValueAsDouble(ResultsTableMt.X0, i) * rebuiltScale) - param.roi.x
							* rebuiltScale);

			if (param.rebuilding == Params.Rebuilding.BooleanReplacePSF)
				pixels[index] = (short) 65000;
			else if (param.rebuilding == Params.Rebuilding.BooleanAddPSF)
				pixels[index] += 1000;
			else {
				double psfSigma = 1.323 / (param.pixelSize * 2.0 * Math.PI * param.na / param.lambda) * 2.0;
				int psfSigmaInt = (int) Math.max(Math.ceil(psfSigma), 2);

				double psfSigmaRebuilt = sigmaRebuilt(param, rtGrouped, i);
				int widthRebuiltBox = 3;
				if (param.rebuilding == Params.Rebuilding.DiffractionLimited)
					widthRebuiltBox *= 20;
				int xmin0 = (int) (rtGrouped.getValueAsDouble(ResultsTableMt.X0, i) - widthRebuiltBox
						* psfSigmaInt)
						* rebuiltScale;
				int ymin0 = (int) (rtGrouped.getValueAsDouble(ResultsTableMt.Y0, i) - widthRebuiltBox
						* psfSigmaInt)
						* rebuiltScale;
				int length0 = (2 * widthRebuiltBox * psfSigmaInt + 1) * rebuiltScale;
				double[] paramTemp = { rtGrouped.getValueAsDouble(ResultsTableMt.X0, i) * rebuiltScale,
						rtGrouped.getValueAsDouble(ResultsTableMt.Y0, i) * rebuiltScale,
						psfSigmaRebuilt * rebuiltScale, psfSigmaRebuilt * rebuiltScale,
						rtGrouped.getValueAsDouble(ResultsTableMt.I0, i) * Math.pow(rebuiltScale, 2), 0 };
				if (param.rebuilding == Params.Rebuilding.IdenticalGaussianPSF)
					paramTemp[4] = 5000;
				if (param.rebuilding == Params.Rebuilding.DiffractionLimited) {
					paramTemp[2] = rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAX, i) * rebuiltScale;
					paramTemp[3] = rtGrouped.getValueAsDouble(ResultsTableMt.SIGMAY, i) * rebuiltScale;
				}

				if (param.rebuilding == Params.Rebuilding.ColouringCirclePSF
						|| param.rebuilding == Params.Rebuilding.IdenticalCirclePSF
						|| param.rebuilding == Params.Rebuilding.ColorCodedPerGroup) {
					paramTemp[2] = param.fixedPSFSizeRebuilt;
					paramTemp[3] = param.fixedPSFSizeRebuilt;
					double add = 0;
					if (param.rebuilding == Params.Rebuilding.ColouringCirclePSF)
						add = rtGrouped.getValueAsDouble(ResultsTableMt.FRAME, i) / frameLength
								* (65535.0 - 12000.0) + 12000.0;
					else if (param.rebuilding == Params.Rebuilding.ColorCodedPerGroup)
						add = rtGrouped.getValueAsDouble(ResultsTableMt.GROUP, i) / maxGroup
								* (65535.0 - 12000.0) + 12000.0;
					else
						add = rtGrouped.getValueAsDouble(ResultsTableMt.I0, i);
					ipRebuilt.setColor((int) (add + 0.5D));
					OvalRoi roi = new OvalRoi((int) (paramTemp[0] + 0.5), (int) (paramTemp[1] + 0.5),
							(int) (2.0 * (paramTemp[2] + 0.5)), (int) (2.0 * (paramTemp[3] + 0.5)));
					ipRebuilt.fill(roi);

				} else {

					for (int ii = Math.max(xmin0, param.roi.x * rebuiltScale); ii < Math.min(xmin0 + length0,
							(param.roi.x + param.roi.width) * rebuiltScale); ii++) {
						for (int jj = Math.max(ymin0, param.roi.y * rebuiltScale); jj < Math.min(ymin0
								+ length0, (param.roi.y + param.roi.height) * rebuiltScale); jj++) {

							double add = llhFunct.PSFintegration(ii, jj, paramTemp) * scaleIntensity;
							index = (jj - param.roi.y * rebuiltScale) * ipRebuilt.getWidth() + ii
									- param.roi.x * rebuiltScale;

							if (ipRebuilt.getBitDepth() < 32) {
								while (pixels[index] + (int) add > MAX) {
									ipRebuilt.setPixels(pixels);
									ipRebuilt.multiply(0.5);
									pixels = (short[]) ipRebuilt.getPixels();

									scaleIntensity *= 0.5;
									add *= 0.5;
								}
							}

							pixels[index] += (int) add;
						}
					}
				}
			}
			ipRebuilt.setPixels(pixels);
			if (refresh)
				impRebuilt.updateAndDraw();
		}

		return impRebuilt;

	}

	// FILTER ON IMAGE
	// Gaussian filter Cf Matlab bpass.m
	public static ImageProcessor gaussianFilter(ImageProcessor ip, double lowNoise, double highNoise,
			Params param) {
		ImageProcessor ipR = ip.duplicate();

		// String newLineToFile = "";
		int kernel_size = 2 * 5 * (int) lowNoise + 1;
		int squaredKernel_size = (int) (Math.pow(kernel_size, 2));
		float[] gaussian_kernel = new float[squaredKernel_size];
		float sum = 0;
		for (int i = 0; i < kernel_size; i++) {
			for (int j = 0; j < kernel_size; j++) {
				gaussian_kernel[j + i * kernel_size] = (float) (Math.exp(-Math.pow((-5 * (int) lowNoise + j)
						/ (2.0 * lowNoise), 2.0)))
						* (float) (Math.exp(-Math.pow((-5 * (int) lowNoise + i) / (2.0 * lowNoise), 2.0)));
				sum += gaussian_kernel[j + i * kernel_size];
			}
		}
		for (int i = 0; i < kernel_size; i++) {
			for (int j = 0; j < kernel_size; j++) {
				gaussian_kernel[j + i * kernel_size] /= sum;
			}
		}

		ImageProcessor ip2 = ip.duplicate();
		ip2.convolve(gaussian_kernel, kernel_size, kernel_size);

		kernel_size = 2 * (int) Math.round(highNoise) + 1;
		squaredKernel_size = (int) (Math.pow(kernel_size, 2));
		float[] boxcar_kernel = new float[squaredKernel_size];
		for (int i = 0; i < squaredKernel_size; i++) {
			boxcar_kernel[i] = (float) (1.0D / squaredKernel_size);
		}

		ImageProcessor ip3 = ip.duplicate();
		ip3.convolve(boxcar_kernel, kernel_size, kernel_size);

		if (param.fitting == Params.Fitting.CountingFit)
			ipR = ip2; // TODO: DEAL WITH THAT!!!
		else {
			for (int i = 0; i < ip.getWidth(); i++) {
				for (int j = 0; j < ip.getHeight(); j++) {
					ipR.set(i, j, Math.max(ip2.get(i, j) - ip3.get(i, j), 0));
				}
			}
		}

		// Zero out the values on the edges to signal that they're not useful.
		int lzero = Math.max((int) Math.round(highNoise), 5 * (int) lowNoise);
		for (int i = 0; i < lzero; i++) {
			for (int j = 0; j < ip.getWidth(); j++) {
				ipR.set(j, i, 0);
				ipR.set(j, ip.getHeight() - (i + 1), 0);
			}
			for (int j = 0; j < ip.getHeight(); j++) {
				ipR.set(i, j, 0);
				ipR.set(ip.getWidth() - (i + 1), j, 0);
			}
		}

		return ipR;
	}

	/**
	 * Puts imageprocessor (ROI) into a new imageprocessor of size width x
	 * height y at position (x,y). The image is mirrored around its edges to
	 * avoid wrap around effects of the FFT.
	 */
	public static ImageProcessor tileMirror(ImageProcessor ip, int width, int height, int x, int y) {

		if (x < 0 || x > (width - 1) || y < 0 || y > (height - 1)) {
			IJ.log("Image to be tiled is out of bounds.");
			return null;
		}

		ImageProcessor ipout = ip.createProcessor(width, height);

		ImageProcessor ip2 = ip.crop();
		int w2 = ip2.getWidth();
		int h2 = ip2.getHeight();

		// how many times does ip2 fit into ipout?
		int i1 = (int) Math.ceil(x / (double) w2);
		int i2 = (int) Math.ceil((width - x) / (double) w2);
		int j1 = (int) Math.ceil(y / (double) h2);
		int j2 = (int) Math.ceil((height - y) / (double) h2);

		// tile
		if ((i1 % 2) > 0.5)
			ip2.flipHorizontal();
		if ((j1 % 2) > 0.5)
			ip2.flipVertical();

		for (int i = -i1; i < i2; i += 2) {
			for (int j = -j1; j < j2; j += 2) {
				ipout.insert(ip2, x - i * w2, y - j * h2);
			}
		}

		ip2.flipHorizontal();
		for (int i = -i1 + 1; i < i2; i += 2) {
			for (int j = -j1; j < j2; j += 2) {
				ipout.insert(ip2, x - i * w2, y - j * h2);
			}
		}

		ip2.flipVertical();
		for (int i = -i1 + 1; i < i2; i += 2) {
			for (int j = -j1 + 1; j < j2; j += 2) {
				ipout.insert(ip2, x - i * w2, y - j * h2);
			}
		}

		ip2.flipHorizontal();
		for (int i = -i1; i < i2; i += 2) {
			for (int j = -j1 + 1; j < j2; j += 2) {
				ipout.insert(ip2, x - i * w2, y - j * h2);
			}
		}

		return ipout;
	}

	// OPERATIONS ON CELLS

	public static Polygon[] getCells(String localRootDirName) {
		ImagePlus impWL = openImage(localRootDirName, "_WL.tif");
		impWL.show();

		Polygon[] cells = getCellContours(impWL, 1, localRootDirName, true);

		impWL.close();
		impWL.flush();

		return cells;
	}

	public static Polygon[] getCellContours(ImagePlus rebuilt, int rebuiltScale, String rootDirName,
			boolean getNewIfNone) {

		Polygon[] cells = new Polygon[20];
		int cell = 0;

		if (new File(rootDirName + File.separator + "0_Cell1-Contour.txt").isFile()
				|| new File(rootDirName.substring(0, rootDirName.lastIndexOf(File.separator)).substring(
						0,
						rootDirName.substring(0, rootDirName.lastIndexOf(File.separator))
								.substring(0, rootDirName.lastIndexOf(File.separator))
								.lastIndexOf(File.separator))
						+ File.separator + "0_Cell1-Contour.txt").isFile()) {
			String rootDirNameTemp;
			if (new File(rootDirName + File.separator + "0_Cell1-Contour.txt").isFile())
				rootDirNameTemp = rootDirName;
			else
				rootDirNameTemp = rootDirName.substring(0, rootDirName.lastIndexOf(File.separator))
						.substring(
								0,
								rootDirName.substring(0, rootDirName.lastIndexOf(File.separator))
										.substring(0, rootDirName.lastIndexOf(File.separator))
										.lastIndexOf(File.separator));
			while (new File(rootDirNameTemp + File.separator + "0_Cell" + (cell + 1) + "-Contour.txt")
					.isFile()) {
				ResultsTableMt saveCell = ResultsTableMt.open2(rootDirNameTemp + File.separator + "0_Cell"
						+ (cell + 1) + "-Contour.txt");

				if (getSum(saveCell.getColumnAsDoubles(saveCell.getColumnIndex("Xi")))
						+ getSum(saveCell.getColumnAsDoubles(saveCell.getColumnIndex("Yi"))) % 1 == 0) {
					cells[cell] = new Polygon(doubleToInt(saveCell.getColumnAsDoubles(saveCell
							.getColumnIndex("Xi"))), doubleToInt(saveCell.getColumnAsDoubles(saveCell
							.getColumnIndex("Yi"))), saveCell.getCounter());
				} else {
					int[] xPol = new int[saveCell.getCounter()];
					int[] yPol = new int[xPol.length];
					for (int pol = 0; pol < xPol.length; pol++) {
						xPol[pol] = (int) Math.floor(saveCell.getValue("Xi", pol) + 0.5);
						yPol[pol] = (int) Math.floor(saveCell.getValue("Yi", pol) + 0.5);
					}
					cells[cell] = new Polygon(xPol, yPol, xPol.length);
				}

				cell++;
			}
		} else if (getNewIfNone) {

			ij.gui.WaitForUserDialog waitingGUI = new ij.gui.WaitForUserDialog("Roi selector",
					"Please select a roi in the new window and press 'OK'.");
			waitingGUI.show();
			int cellNumber = 0;
			while (!waitingGUI.escPressed() && cellNumber < cells.length) {
				cells[cellNumber] = rebuilt.getRoi().getPolygon();
				cellNumber++;

				// ImageProcessor cell = Rebuilt.getRoi().getMask();
				// cell.add(Rebuilt.getStatistics().mean);
				Overlay ol = rebuilt.getOverlay();
				if (ol == null)
					ol = new Overlay();

				Roi roiTemp = rebuilt.getRoi();
				roiTemp.setFillColor(new Color(cellNumber * 10, cellNumber * 10, cellNumber * 10));
				roiTemp.setStrokeColor(new Color(cellNumber * 10, cellNumber * 10, cellNumber * 10));
				ol.add(roiTemp);

				TextRoi label1 = new TextRoi(rebuilt.getRoi().getBounds().x
						+ (int) ((double) rebuilt.getRoi().getBounds().width / 2.0), rebuilt.getRoi()
						.getBounds().y + (int) ((double) rebuilt.getRoi().getBounds().height / 2.0), ""
						+ cellNumber, new Font(Font.SANS_SERIF, Font.PLAIN, 8));
				label1.setStrokeColor(Color.BLACK);
				label1.setNonScalable(false);
				ol.add(label1);

				rebuilt.setOverlay(ol);
				rebuilt.setHideOverlay(false);
				rebuilt.updateAndDraw();
				rebuilt.show();

				// rebuilt.setOverlay(Rebuilt.getRoi(), Color.WHITE, 0,
				// Color.GRAY);
				// rebuilt.updateImage();
				rebuilt.deleteRoi();

				waitingGUI = new ij.gui.WaitForUserDialog("Roi selector",
						"Please select a roi in the new window and press 'OK'.");
				waitingGUI.show();
			}

			for (cell = 0; cell < cells.length && cells[cell] != null; cell++) {
				Polygon polygon = cells[cell];

				ResultsTableMt saveCell = new ResultsTableMt();
				for (int i = 0; i < polygon.npoints; i++) {
					saveCell.incrementCounter();
					saveCell.addValue("Xi",
							(int) (((double) polygon.xpoints[i]) / ((double) rebuiltScale) + 0.5));
					saveCell.addValue("Yi",
							(int) (((double) polygon.ypoints[i]) / ((double) rebuiltScale) + 0.5));
				}

				try {
					saveCell.saveAs(rootDirName + File.separator + "0_Cell" + (cell + 1) + "-Contour.txt");
				} catch (IOException e) {
					e.printStackTrace();
				}

				cells[cell] = new Polygon(doubleToInt(saveCell.getColumnAsDoubles(saveCell
						.getColumnIndex("Xi"))), doubleToInt(saveCell.getColumnAsDoubles(saveCell
						.getColumnIndex("Yi"))), saveCell.getCounter());
			}
		}

		Polygon[] retour = new Polygon[cell];
		for (int i = 0; i < cell; i++)
			retour[i] = cells[i];

		return retour;
	}

	// Cf. http://introcs.cs.princeton.edu/java/35purple/Polygon.java.html &&
	// http://mathworld.wolfram.com/PolygonArea.html
	// return area of polygon
	public static double area(Polygon polygon) {
		return Math.abs(signedArea(polygon));
	}

	public static double area(Point[] points) {
		Polygon polygon = new Polygon();
		for (int i = 0; i < points.length; i++) {
			polygon.addPoint((int) (points[i].x + 0.5D), (int) (points[i].y + 0.5D));
		}

		return area(polygon);
	}

	public static double area(Rectangle rectangle) {
		return rectangle.getHeight() * rectangle.getWidth();
	}

	// return signed area of polygon
	public static double signedArea(Polygon polygon) {
		double sum = 0;
		for (int i = 0; i < polygon.npoints - 1; i++) {
			sum += (polygon.xpoints[i] * polygon.ypoints[i + 1])
					- (polygon.ypoints[i] * polygon.xpoints[i + 1]);
		}
		return 0.5 * sum;
	}

	public static boolean intersects(float cx, float cy, float radius, float left, float top, float right,
			float bottom) {
		float closestX = (cx < left ? left : (cx > right ? right : cx));
		float closestY = (cy < top ? top : (cy > bottom ? bottom : cy));
		float dx = closestX - cx;
		float dy = closestY - cy;

		return (dx * dx + dy * dy) <= radius * radius;
	}

	// COMPUTER

	public static class Time {
		private long start, end;

		public Time() {
			start = System.currentTimeMillis();
		}

		public void start() {
			start = System.currentTimeMillis();
		}

		public long stop(boolean print) {
			end = System.currentTimeMillis();
			if (print)
				IJ.log("Time taken by the algorithm: " + (int) ((end - start) / 1000.0) + " s.");
			return end - start;
		}
	}

	public static double[] reverse(double[] vector) {
		double[] retour = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			retour[i] = vector[vector.length - 1 - i];
		}
		return retour;
	}

}
