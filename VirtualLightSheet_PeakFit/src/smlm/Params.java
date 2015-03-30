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

package smlm;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.Date;

import smlm.util.SRutil;
import smlm.util.ResultsTableMt;

import ij.IJ;
import ij.gui.GenericDialog;


public class Params {

	public final static double version = 1.01;
	public String rootDirName;
	public String fileName;
	public String fileDirName;
	public Rectangle roi;
	public Polygon[] cells;

	public boolean hasBeenInitialised = false;

	public double pixelSize = 110; // nm
	public double pixelSizeZ = 20; // nm
	public double na = 1.49;
	public double lambda = 561; // nm

	public double psfSigma = 1.323 / (pixelSize * 2.0 * Math.PI * na / lambda) * 2.0;
	public int psfSigmaInt = (int) Math.max(Math.ceil(psfSigma), 1);

	public enum Fitting {
		NoFit(0), FullGaussianFit(1), FastGaussianFit(2), CentroidFit(3), CountingFit(4), CountingBinaryFit(5)
		/* Not working correctly */, PeakFit(6);

		private int code;

		private Fitting(int code) {
			this.code = code;
		}

		int getCode() {
			return code;
		};
	}

	public Fitting fitting = Fitting.NoFit;

	public boolean photonCountMode = false;
	public double emGain = 250;
	public int gain = 0;
	public int camera = 0;
	public final static String[] CameraNames = { "Arabidopsis", "Gerbera", "Agave_1", "Agave_2", "Aloe_Vera" };
	public static double[][] ConversionGain = new double[][] { new double[] { 11.5, 6.8, 3.1 }, // Arabidopsis
			new double[] { 4.5, 2.7 }, // Gerbera
			new double[] { 8.4, 4.8 }, // Agave 1
			new double[] { 4.2, 2.6 }, // Agave 2
			new double[] { 4.4, 2.5 } // Aloe Vera
	};
	public double offset = 500;

	public double sigmaThreshold = 1.7; // Threshold applied on both SigmaX and
										// SigmaY (in pixel) (e.g. 1.7)
	public double i0Threshold = 200; // Threshold applied on I0 (in photons)
										// (e.g. 200)
	public double ellipticityThreshold = -1; // Threshold applied on the
												// ellipticity, i.e.
												// max(SigmaX/SigmaY,
												// SigmaY/SigmaX) (e.g. 1.7)

	public enum Rebuilding {
		NoRebuilt(0), ClassicalGaussianPSF(1), IdenticalGaussianPSF(2), IdenticalCirclePSF(3), ColouringCirclePSF(
				4), BooleanReplacePSF(5), BooleanAddPSF(6), ColorCodedPerGroup(7), DiffractionLimited(8);

		private int code;

		private Rebuilding(int code) {
			this.code = code;
		}

		int getCode() {
			return code;
		};
	}

	public Rebuilding rebuilding = Rebuilding.NoRebuilt;
	public int rebuiltScale = 10;
	public boolean theoreticalRebuilt = true; // SigmaRebuilt calculated through
												// the theoretical or
												// experimental equation
	public int fixedPSFSizeRebuilt = 20; // in nm

	public boolean scaleBar = true;

	public enum Grouping {
		NoGrouping(0), ClassicalGrouping(1), DBSCAN(2), GroupInTraj(3);

		private int code;

		private Grouping(int code) {
			this.code = code;
		}

		int getCode() {
			return code;
		};
	}

	public Grouping group = Grouping.NoGrouping;

	public double groupXYnm = 150.0;
	public double groupTms = 52;
	public int exposure = 50;

	public boolean debug = false;
	public boolean checkMovie = false;
	public boolean preBleaching = true;

	public Params(String rootDirName, String fileName) {
		this(rootDirName);
		this.fileName = fileName;
		fileDirName = rootDirName + File.separator + fileName.substring(0, fileName.lastIndexOf(".tif"));
	}

	public Params(String rootDirName) {
		this.rootDirName = rootDirName;
		this.fileDirName = rootDirName;
	}

	public void getNewParameters() {
		getNewParameters(true, true, true, true, true, false);
	}

	public void getNewParameters(boolean askForParameters, boolean askForFitting, boolean askForRebuilding,
			boolean askForMore, boolean askForGroup, boolean askForGain) {

		if (SRutil.getFileName(rootDirName, "metaData.txt") != null) {
			try {
				ResultsTableMt metaData = ResultsTableMt.open(rootDirName + File.separator
						+ SRutil.getFileName(rootDirName, "metadata.txt"));

				photonCountMode = metaData.getValue("QuantView", 0) == 1;
				emGain = metaData.getValue("EMGain", 0);
				gain = (int) metaData.getValue("Gain", 0);
				camera = (int) metaData.getValue("Camera", 0);

				askForGain = false;
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else
			askForGain = true;

		// Asking the parameters to the user
		GenericDialog gui = new GenericDialog("Super-Resolution Algorithm");

		if (askForParameters) {
			gui.addMessage("Parameters:");
			gui.addNumericField("Size_of_a_pixel (nm):", 110, 0);
			gui.addNumericField("Size_of_a_z-step (nm) (if applicable):", 100, 0);
			gui.addNumericField("NA:", 1.49, 2);
			gui.addNumericField("Emission peak (nm):", 561, 0);
		}

		if (askForFitting) {
			gui.addMessage("Fitting:");
			String[] Methods = { "No Fit", "Full Gaussian Fit", "Fast Gaussian Fit", "Centroid Fit",
					"Counting Fit", "(CountingBinary Fit)", "PeakFit" };
			gui.addChoice("Method:", Methods, "Full");
		}

		if (askForRebuilding) {
			gui.addMessage("Rebuilding picture:");
			String[] rebuildingMethods = { "No rebuilding", "Classical Gaussian PSF",
					"Identical Gaussian PSF", "Identical Circle PSF", "Colouring Circle PSF",
					"Boolean Replace PSF", "Boolean Add PSF", "Colouring per group", "Diffraction limited" };
			gui.addChoice("Size_of_the_rebuilt_PSFs:", rebuildingMethods,
					"Depending on their intensity and fit");
			gui.addNumericField("Sigma threshold (pix)", 1.7, 1);
			gui.addNumericField("Ellipticity threshold", -1, 1);
			gui.addNumericField("Intensity threshold (photons)", 200, 0);
			gui.addNumericField("Number of sub-pixels:", 20, 0);
			gui.addCheckbox("Display 1um scale bar:", true);
		}

		if (askForMore) {
			gui.addCheckbox("QuantView", photonCountMode);
			String[] groupingMethods = { "No grouping", "Classical grouping", "DBSCAN",
					"Group in trajectories" };
			gui.addChoice("Group_Method:", groupingMethods, "No grouping");
			gui.addCheckbox("Rebuilt the method movie", false);
			gui.addCheckbox("Pre-bleaching", true);
		}

		if (askForGroup) {
			gui.addMessage("Grouping parameters:");
			gui.addNumericField("Radius (nm):", 150, 1);
			gui.addNumericField("Time window (ms):", 52, 0);
			gui.addNumericField("Exposure (ms/frame):", 50, 0);
		}

		if (askForGain) {
			gui.addMessage("Gain used for imaging:");
			gui.addChoice("Camera used:", CameraNames, CameraNames[camera]);
			String[] gains = { "1", "2", "3" };
			gui.addChoice("Gain", gains, "" + gain);
			gui.addNumericField("EMGain:", emGain, 0);
		}

		gui.showDialog();

		if (askForParameters) {
			pixelSize = gui.getNextNumber();
			pixelSizeZ = gui.getNextNumber();
			na = gui.getNextNumber();
			lambda = gui.getNextNumber();
		}

		if (askForFitting)
			fitting = Fitting.values()[gui.getNextChoiceIndex()];

		if (askForRebuilding) {
			rebuilding = Rebuilding.values()[gui.getNextChoiceIndex()];
			sigmaThreshold = gui.getNextNumber();
			ellipticityThreshold = gui.getNextNumber();
			i0Threshold = gui.getNextNumber();
			rebuiltScale = (int) gui.getNextNumber();
			scaleBar = gui.getNextBoolean();
		}

		if (askForMore) {
			photonCountMode = gui.getNextBoolean();
			group = Grouping.values()[gui.getNextChoiceIndex()];
			checkMovie = gui.getNextBoolean();
			preBleaching = gui.getNextBoolean();
		}

		if (askForGroup) {
			groupXYnm = gui.getNextNumber();
			groupTms = gui.getNextNumber();
			exposure = (int) gui.getNextNumber();
		}

		if (askForGain) {
			camera = gui.getNextChoiceIndex();
			gain = gui.getNextChoiceIndex();
			emGain = gui.getNextNumber();
		}

		save();

		update();

		hasBeenInitialised = true;

	}

	public void save() {
		save(fileDirName + File.separator + "Parameters.txt");
	}

	public void save(String filePath) {

		ResultsTableMt params = new ResultsTableMt();
		params.incrementCounter();
		params.addValue("PixelSize", pixelSize);
		params.addValue("PixelSizeZ", pixelSizeZ);
		params.addValue("NA", na);
		params.addValue("lambda", lambda);

		params.addValue("Fitting", fitting.getCode());

		params.addValue("Rebuilding", rebuilding.getCode());
		params.addValue("SigmaThreshold", sigmaThreshold);
		params.addValue("EllipticityThreshold", ellipticityThreshold);
		params.addValue("I0Threshold", i0Threshold);
		params.addValue("PhotonCountMode", photonCountMode ? 1.0 : 0.0);
		params.addValue("RebuiltScale", rebuiltScale);
		params.addValue("ScaleBar", scaleBar ? 1.0 : 0.0);

		params.addValue("Group", group.getCode());
		params.addValue("CheckMovie", checkMovie ? 1.0 : 0.0);
		params.addValue("PreBleaching", preBleaching ? 1.0 : 0.0);

		params.addValue("GroupXYnm", groupXYnm);
		params.addValue("GroupTms", groupTms);
		params.addValue("exposure", exposure);

		params.addValue("Camera", camera);
		params.addValue("Gain", gain);
		params.addValue("EMGain", emGain);

		params.addValue("Date", new Date().getTime());
		params.addValue("Version", version);

		try {
			params.saveAs(filePath);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void load() {
		load(fileDirName + File.separator + "Parameters.txt");
	}

	public void load(String filePath) {
		try {
			ResultsTableMt params = ResultsTableMt.open(filePath);

			if (params.getColumnIndex("Version") != -1 && params.getValue("Version", 0) >= version) {
				pixelSize = params.getValue("PixelSize", 0);
				pixelSizeZ = params.getValue("PixelSizeZ", 0);
				na = params.getValue("NA", 0);
				lambda = params.getValue("lambda", 0);

				fitting = Fitting.values()[(int) params.getValue("Fitting", 0)];

				rebuilding = Rebuilding.values()[(int) params.getValue("Rebuilding", 0)];
				sigmaThreshold = params.getValue("SigmaThreshold", 0);
				ellipticityThreshold = params.getValue("EllipticityThreshold", 0);
				i0Threshold = params.getValue("I0Threshold", 0);
				photonCountMode = params.getValue("PhotonCountMode", 0) == 1;
				rebuiltScale = (int) params.getValue("RebuiltScale", 0);
				scaleBar = params.getValue("ScaleBar", 0) == 1.0;

				group = Grouping.values()[(int) params.getValue("Group", 0)];
				checkMovie = params.getValue("CheckMovie", 0) == 1;
				preBleaching = params.getValue("PreBleaching", 0) == 1;

				groupXYnm = params.getValue("GroupXYnm", 0);
				groupTms = params.getValue("GroupTms", 0);
				exposure = (int) params.getValue("exposure", 0);

				camera = (int) params.getValue("Camera", 0);
				gain = (int) params.getValue("Gain", 0);
				emGain = params.getValue("EMGain", 0);

				update();

				hasBeenInitialised = true;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		if (SRutil.getFileName(rootDirName, "metaData.txt") != null) {
			try {
				ResultsTableMt MetaData = (ResultsTableMt) ResultsTableMt.open(rootDirName + File.separator
						+ SRutil.getFileName(rootDirName, "metadata.txt"));

				photonCountMode = MetaData.getValue("QuantView", 0) == 1;
				camera = (int) MetaData.getValue("Camera", 0);
				emGain = MetaData.getValue("EMGain", 0);
				gain = (int) MetaData.getValue("Gain", 0);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void update() {
		psfSigma = 1.323 / (pixelSize * 2.0 * Math.PI * na / lambda) * 2.0;
		psfSigmaInt = (int) Math.max(Math.ceil(psfSigma), 1);
	}

	@Override
	public Params clone() {
		Params retour = null;
		try {
			retour = (Params) super.clone();

			if (roi != null)
				retour.roi = new Rectangle(roi.x, roi.y, roi.width, roi.height);
			if (cells != null)
				retour.cells = cells.clone();

			return retour;
		} catch (CloneNotSupportedException e) {
			IJ.log(e.getLocalizedMessage());
			return null;
		}
	}
}
