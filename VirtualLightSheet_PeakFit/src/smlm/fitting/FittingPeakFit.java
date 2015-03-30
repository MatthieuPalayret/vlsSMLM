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

package smlm.fitting;

import java.awt.Rectangle;
import java.io.File;

import smlm.util.SRutil;
import smlm.Params;
import smlm.util.ResultsTableMt;


import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitJob;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitCriteria;
import gdsc.smlm.fitting.FitFunction;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.fitting.function.CCDCameraNoiseModel;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.utils.NoiseEstimator.Method;

public class FittingPeakFit {

	boolean fitted = false;

	Params param;
	FitEngineConfiguration config = null;
	MemoryPeakResults results = new MemoryPeakResults();

	static final public int USUAL = 1;
	static final public int VLS = 2;
	static final public int FIBRILS = 3;

	public FittingPeakFit() {
		super();
	}

	public void setConfig(Params param, int methods) {
		this.param = param;

		FitConfiguration fitConf = new FitConfiguration();

		fitConf.setNmPerPixel(param.pixelSize);
		fitConf.setGain((float) (param.emGain / Params.ConversionGain[param.camera][param.gain]));
		fitConf.setInitialPeakStdDev((float) param.psfSigma / 2.0f);
		fitConf.setInitialAngle(0f);

		fitConf.setFitSolver(FitSolver.LVM);
		fitConf.setFitFunction(FitFunction.FREE);
		if (methods == VLS)
			fitConf.setFitFunction(FitFunction.FREE);
		if (methods == FIBRILS)
			fitConf.setFitFunction(FitFunction.FREE);
		fitConf.setFitCriteria(FitCriteria.LEAST_SQUARED_ERROR);
		fitConf.setNoiseModel(new CCDCameraNoiseModel(805, 85, true));// param.Readnoise,
																		// param.Offset,
																		// true));
		fitConf.setSignificantDigits(5);
		fitConf.setMaxIterations(20);

		fitConf.setCoordinateShiftFactor(1.5f);
		if (methods == VLS)
			fitConf.setCoordinateShiftFactor(3f);
		fitConf.setSignalStrength(20f);
		if (methods == VLS)
			fitConf.setSignalStrength(20f);
		if (methods == FIBRILS)
			fitConf.setSignalStrength(5f);
		fitConf.setWidthFactor(5f);
		if (methods == VLS)
			fitConf.setWidthFactor(30f);
		if (methods == FIBRILS)
			fitConf.setWidthFactor(30f);
		fitConf.setPrecisionThreshold(0f);

		fitConf.setDuplicateDistance(0.5f);

		config = new FitEngineConfiguration(fitConf);
		config.setSmooth(0.5);
		if (param.pixelSize > 150)
			config.setSmooth(-1);
		config.setSmooth2(3);
		if (methods == VLS)
			config.setSmooth2(5);
		if (methods == FIBRILS)
			config.setSmooth2(5);
		config.setSearch(3);
		if (methods == VLS)
			config.setSearch(8);
		if (methods == FIBRILS)
			config.setSearch(5);

		config.setFailuresLimit(10);
		config.setIncludeNeighbours(true);
		config.setNeighbourHeightThreshold(0.3);
		config.setResidualsThreshold(1);

		config.setNoiseMethod(Method.QuickResidualsLeastMeanOfSquares);
		config.initialiseState();
	}

	public void loadConfig(String path) {
		GlobalSettings settings = SettingsManager.loadSettings(path);
		config = settings.getFitEngineConfiguration();
		config.initialiseState();
	}

	public void getNewConfig() {
		config = SettingsManager.loadSettings().getFitEngineConfiguration();
		config.initialiseState();
	}

	public void saveConfig(String path) {
		SettingsManager.saveFitEngineConfiguration(config, path);
	}

	public MemoryPeakResults fitImage(ImagePlus imp) {

		// Check if the image is open
		ImageSource source = new IJImageSource(imp);
		if (!source.open())
			return results;

		// Check the configuration
		if (config == null)
			return results;

		// Create a fit engine
		FitEngine engine = new FitEngine(config, results, Prefs.getThreads(), FitQueue.BLOCKING, 0);

		results.begin();

		boolean shutdown = false; // Flag to allow escape to shutdown the
									// fitting

		// Show fitting progress
		int slice = 0;
		int totalFrames = imp.getStackSize();
		final int step = (totalFrames > 400) ? totalFrames / 200 : 2;

		// Extract a region if necessary
		Rectangle bounds = (imp.getRoi() != null) ? imp.getRoi().getBounds() : null;

		while (!shutdown) {
			float[] data = source.next(bounds);
			if (data == null)
				break;

			if (++slice % step == 0) {
				IJ.showProgress(slice, totalFrames);
				IJ.showStatus("Slice: " + slice + " / " + totalFrames);
			}

			engine.run(new FitJob(slice, data, bounds));

			if (IJ.escapePressed())
				shutdown = true;
		}

		engine.end(shutdown);
		results.end();
		fitted = true;

		saveConfig(param.fileDirName + File.separator + "gdsc.smlm.settings.xml");

		return results;
	}

	public ResultsTableMt getResults() {
		ResultsTableMt rtFit = new ResultsTableMt();

		for (PeakResult peak : results.getResults()) {
			rtFit.incrementCounter();
			rtFit.addValue(ResultsTableMt.FRAME, peak.getEndFrame());
			rtFit.addValue(ResultsTableMt.X0, peak.getXPosition());
			rtFit.addValue(ResultsTableMt.Y0, peak.getYPosition());
			rtFit.addValue(ResultsTableMt.I0, peak.getSignal());
			rtFit.addValue(ResultsTableMt.SIGMAX, peak.getXWidth());
			rtFit.addValue(ResultsTableMt.SIGMAY, peak.getYWidth());
			rtFit.addValue(ResultsTableMt.NOISE, peak.noise);
			rtFit.addValue(ResultsTableMt.OFFSET, peak.getBackground());
			rtFit.addValue(ResultsTableMt.IS_FITTED, peak.error);
			if (param != null)
				rtFit.addValue(ResultsTableMt.SIGMArebuilt, SRutil.sigmaRebuilt(param, rtFit, rtFit.getCounter() - 1));
			rtFit.addValue(ResultsTableMt.AMPLITUDE, SRutil.getAmplitude(rtFit, rtFit.getCounter() - 1));
			if (config.getFitConfiguration().getFitFunction() == FitFunction.FREE)
				rtFit.addValue(ResultsTableMt.THETA, peak.getAngle());
		}

		return rtFit;
	}

}
