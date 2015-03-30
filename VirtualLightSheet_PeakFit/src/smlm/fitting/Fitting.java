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

import org.apache.commons.math3.analysis.MultivariateFunction;

import smlm.Params;
import smlm.util.ResultsTableMt;


import ij.IJ;
import ij.process.ImageProcessor;

public abstract class Fitting {

	// Internal parameters
	public boolean[] isArgumentFixed = new boolean[6];
	public double[] fixedArguments = new double[6];
	protected int pixelPrecision;

	private boolean fast = false;

	// Initial parameters
	protected Params param;
	public ImageProcessor roi;
	protected double xMax;
	protected double yMax;
	protected double i0Max;
	protected double backgroundLevel;
	protected double stdBackground;
	protected int roiWidth = 2;

	protected int realFrame;
	protected int cycle;
	protected int maxNumber;

	// Results of the fit
	public ResultsTableMt fit = new ResultsTableMt();
	public double[] results;

	public Fitting(Params param, ImageProcessor ip, double xMax, double yMax, double i0Max,
			double backgroundLevel, double stdBackground, int maxNumber, int realFrame, int cycle) {
		this.param = param;
		int xmin = (int) Math.max(xMax - roiWidth * param.psfSigmaInt, 0);
		int ymin = (int) Math.max(yMax - roiWidth * param.psfSigmaInt, 0);
		int width = (int) Math.min(2 * roiWidth * param.psfSigmaInt + 1, ip.getWidth() - xmin);
		int height = (int) Math.min(2 * roiWidth * param.psfSigmaInt + 1, ip.getHeight() - ymin);
		ip.setRoi(xmin, ymin, width, height);
		try {
			this.roi = ip.crop();
		} catch (ArrayIndexOutOfBoundsException e) {
			IJ.log("xmin: " + xmin + " ymin: " + ymin + " / width: " + width + " height: " + height
					+ "   // " + ip.getWidth() + "  " + ip.getHeight() + " // " + e);
		}
		this.xMax = xMax;
		this.yMax = yMax;
		this.i0Max = i0Max;
		this.backgroundLevel = backgroundLevel;
		this.stdBackground = stdBackground;
		this.maxNumber = maxNumber;
		this.realFrame = realFrame;
		this.cycle = cycle;
	}

	public Fitting(Params param) {
		this.param = param;
	}

	public abstract class LLH implements MultivariateFunction {

		abstract public double value(double[] x);

	}

	private double PSF(double X, double Y, double[] x) {
		return Math.exp(-0.5D * Math.pow((X - x[0]) / x[2], 2.0D) - 0.5D * Math.pow((Y - x[1]) / x[3], 2.0D));
	}

	public double PSFintegration(int X, int Y, double[] x) {
		if (isFast()) {
			return x[4] / (2.0D * Math.PI * x[2] * x[3]) * PSF(((double) X), ((double) Y), x);
		} else {
			double retour = 0.0;
			for (int i = 0; i <= pixelPrecision; i++) {
				for (int j = 0; j <= pixelPrecision; j++) {
					retour += PSF((double) X - 0.5D + (0.5D + (double) i) / (double) pixelPrecision,
							(double) Y - 0.5D + (0.5D + (double) j) / (double) pixelPrecision, x);
				}
			}
			// IJ.write(""+retour);
			return retour * x[4] / (2.0D * Math.PI * x[2] * x[3]) / Math.pow(pixelPrecision + 1.0D, 2.0D);

		}
	}

	public abstract boolean FitThis();

	public abstract ImageProcessor DrawInitialRoi(boolean save);

	public abstract ImageProcessor DrawFit(boolean save);

	public abstract ImageProcessor DrawFitBoolean(boolean save);

	public boolean isFast() {
		return fast;
	}

	public void setFast(boolean fast) {
		this.fast = fast;
	}

}
