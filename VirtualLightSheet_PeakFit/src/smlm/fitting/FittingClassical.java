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
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import smlm.util.SRutil;
import smlm.Params;
import smlm.util.ResultsTableMt;

import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;

public class FittingClassical extends Fitting {

	boolean LLHorLSE = true;

	public class LLH implements MultivariateFunction {

		public double value(double[] x) {
			for (int i = 0; i < isArgumentFixed.length; i++) {
				if (isArgumentFixed[i]) {
					x[i] = fixedArguments[i];
				}
			}

			double LLHValue = 0;
			double temp = 0;
			for (int i = 0; i < roi.getWidth(); i++) {
				for (int j = 0; j < roi.getHeight(); j++) {
					temp = PSFintegration(i, j, x);
					if (LLHorLSE)
						LLHValue += ((double) roi.get(i, j)) * Math.log(temp + x[5]) - (temp + x[5]);
					else
						LLHValue += Math.pow(((double) roi.get(i, j)) - (temp + x[5]), 2.0D);
				}
			}
			return LLHValue;
		}
	}

	public FittingClassical(Params param, ImageProcessor ip, double xMax, double yMax, double i0Max,
			double backgroundLevel, double stdBackground, int maxNumber, int realFrame, int cycle) {
		super(param, ip, xMax, yMax, i0Max, backgroundLevel, stdBackground, maxNumber, realFrame, cycle);
		pixelPrecision = 3;
	}

	public FittingClassical(Params param) {
		super(param);
	}

	@SuppressWarnings("unused")
	@Override
	public boolean FitThis() {

		// Initial estimates (your initial x)
		// Calculating the centroid
		ImageStatistics stat = roi.getStatistics();
		double x0 = stat.xCenterOfMass;
		double y0 = stat.yCenterOfMass;
		double[] start = { x0, y0, param.psfSigma, param.psfSigma, i0Max, backgroundLevel };

		PointValuePair solutionMult = null;
		MultivariateFunctionMappingAdapter fitFunc = null;

		if (param.fitting != Params.Fitting.CentroidFit) {
			// initial step sizes (take a good guess)
			double[] step = { 1, 1, 0.1, 0.1, stdBackground, stdBackground };
			// convergence tolerance
			double ftol = 0.0001;// 0.000001;
			pixelPrecision = 3;// 5;

			if (param.fitting == Params.Fitting.FastGaussianFit) {
				ftol = 0.1;
				pixelPrecision = 3;
			}

			SimplexOptimizer Fit = new SimplexOptimizer(ftol, ftol * ftol);
			double[] low = new double[] { 0, 0, 0, 0, 0, 0 };
			double[] up = new double[] { roi.getWidth(), roi.getHeight(), Double.POSITIVE_INFINITY,
					Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
			fitFunc = new MultivariateFunctionMappingAdapter(new LLH(), low, up);

			// maximal number of iterations
			int maxIter = 5000;
			// Nelder and Mead maximisation procedure
			// x0 e [0, xmax]
			// Fit.addConstraint(0, -1, 0);
			// Fit.addConstraint(0, 1, roi.getWidth());
			// y0 e [0, ymax]
			// Fit.addConstraint(1, -1, 0);
			// Fit.addConstraint(1, 1, roi.getHeight());
			/*
			 * // sigmax e [PSFSigma/3, 3*PSFSigma] Fit.addConstraint(2, -1,
			 * PSFSigmaInt-PSFSigmaInt/2); Fit.addConstraint(2, 1,
			 * 2*PSFSigmaInt); // sigmay e [PSFSigma/3, 3*PSFSigma]
			 * Fit.addConstraint(3, -1, PSFSigmaInt-PSFSigmaInt/2);
			 * Fit.addConstraint(3, 1, 2*PSFSigmaInt); // I0 e [StdBackground,
			 * 50*Intensities[i]] Fit.addConstraint(4, -1, StdBackground);
			 * Fit.addConstraint(4, 1, 50*XY[3][i]); // PoissonNoise e
			 * [BackgroundLevel/3, 3*BackgroundLevel] Fit.addConstraint(5, -1,
			 * BackgroundLevel/3); Fit.addConstraint(5, 1, 3*BackgroundLevel);
			 */

			solutionMult = Fit.optimize(new MaxEval(maxIter), new ObjectiveFunction(fitFunc),
					GoalType.MAXIMIZE, new InitialGuess(fitFunc.boundedToUnbounded(start)),
					new MultiDirectionalSimplex(step));
		}

		// Result of minimisation
		// Save the fit results
		this.fit.incrementCounter();
		if (param.fitting == Params.Fitting.CentroidFit) {
			results = start;
		} else {
			results = fitFunc.unboundedToBounded(solutionMult.getPoint());
		}
		for (int i = 0; i < isArgumentFixed.length; i++) {
			if (isArgumentFixed[i]) {
				results[i] = fixedArguments[i];
			}
		}

		if (cycle != -1)
			this.fit.addValue(ResultsTableMt.CYCLE, cycle);
		this.fit.addValue(ResultsTableMt.FRAME, realFrame);
		this.fit.addValue(ResultsTableMt.X0, results[0] + xMax - ((double) roiWidth * param.psfSigmaInt));
		this.fit.addValue(ResultsTableMt.Y0, results[1] + yMax - ((double) roiWidth * param.psfSigmaInt));
		this.fit.addValue(ResultsTableMt.SIGMAX, results[2]);
		this.fit.addValue(ResultsTableMt.SIGMAY, results[3]);
		this.fit.addValue(ResultsTableMt.I0, results[4]);
		this.fit.addValue(ResultsTableMt.NOISE, results[5]);

		if (param.fitting == Params.Fitting.CentroidFit || true) {
			this.fit.addValue(ResultsTableMt.IS_FITTED, 1);
			if (param.fitting != Params.Fitting.CentroidFit) {
				this.fit.addValue("MinFit", solutionMult.getValue());
			}
		} else
			this.fit.addValue(ResultsTableMt.IS_FITTED, 0);

		// Save results
		if (param.debug) {
			DrawInitialRoi(true);
			DrawFit(true);
		}
		return true;
	}

	@Override
	public ImageProcessor DrawInitialRoi(boolean save) {
		if (save)
			SRutil.saveTiff(new ImagePlus("", roi), param.fileDirName + "/Temp/Frame_" + realFrame
					+ "_CountingToFit_" + maxNumber + ".tif", true);
		return roi;
	}

	@Override
	public ImageProcessor DrawFit(boolean save) {
		ImageProcessor roiFitted = new ShortProcessor(roi.getWidth(), roi.getHeight());

		for (int i = 0; i < roi.getWidth(); i++) {
			for (int j = 0; j < roi.getHeight(); j++) {
				roiFitted.set(i, j, (int) PSFintegration(i, j, results));
			}
		}

		SRutil.saveTiff(new ImagePlus("", roiFitted), param.fileDirName + "/Temp/Frame_" + realFrame
				+ "_CountingFitted_" + maxNumber + ".tif", true);
		return roiFitted;
	}

	@Override
	public ImageProcessor DrawFitBoolean(boolean save) {
		return null;
	}

}
