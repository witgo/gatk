package org.broadinstitute.hellbender.utils.solver;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class GradientDescentSolver {

    private final Function<double[], Double> energyFunction;
    private final double[] gradient;
    private final double[] delta;
    private final double[] mt;
    private final double[] vt;
    private final double[] x;
    private double[] lastSolution;
    private final int size;
    private final double learningRate;
    private final double gradientStep;
    private final Logger logger = LogManager.getLogger(this.getClass());

    public GradientDescentSolver(final Function<double[], Double> energyFunction, final int size, final double learningRate, final double gradientStep) {
        this.energyFunction = energyFunction;
        this.size = size;
        this.learningRate = learningRate;
        this.gradientStep = gradientStep;
        gradient = new double[size];
        delta = new double[size];
        mt = new double[size];
        vt = new double[size];
        x = new double[size];
    }

    public double[] getSolution() {
        return lastSolution;
    }

    final String stateString(final double[] x) {
        return "[" + String.join(", ", Arrays.stream(x).boxed().map(String::valueOf).collect(Collectors.toList())) + "]";
    }

    public double solve(final double[] x0, final int numSteps, final double maxStepSize, final int loggingInterval) {
        copyStates(x0, x);
        final double[] minX = new double[size];
        double currentEnergy = computeEnergy(x);
        if (loggingInterval > 0) {
            logger.info("Solving " + size + " parameters over " + numSteps + " steps; f0 = " + currentEnergy);
        }
        double minEnergy = currentEnergy;
        int tMin = 0;
        copyStates(x, minX);
        double testGradient;
        for (int step = 0; step < numSteps; step++) {
            testGradient = computeGradient();
            computeDelta(step, maxStepSize);
            step();
            if (loggingInterval > 0 && step % loggingInterval == 0) {
                logger.info("Iteration " + step + " : f = " + currentEnergy + ", fmin = " + minEnergy + ", gradient = " + testGradient / Math.sqrt(size*gradientStep*gradientStep));
                logger.info("\tx = " + stateString(x));
            }
        }
        if (loggingInterval > 0) {
            logger.info("Solved after " + numSteps + " : f = " + currentEnergy + ", fmin = " + minEnergy + ", tmin = " + tMin);
            logger.info("\tx_min = " + stateString(minX));
        }
        lastSolution = minX;
        return minEnergy;
    }

    private void computeDelta(final int stepNumber, final double maxStepSize) {
        for (int i = 0; i < gradient.length; i++) {

            final double beta1 = 0.9; //0.99;
            final double beta2 = 0.999; //0.9999;

            //Basic
            //deltaX[i] = gradientX[i] * learningRate;
            //if (Math.abs(deltaX[i]) > maxStepSize) {
            //    deltaX[i] = maxStepSize * Math.signum(deltaX[i]);
            //}

            //Simple momentum
            //delta[i] = mt[i] * 0.99 + gradient[i] * learningRate;
            //if (Math.abs(delta[i]) > maxStepSize) {
            //    delta[i] = maxStepSize * Math.signum(delta[i]);
            //}
            //mt[i] = delta[i];

            //ADAM
            if (stepNumber == 0) {
                mt[i] = gradient[i];
                vt[i] = gradient[i] * gradient[i];
            }
            mt[i] = (beta1 * mt[i] +  (1 - beta1) * gradient[i]); // / (1 - Math.pow(beta1, numIter + 1));
            vt[i] = (beta2 * vt[i] +  (1 - beta2) * gradient[i] * gradient[i]); // / (1 - Math.pow(beta2, numIter + 1));
            delta[i] = - mt[i] * learningRate / (Math.sqrt(vt[i]) + 1e-8);

            //ADAMAX
            //mtX[i] = (beta1 * mtX[i] +  (1 - beta1) * gradientX[i]) / (1 - Math.pow(beta1, numIter + 1));
            //vtX[i] = Math.max(beta2 * vtX[i], Math.abs(gradientX[i]));
            //deltaX[i] = mtX[i] * learningRate / (vtX[i] + 1e-8);

        }
    }

    private double step() {
        double deltaTotal = 0;
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] + delta[i];
            deltaTotal += delta[i]*delta[i];
        }
        return Math.sqrt(deltaTotal);
    }

    private void copyStates(final double[] fromX, final double[] toX) {
        IntStream.range(0, size).forEach(i -> toX[i] = fromX[i]);
    }

    private double computeEnergy(final double[] x) {
        return energyFunction.apply(x);
    }

    private double computeGradient() {
        double total = 0;
        for (int i = 0; i < x.length; i++) {
            double xTemp = x[i];
            x[i] += gradientStep;
            final double testLogPosterior1 = computeEnergy(x);
            x[i] = xTemp - gradientStep;
            final double testLogPosterior2 = computeEnergy(x);
            x[i] = xTemp;
            gradient[i] = (testLogPosterior1 - testLogPosterior2) / (2 * gradientStep);
            total += gradient[i]*gradient[i];
        }
        return Math.sqrt(total);
    }
}