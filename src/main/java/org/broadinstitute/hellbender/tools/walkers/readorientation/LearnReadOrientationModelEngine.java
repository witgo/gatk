package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.Histogram;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;


public class LearnReadOrientationModelEngine {

    // When the increase in likelihood falls under this value, we call the algorithm converged
    private final double convergenceThreshold;

    // If the EM does not converge in a few steps we should suspect that something went wrong
    private final int maxEMIterations;

    static final double EPSILON = 1e-3;

    private final String referenceContext;

    private final Nucleotide refAllele;

    private final Histogram<Integer> refHistogram;

    private final List<Histogram<Integer>> altDepthOneHistograms;

    private final List<AltSiteRecord> altDesignMatrix;

    /**
     * N by K matrix of posterior probabilities of latent variable z, where N is the number of alt sites,
     * evaluated at the current estimates of the mixture weights artifactPriors
     */
    private final double[][] altResponsibilities;

    /**
     * Maps a (Depth, Alt allele, F1R2/F2R1) 3-tuple to the responsibilities, where the alt depth is 1
     */
    private final Map<Triple<Integer, Nucleotide, ReadOrientation>, double[]> responsibilitiesOfAltDepth1Sites;

    /**
     * MAX_COVERAGE by K matrix of responsibilities of a ref site (i.e. ALT Depth = 0, ALT F1R2 = 0)
     * for ref sites with coverage 1, 2,..., maxDepth
     */
    private final double[][] refResponsibilities;

    private final int numAltExamples;

    private final int numRefExamples;

    private final int numExamples;

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the altResponsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    private double[] effectiveCounts = new double[F1R2FilterConstants.NUM_STATES];

    // K-dimensional vector of beta-binomial parameters alpha and beta for
    // the conditional distribution over the alt count m given artifact state z
    private final static Map<ArtifactState, Pair<Double, Double>> alleleFractionPseudoCounts = getPseudoCountsForAlleleFraction();

    // K-dimensional vector of beta-binomial parameters alpha and beta for
    // the conditional distribution over the alt F1R2 count x given alt count m
    private final static Map<ArtifactState, Pair<Double, Double>> altF1R2FractionPseudoCounts = getPseudoCountsForAltF1R2Fraction();

    private int numIterations = 0;

    private final Logger logger;


    /**
     * Contract: the reference contexts must be combined with its reverse complements prior to instantiating this class
     */
    public LearnReadOrientationModelEngine(final Histogram<Integer> refHistogram, final List<Histogram<Integer>> altDepthOneHistograms,
                                           final List<AltSiteRecord> altDesignMatrixForContext,
                                           final double convergenceThreshold, final int maxEMIterations, final Logger logger) {
        Utils.nonNull(refHistogram);
        Utils.nonNull(altDepthOneHistograms);
        Utils.nonNull(altDesignMatrixForContext);

        referenceContext = refHistogram.getValueLabel();
        Utils.validateArg(referenceContext.length() == F1R2FilterConstants.REFERENCE_CONTEXT_SIZE,
                String.format("reference context must have length %d but got %s", F1R2FilterConstants.REFERENCE_CONTEXT_SIZE, referenceContext));
        Utils.validateArg(F1R2FilterConstants.CANONICAL_KMERS.contains(referenceContext),
                referenceContext + " is not in the set of canonical kmers");

        this.refHistogram = refHistogram;
        this.altDepthOneHistograms = altDepthOneHistograms;
        this.altDesignMatrix = altDesignMatrixForContext;
        this.numAltExamples = altDesignMatrix.size() + altDepthOneHistograms.stream().mapToInt(h -> (int) h.getSumOfValues()).sum();
        this.numRefExamples = (int) refHistogram.getSumOfValues();
        this.numExamples = numAltExamples + numRefExamples;

        // Responsibilities of ref sites with equal depth are the same so we can compute it for each depth and
        // multiply by the number of counts for that depth
        this.refResponsibilities = new double[F1R2FilterConstants.maxDepth][F1R2FilterConstants.NUM_STATES];

        this.altResponsibilities = new double[altDesignMatrix.size()][F1R2FilterConstants.NUM_STATES];

        // Store responsibilities for each depth and the F1R2/F2R1 of the one alt read
        this.responsibilitiesOfAltDepth1Sites = new HashMap<>();
        this.refAllele = F1R2FilterUtils.getMiddleBase(referenceContext);
        this.convergenceThreshold = convergenceThreshold;
        this.maxEMIterations = maxEMIterations;
        this.logger = logger;
    }

    // Learn the prior probabilities for the artifact states by the EM algorithm
    public ArtifactPrior learnPriorForArtifactStates() {
        boolean converged = false;
        double[] l2distancesOfParameters = new double[maxEMIterations];

        // Initialize the prior for artifact
        double[] statePrior = getFlatPrior(refAllele);
        double[] oldStatePrior = getFlatPrior(refAllele);

        while (!converged && numIterations < maxEMIterations) {
            // Responsibilities are updated by side effect to save space
            takeEstep(statePrior);
            statePrior = takeMstep();
            Utils.validate(Math.abs(MathUtils.sum(statePrior) - 1.0) < EPSILON, "artifactPriors must be normalized");

            // TODO: make sure EM increases the likelihood
            // newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            final double l2Distance = MathArrays.distance(oldStatePrior, statePrior);
            converged = l2Distance < convergenceThreshold;

            l2distancesOfParameters[numIterations] = l2Distance;
            oldStatePrior = Arrays.copyOf(statePrior, F1R2FilterConstants.NUM_STATES);

            numIterations++;
        }

        if (numIterations == maxEMIterations){
            logger.info(String.format("Context %s: with %s ref and %s alt examples, EM failed to converge within %d steps",
                    referenceContext, numRefExamples, numAltExamples, maxEMIterations));
        } else {
            logger.info(String.format("Context %s: with %s ref and %s alt examples, EM converged in %d steps",
                    referenceContext, numRefExamples, numAltExamples, numIterations));
        }

        logger.info("Changes in L2 distance of artifactPriors between iterations: " + ArtifactPrior.doubleArrayToString(l2distancesOfParameters));
        return new ArtifactPrior(referenceContext, statePrior, numExamples, numAltExamples);
    }

    /**
     * Given the current estimates of @{link artifactPriors}, compute the responsibilities, which is
     * the posterior probabilities of artifact states, for each example
     **/
    private void takeEstep(double[] artifactPriors) {
        /**
         * Compute the responsibilities of ref examples.
         * Given that for moderate to high depths we will always have P(HOM REF) = 1, this is largely overkill
         *
         * Ref sites with the same depth have the same alt and alt F1R2 depths (i.e. zero) so avoid repeated computations
         */
        for (int i = 0; i < F1R2FilterConstants.maxDepth; i++) {
            final int depth = i + 1;
            final double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(refAllele, refAllele, 0, 0, depth, artifactPriors);
            refResponsibilities[i] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
            Utils.validate(Math.abs(MathUtils.sum(refResponsibilities[i]) - 1.0) < EPSILON,
                    String.format("ref responsibility for depth = %d added up to %f", depth, MathUtils.sum(refResponsibilities[i])));
        }

        // Compute the responsibilities of alt sites
        for (int n = 0; n < altDesignMatrix.size(); n++) {
            final AltSiteRecord example = altDesignMatrix.get(n);
            final int depth = example.getDepth();
            final int altDepth = example.getAltCount();
            final int altF1R2 = example.getAltF1R2();

            // K-dimensional array of one of the terms that comprises gamma*_{nk}
            double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(refAllele, example.getAltAllele(), altDepth, altF1R2, depth, artifactPriors);

            // Normalize responsibilities here because the M-step uses normalized responsibilities
            altResponsibilities[n] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);

            Utils.validate(Math.abs(MathUtils.sum(altResponsibilities[n]) - 1.0) < EPSILON,
                    String.format("responsibility for %dth example added up to %f", n, MathUtils.sumLog10(altResponsibilities[n])));
        }

        // Compute the responsibilities of sites of ref sites (alt depth = 0) and alt sites with depth=1
        for (int i = 0; i < F1R2FilterConstants.maxDepth; i++){
            final int depth = i+1;
            for (Nucleotide altAllele : F1R2FilterConstants.REGULAR_BASES){
                for (ReadOrientation orientation : ReadOrientation.values()){
                    if (altAllele == refAllele){
                        continue;
                    }

                    final int f1r2Depth = orientation == ReadOrientation.F1R2 ? 1 : 0;

                    double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(
                            refAllele, altAllele, 1, f1r2Depth, depth, artifactPriors);
                    final Triple<Integer, Nucleotide, ReadOrientation> key = createKey(depth, altAllele, orientation);
                    responsibilitiesOfAltDepth1Sites.put(key, MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities));
                }
            }
        }
    }

    /**
     * Given the posterior distributions over the artifact states (ie responsibilities) under the current estimates for the prior,
     * update the prior weights such that they maximize the lower bound on the marginal likelihood P(Data).
     * We do so by maximizing the expectation of the complete data likelihood with respect to the posterior
     * for the artifact states from the E-step
     */
    private double[] takeMstep() {
        // First we compute the effective counts of each state, N_k in the docs. We do this separately over alt and ref sites
        final double[] effectiveAltCountsFromDesignMatrix = GATKProtectedMathUtils.sumArrayFunction(0, altDesignMatrix.size(), n -> altResponsibilities[n]);
        double[] effectiveAltCountsFromHistograms = new double[F1R2FilterConstants.NUM_STATES];

        for (Histogram<Integer> histogram : altDepthOneHistograms){
            final Triple<String, Nucleotide, ReadOrientation> triplet = F1R2FilterUtils.labelToTriplet(histogram.getValueLabel());
            final Nucleotide altAllele = triplet.getMiddle();
            final ReadOrientation orientation = triplet.getRight();


            final double[] effectiveAltCountsFromHistogram = GATKProtectedMathUtils.sumArrayFunction(0, F1R2FilterConstants.maxDepth,
                    i -> MathArrays.scale(histogram.get(i + 1).getValue(), responsibilitiesOfAltDepth1Sites.get(createKey(i+1, altAllele, orientation))));
            effectiveAltCountsFromHistograms = MathArrays.ebeAdd(effectiveAltCountsFromHistograms, effectiveAltCountsFromHistogram);
        }

        final double[] effectiveAltCounts = MathArrays.ebeAdd(effectiveAltCountsFromDesignMatrix, effectiveAltCountsFromHistograms);

        Utils.validate(Math.abs(MathUtils.sum(effectiveAltCounts) - numAltExamples) < EPSILON,
                String.format("effective alt counts must add up to %d but got %f", numAltExamples, MathUtils.sum(effectiveAltCounts)));

        // TODO: at some depth, the responsibilities must be 1 for z = HOM_REF and 0 for everything else, we could probably save some time there
        // Over ref sites, we have a histogram of sites over different depths. At each depth we simply multiply the responsibilities by the number of sites,
        // and sum them over all of depths. Because we cut off the depth histogram at {@code MAX_COVERAGE}, we underestimate the ref effective counts by design
        final double[] effectiveRefCounts = GATKProtectedMathUtils.sumArrayFunction(0, F1R2FilterConstants.maxDepth,
                i -> MathArrays.scale(refHistogram.get(i + 1).getValue(), refResponsibilities[i]));

        Utils.validate(Math.abs(MathUtils.sum(effectiveRefCounts) - numRefExamples) < EPSILON,
                String.format("effective ref counts must add up to %d but got %f", numRefExamples, MathUtils.sum(effectiveRefCounts)));

        effectiveCounts = MathArrays.ebeAdd(effectiveAltCounts, effectiveRefCounts);
        Utils.validate(Math.abs(MathUtils.sum(effectiveCounts) - numExamples) < EPSILON,
                String.format("effective counts must add up to number of examples %d but got %f", numExamples, MathUtils.sum(effectiveCounts)));

        return MathArrays.scale(1.0 / numExamples, effectiveCounts);
    }

    public static double[] computeLog10Responsibilities(final Nucleotide refAllele, final Nucleotide altAllele,
                                                        final int altDepth, final int f1r2AltCount, final int depth,
                                                        final double[] artifactPrior) {
        final double[] log10UnnormalizedResponsibilities = new double[F1R2FilterConstants.NUM_STATES];
        List<ArtifactState> refToRefArtifacts = ArtifactState.getRefToRefArtifacts(refAllele);

        for (ArtifactState state : ArtifactState.values()){
            final int stateIndex = state.ordinal();
            if (refToRefArtifacts.contains(state)) {
                // This state is really just hom ref so give it zero probability and skip
                log10UnnormalizedResponsibilities[stateIndex] = Double.NEGATIVE_INFINITY;
                continue;
            }

            if (ArtifactState.artifactStates.contains(state) && ArtifactState.getAltAlleleOfArtifact(state) != altAllele) {
                // The indicator function is 0
                log10UnnormalizedResponsibilities[stateIndex] = Double.NEGATIVE_INFINITY;
                continue;
            }

            // If we get here, we have a non-artifact state i.e. { germline het, hom ref, hom var, somatic het }
            // or an artifact state whose transitions match the observed alt allele (e.g. alt allele = A, z = F1R2_A, F2R1_A)
            log10UnnormalizedResponsibilities[stateIndex] = computePosterior(altDepth, f1r2AltCount, depth, artifactPrior[stateIndex],
                    alleleFractionPseudoCounts.get(state), altF1R2FractionPseudoCounts.get(state));
        }

        return log10UnnormalizedResponsibilities;
    }


    /**
     * Compute the posterior probability of the state z given data. The caller is responsible for not calling
     * this method on inconsistent states e.g. z = F1R2_C where the reference context is ACT
     */
    private static double computePosterior(final int altDepth, final int altF1R2Depth, final int depth,
                                           final double statePrior, final Pair<Double, Double> afPseudoCounts,
                                           final Pair<Double, Double> f1r2PseudoCounts){
        Utils.validateArg(MathUtils.isAProbability(statePrior), String.format("artifactPriors must be a probability but got %f", statePrior));
        Utils.validateArg(afPseudoCounts.getFirst() > 0 && afPseudoCounts.getSecond() > 0,
                String.format("pseudocounts for allele fraction must be greater than 0 but got %f and %f",
                        afPseudoCounts.getFirst(), afPseudoCounts.getSecond()));
        Utils.validateArg(f1r2PseudoCounts.getFirst() > 0 && f1r2PseudoCounts.getSecond() > 0,
                String.format("pseudocounts for alt F1R2 fraction must be greater than 0 but got %f and %f",
                        f1r2PseudoCounts.getFirst(), f1r2PseudoCounts.getSecond()));

        return Math.log10(statePrior) +
                MathUtils.log10BetaBinomialProbability(altDepth, depth, afPseudoCounts.getFirst(), afPseudoCounts.getSecond()) +
                MathUtils.log10BetaBinomialProbability(altF1R2Depth, altDepth, f1r2PseudoCounts.getFirst(), f1r2PseudoCounts.getSecond());
    }

    /**
     * For each state z, define allele fraction distribution f_z and alt f1r2 fraction theta_z
     * They are both beta binomial distributions and are therefore parameterized by the pseudocounts alpha and beta
     */
    private static Map<ArtifactState, Pair<Double, Double>> getPseudoCountsForAlleleFraction(){
        final Map<ArtifactState, Pair<Double, Double>> alleleFractionPseudoCounts = new HashMap<>(ArtifactState.values().length);

        final double altPseudocount = 1.0;
        final double refPseudocount = 9.0;

        // for home ref sites assume Q35, but maintaining some width
        final double pseudoCountOfHomLikely = 10000.0;
        final double pseudoCountOfHomUnlikely = 3;

        final double balancedPseudocount = 5;

        // These variables define the distribution of allele fraction in a somatic variant. It should be learned from
        // the data in the future. In the meantime, one should tweak this by hand when e.g. applying the read orientation
        // filter on low allele fraction samples such as the blood biopsy
        final double pseudoCountOfSomaticAlt = 2;
        final double pseudoCountOfSomaticRef = 5;

        // The allele fraction distribution, which is not aware of the read orientation, should be the same between
        // F1R2 and F2R1 artifacts
        ArtifactState.getF1R2ArtifactStates().forEach(s -> alleleFractionPseudoCounts.put(s, new Pair<>(altPseudocount, refPseudocount)));
        ArtifactState.getF2R1ArtifactStates().forEach(s -> alleleFractionPseudoCounts.put(s, new Pair<>(altPseudocount, refPseudocount)));


        for (ArtifactState z : ArtifactState.getNonArtifactStates()){
            switch (z) {
                case HOM_REF:
                    alleleFractionPseudoCounts.put(z, new Pair<>(pseudoCountOfHomUnlikely, pseudoCountOfHomLikely));
                    break;
                case GERMLINE_HET:
                    alleleFractionPseudoCounts.put(z, new Pair<>(balancedPseudocount, balancedPseudocount));
                    break;
                case SOMATIC_HET:
                    alleleFractionPseudoCounts.put(z, new Pair<>(pseudoCountOfSomaticAlt, pseudoCountOfSomaticRef));
                    break;
                case HOM_VAR:
                    alleleFractionPseudoCounts.put(z, new Pair<>(pseudoCountOfHomLikely, pseudoCountOfHomUnlikely));
                    break;
            }
        }
        return alleleFractionPseudoCounts;
    }

    private static Map<ArtifactState, Pair<Double, Double>> getPseudoCountsForAltF1R2Fraction(){
        final Map<ArtifactState, Pair<Double, Double>> altF1R2FractionPseudoCounts = new HashMap<>(ArtifactState.values().length);
        final double pseudoCountOfLikelyOutcome = 100.0;
        final double pseudoCountOfRareOutcome = 1.0;
        final double balancedPrior = 10;

        ArtifactState.getF1R2ArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new Pair<>(pseudoCountOfLikelyOutcome, pseudoCountOfRareOutcome)));
        ArtifactState.getF2R1ArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new Pair<>(pseudoCountOfRareOutcome, pseudoCountOfLikelyOutcome)));
        ArtifactState.getNonArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new Pair<>(balancedPrior, balancedPrior)));

        return altF1R2FractionPseudoCounts;
    }

    public double[] getEffectiveCounts(){
        return effectiveCounts;
    }

    public double getEffectiveCounts(ArtifactState state){
        return effectiveCounts[state.ordinal()];
    }

    public static double[] getFlatPrior(final Nucleotide refAllele){
        // We skip the artifact states in which ref transitions to ref e.g. under the ref context AGT, F1R2_G and F2R1_G
        final List<ArtifactState> refToRefStates = ArtifactState.getRefToRefArtifacts(refAllele);

        double[] prior = new double[F1R2FilterConstants.NUM_STATES];

        Arrays.fill(prior, 1.0 / (F1R2FilterConstants.NUM_STATES - refToRefStates.size()));
        for (ArtifactState s : refToRefStates){
            prior[s.ordinal()] = 0;
        }

        return prior;
    }

    private Triple<Integer, Nucleotide, ReadOrientation> createKey(final int depth, final Nucleotide altAllele, final ReadOrientation orientation){
        return new ImmutableTriple<>(depth, altAllele, orientation);
    }




}
