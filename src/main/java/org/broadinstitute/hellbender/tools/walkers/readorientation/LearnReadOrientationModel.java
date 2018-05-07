package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Learn the prior probability of read orientation artifact from the output of {@link CollectF1R2Counts}
 * Details of the model may be found in docs/mutect/mutect.pdf.
 */
@CommandLineProgramProperties(
        summary = "Collect counts of F1R2 reads at each locus of a sam/bam/cram",
        oneLineSummary = "Collect counts of F1R2 reads at each locus of a sam/bam/cram",
        programGroup = CoverageAnalysisProgramGroup.class // TODO: check that this is correct
)
public class LearnReadOrientationModel extends CommandLineProgram {
    public static final double DEFAULT_CONVERGENCE_THRESHOLD = 1e-4;
    public static final int DEFAULT_MAX_ITERATIONS = 20;

    public static final String EM_CONVERGENCE_THRESHOLD_LONG_NAME = "convergence-threshold";
    public static final String MAX_EM_ITERATIONS_LONG_NAME = "num-em-iterations";

    @Argument(fullName = CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, doc = "histograms of depths over ref sites for each reference context")
    private File refHistogramTable;

    @Argument(fullName = CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME,  doc = "a table of F1R2 and depth counts")
    private File altDataTable;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "table of artifact priors")
    private File output;

    @Argument(fullName = CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME, doc = "histograms of depth 1 alt sites", optional = true)
    private File altHistogramTable = null;

    @Argument(fullName = EM_CONVERGENCE_THRESHOLD_LONG_NAME, doc = "Stop the EM when the distance between parameters between iterations falls below this value", optional = true)
    private double converagenceThreshold = DEFAULT_CONVERGENCE_THRESHOLD;

    @Argument(fullName = MAX_EM_ITERATIONS_LONG_NAME, doc = "give up on EM after this many iterations", optional = true)
    private int maxEMIterations = DEFAULT_MAX_ITERATIONS;

    List<Histogram<Integer>> refHistograms;
    List<Histogram<Integer>> altHistograms;

    @Override
    protected void onStartup(){
        final MetricsFile<?, Integer> referenceSiteMetrics = readMetricsFile(refHistogramTable);
        refHistograms = referenceSiteMetrics.getAllHistograms();

        if (altHistogramTable != null) {
            final MetricsFile<?, Integer> altSiteMetrics = readMetricsFile(altHistogramTable);
            altHistograms = altSiteMetrics.getAllHistograms();
        } else {
            altHistograms = Collections.emptyList();
        }
    }

    @Override
    public Object doWork(){
        final int defaultInitialListSize = 1_000_000;

        final Map<String, List<AltSiteRecord>> altDesignMatrixByContext =
                AltSiteRecord.readAltSiteRecords(altDataTable, defaultInitialListSize).stream()
                        .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));

        // Since AGT F1R2 is equivalent to ACT F2R1 (in the sense that the order of bases in the original molecule on which
        // the artifact befell and what the base changed to is the same)
        final List<ArtifactPrior> artifactPriorAcrossContexts = new ArrayList<>((int) Math.pow(F1R2FilterConstants.REGULAR_BASES.size(), F1R2FilterConstants.REFERENCE_CONTEXT_SIZE)/2);

        // TODO: extract a method to account for reverse complement, and create a test
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final String reverseComplement = SequenceUtil.reverseComplement(refContext);

            final Histogram<Integer> refHistogram = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(refContext))
                    .findFirst().orElseGet(() -> F1R2FilterUtils.createRefHistogram(refContext));
            final Histogram<Integer> refHistogramRevComp = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(reverseComplement))
                    .findFirst().orElseGet(() -> F1R2FilterUtils.createRefHistogram(reverseComplement));

            final List<Histogram<Integer>> altDepthOneHistogramsForContext = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(refContext))
                    .collect(Collectors.toList());

            final List<Histogram<Integer>> altDepthOneHistogramsRevComp = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(reverseComplement))
                    .collect(Collectors.toList());

            final List<AltSiteRecord> altDesignMatrix = altDesignMatrixByContext
                    .getOrDefault(refContext, new ArrayList<>()); // Cannot use Collections.emptyList() here because we might add to it
            final List<AltSiteRecord> altDesignMatrixRevComp = altDesignMatrixByContext
                    .getOrDefault(reverseComplement, Collections.emptyList());
            // Warning: the below method will mutate the content of {@link altDesignMatrixRevComp} and append to {@code altDesignMatrix}
            mergeDesignMatrices(altDesignMatrix, altDesignMatrixRevComp);

            final Histogram<Integer> combinedRefHistograms = combineRefHistogramWithRC(refContext, refHistogram, refHistogramRevComp);
            final List<Histogram<Integer>> combinedAltHistograms = combineAltDepthOneHistogramWithRC(altDepthOneHistogramsForContext, altDepthOneHistogramsRevComp);

            if (combinedRefHistograms.getSumOfValues() == 0 || altDesignMatrix.isEmpty()) {
                logger.info(String.format("Skipping the reference context %s as we didn't find either the ref or alt table for the context", refContext));
                continue;
            }

            final LearnReadOrientationModelEngine engine = new LearnReadOrientationModelEngine(
                    combinedRefHistograms,
                    combinedAltHistograms,
                    altDesignMatrix,
                    converagenceThreshold,
                    maxEMIterations,
                    logger);
            final ArtifactPrior artifactPrior = engine.learnPriorForArtifactStates();
            artifactPriorAcrossContexts.add(artifactPrior);
        }

        ArtifactPrior.writeArtifactPriors(artifactPriorAcrossContexts, output);
        return "SUCCESS";
    }

    @VisibleForTesting
    public static Histogram<Integer> combineRefHistogramWithRC(final String refContext,
                                                               final Histogram<Integer> refHistogram,
                                                               final Histogram<Integer> refHistogramRevComp){
        Utils.validateArg(refHistogram.getValueLabel()
                .equals(SequenceUtil.reverseComplement(refHistogramRevComp.getValueLabel())),
                "ref context = " + refHistogram.getValueLabel() + ", rev comp = " + refHistogramRevComp.getValueLabel());
        Utils.validateArg(refHistogram.getValueLabel().equals(refContext), "this better match");

        final Histogram<Integer> combinedRefHistogram = F1R2FilterUtils.createRefHistogram(refContext);

        for (final Integer depth : refHistogram.keySet()){
            final double newCount = refHistogram.get(depth).getValue() + refHistogramRevComp.get(depth).getValue();
            combinedRefHistogram.increment(depth, newCount);
        }

        return combinedRefHistogram;
    }

    @VisibleForTesting
    public static List<Histogram<Integer>> combineAltDepthOneHistogramWithRC(final List<Histogram<Integer>> altHistograms,
                                                                             final List<Histogram<Integer>> altHistogramsRevComp){
        if (altHistograms.isEmpty() && altHistogramsRevComp.isEmpty()){
            return Collections.emptyList();
        }

        final String refContext = ! altHistograms.isEmpty() ?
                F1R2FilterUtils.labelToTriplet(altHistograms.get(0).getValueLabel()).getLeft() :
                SequenceUtil.reverseComplement(F1R2FilterUtils.labelToTriplet(altHistogramsRevComp.get(0).getValueLabel()).getLeft());

        // Contract: altHistogram must be of the canonical representation of the kmer
        Utils.validateArg(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "refContext must be the canonical representation but got " + refContext);

        final List<Histogram<Integer>> combinedHistograms = new ArrayList<>(F1R2FilterConstants.numAltHistogramsPerContext);

        for (Nucleotide altAllele : F1R2FilterConstants.REGULAR_BASES){
            // Skip when the alt base is the ref base, which doesn't make sense because this is a histogram of alt sites
            if (altAllele == F1R2FilterUtils.getMiddleBase(refContext)){
                continue;
            }

            final String reverseComplement = SequenceUtil.reverseComplement(refContext);
            final Nucleotide altAlleleRevComp = Nucleotide.valueOf(SequenceUtil.reverseComplement(altAllele.toString()));

            for (ReadOrientation orientation : ReadOrientation.values()) {
                final ReadOrientation oppositeType = orientation == ReadOrientation.F1R2 ? ReadOrientation.F2R1 : ReadOrientation.F1R2;
                final Histogram<Integer> altHistogram = altHistograms.stream()
                        .filter(h -> h.getValueLabel().equals(F1R2FilterUtils.tripletToLabel(refContext, altAllele, orientation)))
                        .findFirst().orElseGet(() -> F1R2FilterUtils.createAltHistogram(refContext, altAllele, orientation));

                final Histogram<Integer> altHistogramRevComp = altHistogramsRevComp.stream()
                        .filter(h -> h.getValueLabel().equals(F1R2FilterUtils.tripletToLabel(reverseComplement, altAlleleRevComp, oppositeType)))
                        .findFirst().orElseGet(() -> F1R2FilterUtils.createAltHistogram(reverseComplement, altAlleleRevComp, oppositeType));

                final Histogram<Integer> combinedHistogram = F1R2FilterUtils.createAltHistogram(refContext, altAllele, orientation);

                // Add the histograms manually - I don't like the addHistogram() in htsjdk method because it does so with side-effect
                for (final Integer depth : altHistogram.keySet()){
                    final double newCount = altHistogram.get(depth).getValue() + altHistogramRevComp.get(depth).getValue();
                    combinedHistogram.increment(depth, newCount);
                }

                combinedHistograms.add(combinedHistogram);
            }
        }

        return combinedHistograms;
    }


    /**
     * Contract: this method must be called after grouping the design matrices by context.
     * That is, {@param altDesignMatrix} and {@param altDesigmRevComp} must each be a list of {@link AltSiteRecord}
     * of a single reference context, and their reference contexts must match up to reverse complement
     */
    @VisibleForTesting
    public static void mergeDesignMatrices(final List<AltSiteRecord> altDesignMatrix, List<AltSiteRecord> altDesignMatrixRevComp){
        if (altDesignMatrix.isEmpty() && altDesignMatrixRevComp.isEmpty()){
            return;
        }

        // Order matters here. Assumes that all elements in the list have the same reference context
        Utils.validateArg(altDesignMatrix.isEmpty() || F1R2FilterConstants.CANONICAL_KMERS.contains(altDesignMatrix.get(0).getReferenceContext()),
                "altDesignMatrix must have the canonical representation");

        final Optional<String> refContext = altDesignMatrix.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrix.get(0).getReferenceContext());
        final Optional<String> revCompContext = altDesignMatrixRevComp.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrixRevComp.get(0).getReferenceContext());
        if (refContext.isPresent() && revCompContext.isPresent()){
            Utils.validateArg(refContext.get().equals(SequenceUtil.reverseComplement(revCompContext.get())),
                    "ref context and its rev comp don't match");
        }

        altDesignMatrix.addAll(altDesignMatrixRevComp.stream()
                .map(AltSiteRecord::getReverseComplementOfRecord)
                .collect(Collectors.toList()));
    }

    private MetricsFile<?, Integer> readMetricsFile(File file){
        final MetricsFile<?, Integer> metricsFile = new MetricsFile<>();
        final Reader reader = IOUtil.openFileForBufferedReading(file);
        metricsFile.read(reader);
        CloserUtil.close(reader);
        return metricsFile;
    }

}
