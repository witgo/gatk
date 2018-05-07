package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.primitives.Ints;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.tools.walkers.readorientation.AltSiteRecord.AltSiteRecordTableWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientation.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientation.F2R1;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2FilterConstants.*;

/**
 * {@link LearnReadOrientationModel} takes in the tsv output of this tool
 */

@CommandLineProgramProperties(
        summary = "Collect data from a tumor bam for Mutect2 Read Orientation Filter",
        oneLineSummary = "Data collection for Mutect2 Read Orientation Filter",
        programGroup = CoverageAnalysisProgramGroup.class
)

public class CollectF1R2Counts extends LocusWalker {
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt-table";
    public static final String ALT_DEPTH1_HISTOGRAM_LONG_NAME = "alt-hist";
    public static final String REF_SITE_METRICS_LONG_NAME = "ref-table";
    public static final String MIN_MEDIAN_MQ_LONG_NAME = "median-mq";
    public static final String MIN_BASE_QUALITY_LONG_NAME = "min-bq";

    @Argument(fullName = MIN_MEDIAN_MQ_LONG_NAME, doc = "skip sites with median mapping quality below this value", optional = true)
    private static int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = MIN_BASE_QUALITY_LONG_NAME, doc = "exclude bases below this quality from pileup", optional = true)
    private static int MINIMUM_BASE_QUALITY = 20;

    @Argument(fullName = ALT_DATA_TABLE_LONG_NAME, doc = "a tab-separated output table of pileup data over alt sites")
    private static File altDataTable = null;

    @Argument(fullName = REF_SITE_METRICS_LONG_NAME, doc = "a metrics file with overall summary metrics and reference context-specific depth histograms")
    private static File refMetricsOutput = null;

    @Argument(fullName = ALT_DEPTH1_HISTOGRAM_LONG_NAME, doc = "a histogram of alt sites with alt depth = 1")
    private static File altMetricsOutput = null;

    // For computational efficiency, for each reference context, we build a depth histogram over ref sites
    private static Map<String, Histogram<Integer>> refSiteHistograms = new HashMap<>(ALL_KMERS.size());

    // Context -> (Alt, F1R2) -> Histogram
    private static Map<String, Map<Pair<Nucleotide, ReadOrientation>, Histogram<Integer>>> depthOneAltHistograms = new HashMap<>(ALL_KMERS.size());

    private AltSiteRecordTableWriter altTableWriter;

    private final MetricsFile<?, Integer> refMetricsFile = getMetricsFile();

    private final MetricsFile<?, Integer> altMetricsFile = getMetricsFile();

    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
        // Initialize for each reference the histogram of the counts of reference sites by depth
        ALL_KMERS.forEach(context -> {
            Histogram<Integer> emptyRefHistogram = F1R2FilterUtils.createRefHistogram(context);
            refSiteHistograms.put(context, emptyRefHistogram);
        });

        // Initialize for each reference context the (Alt Allele, Artifact Type) -> Histogram map
        ALL_KMERS.forEach(context -> {
            depthOneAltHistograms.put(context, new HashMap<>((REGULAR_BASES.size() - 1)*ReadOrientation.values().length));

            for (Nucleotide altAllele : REGULAR_BASES){
                // Skip e.g. AGT -> AGT because G is not an alt allele
                if (altAllele == Nucleotide.valueOf(context.substring(MIDDLE_INDEX, MIDDLE_INDEX+1))){
                    continue;
                }

                for (ReadOrientation artifactType : ReadOrientation.values()){
                    depthOneAltHistograms.get(context).put(new Pair<>(altAllele, artifactType),
                            F1R2FilterUtils.createAltHistogram(context, altAllele, artifactType));
                }
            }
        });

        // Intentionally not use try-with-resources so that the writer stays open outside of the try block
        try {
            altTableWriter = new AltSiteRecordTableWriter(altDataTable);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception creating a writer for %s", altDataTable), e);
        }

    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final int position = referenceContext.getInterval().getStart();
        final String refContext = referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING);
        // final String refContext = Arrays.toString(referenceContext.getKmerAround(position, F1R2FilterConstants.REF_CONTEXT_PADDING));
        final Nucleotide refBase = F1R2FilterUtils.getMiddleBase(refContext);

        // START deprecated code
        referenceContext.setWindow(F1R2FilterConstants.REF_CONTEXT_PADDING, F1R2FilterConstants.REF_CONTEXT_PADDING);
        final String oldRefContext = new String(referenceContext.getBases());
        Utils.validate(refContext.equals(oldRefContext), String.format("old ref context %s and new ref context %s disagree", oldRefContext, refContext));
        // END deprecated code

        if (refContext.contains("N") || refContext.length() != F1R2FilterConstants.REFERENCE_CONTEXT_SIZE) {
            return;
        }

        if (refContext == null) {
            logger.warn(String.format("Skipped a site with null reference at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        final ReadPileup pileup = alignmentContext.getBasePileup().makeFilteredPileup(pe -> pe.getQual() > MINIMUM_BASE_QUALITY);
        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        if (!isPileupGood(pileup, referenceContext)) {
            return;
        }

        // Make a copy of base counts and update the counts of ref to -1. Now the maxElementIndex of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = -1;
        final int altBaseIndex = MathUtils.maxElementIndex(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        // If the site is ref, we simply update the coverage histogram
        if (referenceSite) {
            refSiteHistograms.get(refContext).increment(depth <= F1R2FilterConstants.maxDepth ?
                    depth : F1R2FilterConstants.maxDepth);
            return;
        }

        // If we got here, we have an alt site with a single alt base
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        final int refCount = baseCounts[refBase.ordinal()];
        final int altCount = baseCounts[altBaseIndex];
        Utils.validate(altCount > 0, "We must have a nonzero alt read but got " + altCount);

        final int refF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.valueOf(pe.getBase()) == refBase && ReadUtils.isF1R2(pe.getRead()));
        final int altF1R2 = pileup.getNumberOfElements(pe -> Nucleotide.valueOf(pe.getBase()) == altBase && ReadUtils.isF1R2(pe.getRead()));

        // If there's only one alt read at depth > minimumDepthForOptimization, we store the data differently to save space
        if (altCount == 1) {
            final ReadOrientation type = altF1R2 == 1 ? F1R2 : F2R1;
            depthOneAltHistograms.get(refContext).get(new Pair<>(altBase, type))
                    .increment(depth <= F1R2FilterConstants.maxDepth ? depth : F1R2FilterConstants.maxDepth);
            return;
        }

        try {
            altTableWriter.writeRecord(new AltSiteRecord(refContext, refCount, altCount, refF1R2, altF1R2, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        refSiteHistograms.values().forEach(h -> refMetricsFile.addHistogram(h));
        refMetricsFile.write(refMetricsOutput);

        depthOneAltHistograms.values().forEach(hs -> hs.values().forEach(h -> altMetricsFile.addHistogram(h)));
        altMetricsFile.write(altMetricsOutput);

        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (altTableWriter != null) {
            try {
                altTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the alt table writer", e);
            }
        }
    }

    /**
     * Use a series of heuristics to detect a bad pileup.
     */
    private boolean isPileupGood(final ReadPileup pileup, final ReferenceContext referenceContext){
        // This case should not happen, as AlignmentContext should come filtered, but it does happen once in a while

        final int[] baseCounts = pileup.getBaseCounts();
        final int depth = (int) MathUtils.sum(baseCounts);

        // If observe any bases that are neither ref nor alt, the site is no good so filter
        final long numAlleles = Arrays.stream(baseCounts).filter(i -> i > 0).count();

        List<Integer> mappingQualities = Ints.asList(pileup.getMappingQuals());

        boolean isIndel = pileup.getNumberOfElements(pe -> pe.isDeletion() || pe.isAfterInsertion() || pe.isBeforeDeletionStart()) > 0;

        // If depth (the sum of base counts) is 0 but the pileup is non-empty, that means all the reads
        // have deleted bases at this particular locus
        isIndel = isIndel || depth == 0 && pileup.size() > 0;

        return depth > 0 && ! isIndel && MathUtils.median(mappingQualities) >= MINIMUM_MEDIAN_MQ && numAlleles <= 2;

    }
}
