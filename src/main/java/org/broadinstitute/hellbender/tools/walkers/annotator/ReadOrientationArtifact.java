package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.readorientation.*;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.F1R2FilterConstants.REF_CONTEXT_PADDING;

/**
 * Created by tsato on 10/20/17.
 */
public class ReadOrientationArtifact extends GenotypeAnnotation implements StandardMutectAnnotation {
    private File artifactPriorTable;
    private List<ArtifactPrior> artifactPriors = Collections.emptyList();
    private int minimumBaseQuality = 20;

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.ROF_POSTERIOR_KEY, GATKVCFConstants.ROF_PRIOR_KEY, GATKVCFConstants.ROF_TYPE_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(GATKVCFConstants.ROF_POSTERIOR_KEY, 1,
                        VCFHeaderLineType.Float, "posterior probability of read orientaion artifact"),
                new VCFFormatHeaderLine(GATKVCFConstants.ROF_PRIOR_KEY, 1,
                        VCFHeaderLineType.Float, "prior probability of read oientation artifact under the present reference context"),
                new VCFFormatHeaderLine(GATKVCFConstants.ROF_TYPE_KEY, 1,
                        VCFHeaderLineType.String, "F1R2 or F2R1"));
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        Utils.nonNull(likelihoods);
        final boolean normalSample = g.isHomRef();

        /**
         * Skip the following cases:
         *  1. Normal sample
         *  2. Non-SNP variants
         *  3. The table of prior is not provided
         */
        if (normalSample || !vc.isSNP() || artifactPriors.isEmpty()){
            return;
        }

        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY, () -> null, -1);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
        final Allele altAllele = vc.getAlternateAllele(indexOfMaxTumorLod);
        final Nucleotide altBase = Nucleotide.valueOf(altAllele.toString());


        final String refContext = ref.getKmerAround(vc.getStart(), REF_CONTEXT_PADDING);
        Utils.validate(refContext.length() == 2* REF_CONTEXT_PADDING + 1,
                String.format("kmer must have length %d but got %d", 2* REF_CONTEXT_PADDING + 1, refContext.length()));

        if (refContext.contains("N")) {
            return;
        }

        final Nucleotide refAllele = F1R2FilterUtils.getMiddleBase(refContext);
        Utils.validate(refAllele == Nucleotide.valueOf(vc.getReference().toString().replace("*", "")),
                String.format("ref allele in the kmer, %s, does not match the ref allele in the variant context, %s",
                        refAllele, vc.getReference().toString().replace("*", "")));

        int refCount = 0;
        int altCount = 0;
        int altF1R2 = 0;

        final Collection<ReadLikelihoods<Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(g.getSampleName());
        for (ReadLikelihoods<Allele>.BestAllele bestAllele : bestAlleles){
            final Allele allele = bestAllele.allele;

            if (!bestAllele.isInformative() || allele.length() != 1){
                continue;
            }

            final GATKRead read = bestAllele.read;

            // We will throw out any read that does not overlap the locus (this happens, presumably, when there are other reads
            // that are in phase with this particular read and the engine imputes the missing base according to their haplotype)
            if (read.getEnd() < vc.getStart() || read.getStart() > vc.getStart()){
                continue;
            }

            // Use AlignmentStateMachine to find the base quality of the base at the position
            // TODO: perhaps using the read likelihood is cheaper and serves the same purpose
            final AlignmentStateMachine asm = new AlignmentStateMachine(read);
            while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < vc.getStart()) { }

            final int readOffset = asm.getReadOffset();

            // Throw away bases that is below the minimum quality
            if (asm.getGenomePosition() == vc.getStart() && read.getBaseQuality(readOffset) < minimumBaseQuality){
                continue;
            }

            if (allele.isReference()){
                refCount++;
            } else if (allele.equals(altAllele)){
                altCount++;
                if (ReadUtils.isF1R2(bestAllele.read)){
                    altF1R2++;
                }
            }
        }

        final Optional<ArtifactPrior> hyps = ArtifactPrior.searchByContext(artifactPriors, refContext);

        if (! hyps.isPresent()){
            return;
        }
        
        final double[] artifactPrior = hyps.get().getPi(refContext);
        final int depth = refCount + altCount;

        final double[] log10UnnormalizedPosteriorProbabilities = LearnReadOrientationModelEngine.computeLog10Responsibilities(
                refAllele, altBase, altCount, altF1R2, depth, artifactPrior);

        // We want the posterior of artifacts given that the site is not hom ref
        log10UnnormalizedPosteriorProbabilities[ArtifactState.HOM_REF.ordinal()] = Double.NEGATIVE_INFINITY;

        final double[] posterior = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriorProbabilities);

        final double posteriorOfF1R2 = posterior[ArtifactState.getF1R2StateForAlt(altBase).ordinal()];
        final double posteriorOfF2R1 = posterior[ArtifactState.getF2R1StateForAlt(altBase).ordinal()];

        final double posteriorOfArtifact = Math.max(posteriorOfF1R2, posteriorOfF2R1);
        final ReadOrientation artifactType = posteriorOfF1R2 > posteriorOfF2R1 ? ReadOrientation.F1R2 : ReadOrientation.F2R1;

        gb.attribute(GATKVCFConstants.ROF_POSTERIOR_KEY, posteriorOfArtifact);
        final int indexOfArtifact = artifactType == ReadOrientation.F1R2 ?
                ArtifactState.getF1R2StateForAlt(altBase).ordinal() : ArtifactState.getF2R1StateForAlt(altBase).ordinal();
        gb.attribute(GATKVCFConstants.ROF_PRIOR_KEY, artifactPrior[indexOfArtifact]);
        gb.attribute(GATKVCFConstants.ROF_TYPE_KEY, artifactType.toString());
    }

    /** Part of the hack scheme that is required until we can pass an argument directly ot an annotation class **/
    public void setPriorArtifactTable(final File artifactPriorTable){
        this.artifactPriorTable = artifactPriorTable;
        artifactPriors = ArtifactPrior.readArtifactPriors(artifactPriorTable);
    }

}
