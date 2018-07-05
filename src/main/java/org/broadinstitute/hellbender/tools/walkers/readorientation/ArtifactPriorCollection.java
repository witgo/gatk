package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;

/**
 * Container class for {@link ArtifactPrior} objects. The container will always have {@link F1R2FilterConstants.NUM_KMERS} objects
 */
public class ArtifactPriorCollection {
    final private Map<String, ArtifactPrior> map = new HashMap<>(F1R2FilterConstants.NUM_KMERS);

    public ArtifactPriorCollection(){
        for (final String kmer : F1R2FilterConstants.ALL_KMERS){
            map.put(kmer, new ArtifactPrior(kmer, new double[F1R2FilterConstants.NUM_STATES], 0, 0));
        }
    }

    public Optional<ArtifactPrior> get(final String refContext){
        if (map.containsKey(refContext)){
            return Optional.of(map.get(refContext));
        } else {
            return Optional.empty();
        }
    }

    /**
     *
     * @param artifactPrior
     */
    public void set(final ArtifactPrior artifactPrior){
        final String refContext = artifactPrior.getReferenceContext();
        Utils.validate(map.get(refContext).getNumExamples() == 0,
                "updating an existing ArtifactPrior is not allowed. Ref context = " + refContext);
        Utils.validate(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "set must be called on an artifactPrior object with a canonical representation");

        map.put(refContext, artifactPrior);

        final String revCompRefContext = SequenceUtil.reverseComplement(refContext);
        final ArtifactPrior revCompArtifactPrior = artifactPrior.getReverseComplement();
        map.put(revCompRefContext, revCompArtifactPrior);
    }

    public void writeArtifactPriors(final File output){
        ArtifactPrior.writeArtifactPriors(new ArrayList<>(map.values()), output);
    }

    public static ArtifactPriorCollection readArtifactPriors(final File input){
        final List<ArtifactPrior> priors = ArtifactPrior.readArtifactPriors(input);
        final ArtifactPriorCollection artifactPriorCollection = new ArtifactPriorCollection();

        /**
         * We iterate through the canonical kmers instead of all reference contexts because otherwise we would
         * visit each reference context twice and get an error the second time, as updating artifact prior that already
         * exists in the container class is prohibited
         */
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final Optional<ArtifactPrior> ap = priors.stream().filter(a -> a.getReferenceContext().equals(refContext)).findAny();
            Utils.validate(ap.isPresent(), "ArtifactPrior object isn't present for reference context " + refContext + "in file " + input);
            artifactPriorCollection.set(ap.get());
        }
        return artifactPriorCollection;
    }

    public int getNumUniqueContexts(){
        return (int) map.values().stream().filter(a -> a.getNumExamples() > 0).count()/2;
    }
}
