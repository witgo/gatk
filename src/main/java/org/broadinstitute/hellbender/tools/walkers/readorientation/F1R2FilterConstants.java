package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by tsato on 3/14/18.
 */
public class F1R2FilterConstants {
    // Padding on each side
    public static final int REF_CONTEXT_PADDING = 1;
    public static final int NUM_STATES = ArtifactState.values().length;

    static final int MIDDLE_INDEX = REF_CONTEXT_PADDING;

    static final int REFERENCE_CONTEXT_SIZE = 2 * REF_CONTEXT_PADDING + 1; // aka 3

    public static final List<Nucleotide> REGULAR_BASES = Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T);

    // the list of all possible kmers, where k = REFERENCE_CONTEXT_SIZE
    static final List<String> ALL_KMERS = SequenceUtil.generateAllKmers(REFERENCE_CONTEXT_SIZE).stream()
            .map(String::new).collect(Collectors.toList());

    // For each pair of a K-mer and its reverse complement, pick lexicographically smaller K-mer
    // to be the canonical representation of the pair
    static final List<String> CANONICAL_KMERS = ALL_KMERS.stream()
            .map(context -> new TreeSet<>(Arrays.asList(context, SequenceUtil.reverseComplement(context))))
            .distinct()
            .map(s -> s.first().compareTo(s.last()) < 0 ? s.first() : s.last())
            .collect(Collectors.toList());

    // Delimiter of the alt histogram label
    static final String FIELD_SEPARATOR = "_";

    static final String binName = "depth";

    // We combine all sites of depths above this value in the last bin of the histogram
    static final int maxDepth = 200;

    static final Integer[] emptyBins = IntStream.rangeClosed(1, F1R2FilterConstants.maxDepth).boxed().toArray( Integer[]::new );

    static final int numAltHistogramsPerContext = (REGULAR_BASES.size() - 1) * (ReadOrientation.values().length);

}
