package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * This enum encapsulates the domain of the discrete latent random variable z
 */
public enum ArtifactState {
    // Artifact States
    F1R2_A, F1R2_C, F1R2_G, F1R2_T,
    F2R1_A, F2R1_C, F2R1_G, F2R1_T,

    // Non Artifact States
    HOM_REF, GERMLINE_HET, SOMATIC_HET, HOM_VAR;

    public static List<ArtifactState> getStates(){
        return Arrays.stream(ArtifactState.values()).collect(Collectors.toList());
    }

    static List<ArtifactState> getF1R2ArtifactStates(){
        return Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T);
    }

    static List<ArtifactState> getF2R1ArtifactStates(){
        return Arrays.asList(F2R1_A, F2R1_C, F2R1_G, F2R1_T);
    }

    public static List<ArtifactState> getNonArtifactStates(){
        return Arrays.asList(HOM_REF, GERMLINE_HET, SOMATIC_HET, HOM_VAR);
    }

    // For a given reference base, return the ref to ref artifact states (e.g. AGT -> G), which we want to ignore
    public static List<ArtifactState> getRefToRefArtifacts(final Nucleotide refAllele){
        switch (refAllele){
            case A : return Arrays.asList( F1R2_A, F2R1_A );
            case C : return Arrays.asList( F1R2_C, F2R1_C );
            case G : return Arrays.asList( F1R2_G, F2R1_G );
            case T : return Arrays.asList( F1R2_T, F2R1_T );
            default: throw new UserException(String.format("Invalid nucleotide given: %s", refAllele));
        }
    }

    // Given a state z, return the alt allele of the artifact that the state encodes
    public static Nucleotide getAltAlleleOfArtifact(final ArtifactState z){
        Utils.validateArg(Arrays.asList(ArtifactState.F1R2_A, ArtifactState.F1R2_C, ArtifactState.F1R2_G, ArtifactState.F1R2_T,
                ArtifactState.F2R1_A, ArtifactState.F2R1_C, ArtifactState.F2R1_G, ArtifactState.F2R1_T).contains(z),
                String.format("ArtifactState must be F1R2_a or F2R1_a but got %s", z));
        switch (z){
            case F1R2_A : return Nucleotide.A;
            case F1R2_C : return Nucleotide.C;
            case F1R2_G : return Nucleotide.G;
            case F1R2_T : return Nucleotide.T;
            case F2R1_A : return Nucleotide.A;
            case F2R1_C : return Nucleotide.C;
            case F2R1_G : return Nucleotide.G;
            case F2R1_T : return Nucleotide.T;
            default: throw new UserException(String.format("Invalid state: %s", z));
        }
    }

    static List<ArtifactState> artifactStates = Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T, F2R1_A, F2R1_C, F2R1_G, F2R1_T);


    public static ArtifactState getF1R2StateForAlt(final Nucleotide altAllele) {
        switch (altAllele) {
            case A:
                return ArtifactState.F1R2_A;
            case C:
                return ArtifactState.F1R2_C;
            case G:
                return ArtifactState.F1R2_G;
            case T:
                return ArtifactState.F1R2_T;
            default:
                throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
        }
    }

    public static ArtifactState getF2R1StateForAlt(final Nucleotide altAllele) {
        switch (altAllele) {
            case A:
                return ArtifactState.F2R1_A;
            case C:
                return ArtifactState.F2R1_C;
            case G:
                return ArtifactState.F2R1_G;
            case T:
                return ArtifactState.F2R1_T;
            default:
                throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
        }
    }

    public static ArtifactState getRevCompState(final ArtifactState state){
        switch (state) {
            case F1R2_A:
                return F2R1_T;
            case F1R2_C:
                return F2R1_G;
            case F1R2_G:
                return F2R1_C;
            case F1R2_T:
                return F2R1_A;
            case F2R1_A:
                return F1R2_T;
            case F2R1_C:
                return F1R2_G;
            case F2R1_G:
                return F1R2_C;
            case F2R1_T:
                return F1R2_A;
            default:
                return state;
        }
    }
}
