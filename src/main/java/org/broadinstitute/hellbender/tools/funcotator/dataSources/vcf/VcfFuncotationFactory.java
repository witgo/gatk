package org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A class to create annotations from VCF feature sources.
 * Created by jonn on 3/23/18.
 */
public class VcfFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(VcfFuncotationFactory.class);

    //==================================================================================================================
    // Private Members:

    /**
     * Name of this {@link VcfFuncotationFactory}.
     */
    private final String name;

    /**
     * The file on which this {@link VcfFuncotationFactory} is based.
     */
    private final Path sourceFilePath;

    /**
     * The field names that this {@link VcfFuncotationFactory} supports
     * and default values for each.
     */
    private final LinkedHashMap<String, Object> supportedFieldNamesAndDefaults;

    /**
     * A list of values to use when there are no annotations for an allele.
     */
    private final LinkedHashSet<String> supportedFieldNames;

    /**
     * Should contain metadata only for the fields in supportedFieldNames.
     */
    private final FuncotationMetadata supportedFieldMetadata;

    private final LRUCache<Triple<VariantContext, ReferenceContext, List<Feature>>, List<Funcotation>> cache = new LRUCache<>();
    private int cacheHits = 0;
    private int cacheMisses = 0;

    //==================================================================================================================
    // Constructors:

    public VcfFuncotationFactory(final String name, final String version, final Path sourceFilePath) {
        this(name, version, sourceFilePath, new LinkedHashMap<>());
    }

    public VcfFuncotationFactory(final String name, final String version, final Path sourceFilePath, final LinkedHashMap<String, String> annotationOverridesMap) {
        this.name = name;
        this.version = version;
        this.sourceFilePath = sourceFilePath;

        // Handle the supported field names here:
        supportedFieldNamesAndDefaults = new LinkedHashMap<>();
        supportedFieldNames = new LinkedHashSet<>();
        populateSupportedFieldNamesFromVcfFile();

        // This step has to occur after supported field names and name have been populated.
        supportedFieldMetadata = createFuncotationMetadata(sourceFilePath);

        if (supportedFieldNames.size() == 0) {
            logger.warn("WARNING: VcfFuncotationFactory has nothing to annotate from VCF File: " + sourceFilePath.toUri().toString());
        }
        else {
            // Now check if we have any overrides to take care of:
            this.annotationOverrideMap = new LinkedHashMap<>();
            for ( final Map.Entry<String, String> entry : annotationOverridesMap.entrySet() ) {
                if ( supportedFieldNamesAndDefaults.containsKey(entry.getKey()) ) {
                    annotationOverrideMap.put(entry.getKey(), entry.getValue());
                }
            }
        }
    }

    private FuncotationMetadata createFuncotationMetadata(final Path sourceFilePath) {
        // Read the VCF to just get the header
        try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(sourceFilePath.toString()) ) {
            final Object header = vcfReader.getHeader();
            if ( ! (header instanceof VCFHeader) ) {
                throw new IllegalArgumentException(sourceFilePath + " does not have a valid VCF header");
            }
            final VCFHeader sourceVcfHeader = (VCFHeader) header;
            final List<VCFInfoHeaderLine> metadataVcfInfoHeaderLines = createFuncotationVcfInfoHeaderLines(sourceVcfHeader);
            return VcfFuncotationMetadata.create(metadataVcfInfoHeaderLines);
        }
    }

    @VisibleForTesting
    List<VCFInfoHeaderLine> createFuncotationVcfInfoHeaderLines(final VCFHeader vcfHeader) {
        final List<VCFInfoHeaderLine> supportedVcfInfoHeaderLines = vcfHeader.getInfoHeaderLines().stream()
                .filter(vcfInfoHeaderLine -> supportedFieldNames.contains(createFinalFieldName(name, vcfInfoHeaderLine.getID())))
                .collect(Collectors.toList());

        // Make sure to rename the input VCF field names to the output funcotation field names for this funcotation factory.
        return supportedVcfInfoHeaderLines.stream()
                .map(vcfInfoHeaderLine -> copyWithRename(vcfInfoHeaderLine, name))
                .collect(Collectors.toList());
    }

    private static VCFInfoHeaderLine copyWithRename(final VCFInfoHeaderLine vcfInfoHeaderLine, final String name) {
        if (vcfInfoHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER) {
            return new VCFInfoHeaderLine(createFinalFieldName(name, vcfInfoHeaderLine.getID()),
                    vcfInfoHeaderLine.getCount(), vcfInfoHeaderLine.getType(), vcfInfoHeaderLine.getDescription());
        } else {
            return new VCFInfoHeaderLine(createFinalFieldName(name, vcfInfoHeaderLine.getID()),
                    vcfInfoHeaderLine.getCountType(), vcfInfoHeaderLine.getType(), vcfInfoHeaderLine.getDescription());
        }
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    protected Class<? extends Feature> getAnnotationFeatureClass() {
        return VariantContext.class;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.VCF;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return supportedFieldNames;
    }

    @Override
    protected List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext ) {
        if ( supportedFieldNames.size() != 0 ) {
            return createDefaultFuncotationsOnVariantHelper(variant, referenceContext, Collections.emptySet());
        }
        else {
            return Collections.emptyList();
        }
    }

    @Override
    /**
     * {@inheritDoc}
     *
     * This is really the only entry point for the Vcf FuncotationFactory.
     *
     * {@link VcfFuncotationFactory} can be used with or without Gencode annotations.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {

        final List<Funcotation> outputFuncotations = new ArrayList<>();

        final Triple<VariantContext, ReferenceContext, List<Feature>> cacheKey = createCacheKey(variant, referenceContext, featureList);
        final List<Funcotation> cacheResult = cache.get(cacheKey);
        if (cacheResult != null) {
            cacheHits++;
            return cacheResult;
        }

        // Only create annotations if we have data to annotate:
        if ( supportedFieldNames.size() != 0 ) {

            // Alternate alleles in the query variant.
            final List<Allele> queryAlternateAlleles = variant.getAlternateAlleles();

            // Create a set to put our annotated Alternate alleles in.
            // We'll use this to determine if the alt allele has been annotated.
            final Set<Allele> annotatedAltAlleles = new HashSet<>(queryAlternateAlleles.size());

            // Get rid of any null features.
            // By this point we know the feature type is correct, so we cast it:
            final List<VariantContext> funcotationFactoryVariants = featureList.stream().filter(f -> f != null)
                    .map(f -> (VariantContext) f).collect(Collectors.toList());

            for ( final VariantContext funcotationFactoryVariant : funcotationFactoryVariants ) {
                final List<Pair<Integer, Integer>> matchIndices = matchAlleles(variant, funcotationFactoryVariant);

                for (final Pair<Integer, Integer> matchIndex : matchIndices) {

                    final Allele queryAllele = variant.getAlternateAllele(matchIndex.getLeft());
                    final Allele funcotationFactoryAllele = funcotationFactoryVariant.getAlternateAllele(matchIndex.getRight());

                    final LinkedHashMap<String, Object> annotations = new LinkedHashMap<>(supportedFieldNamesAndDefaults);

                    for ( final Map.Entry<String, Object> entry : funcotationFactoryVariant.getAttributes().entrySet() ) {
                        populateAnnotationMap(funcotationFactoryVariant, variant, matchIndex.getRight(), annotations, entry);
                    }

                    outputFuncotations.add(TableFuncotation.create(annotations, queryAllele, name, supportedFieldMetadata));
                    annotatedAltAlleles.add(funcotationFactoryAllele);
                }
            }

            // If we didn't add funcotations for an allele, we should add in blank funcotations to that allele for each field that can be produced
            // by this VcfFuncotationFactory:
            if ( annotatedAltAlleles.size() != queryAlternateAlleles.size() ) {
                outputFuncotations.addAll( createDefaultFuncotationsOnVariantHelper(variant, referenceContext, annotatedAltAlleles) );
            }
        }
        cacheMisses++;
        cache.put(cacheKey, outputFuncotations);
        return outputFuncotations;
    }

    private void populateAnnotationMap(final VariantContext funcotationFactoryVariant, final VariantContext queryVariant, final int funcotationFactoryAltAlleleIndex, final LinkedHashMap<String, Object> annotations, final Map.Entry<String, Object> attributeEntry) {
        final String valueString;
        final String attributeName = attributeEntry.getKey();

        // Handle collections a little differently:
        if (attributeEntry.getValue() instanceof Collection<?>) {
            @SuppressWarnings("unchecked") final Collection<Object> objectList = ((Collection<Object>) attributeEntry.getValue());
            final VCFHeaderLineCount countType = supportedFieldMetadata.retrieveHeaderInfo(createFinalFieldName(this.name, attributeName)).getCountType();

            if (isBiallelic(funcotationFactoryVariant) && (isBiallelic(queryVariant))) {
                valueString = objectList.stream().map(Object::toString).collect(Collectors.joining(","));
            } else {
                valueString = determineValueStringFromMultiallelicAttributeList(funcotationFactoryVariant, funcotationFactoryAltAlleleIndex, objectList, countType);
            }
        } else {
            valueString = attributeEntry.getValue().toString();
        }

        annotations.put(createFinalFieldName(name, attributeName), valueString);
    }

    private boolean isBiallelic(final VariantContext funcotationFactoryVariant) {
        return funcotationFactoryVariant.getAlternateAlleles().size() == 1;
    }

    private String determineValueStringFromMultiallelicAttributeList(final VariantContext funcotationFactoryVariantContext, final int funcotationFactoryAltAlleleIndex, final Collection<Object> objectList, final VCFHeaderLineCount countType) {
        String result;
        //TODO: Refactor to a switch
        if (countType.equals(VCFHeaderLineCount.A) || countType.equals(VCFHeaderLineCount.R)) {

            // TODO: What about "R"?  Do we want to drop the reference number?
            int idx = funcotationFactoryAltAlleleIndex;
            if (countType.equals(VCFHeaderLineCount.R)) {
                idx++;
            }
            result = objectList.toArray()[idx].toString();
        } else {
            result = objectList.stream().map(Object::toString).collect(Collectors.joining(","));
        }
        return result;
    }

    // Returns indexes into the alternate alleles.  Note that this method assumes that (when biallelic) the variant
    //  contexts are already trimmed.
    private List<Pair<Integer, Integer>> matchAlleles(final VariantContext variant1, VariantContext variant2) {

        // Grab the trivial case:
        if (isBiallelic(variant1) && isBiallelic(variant2)) {
            if (variant1.getAlternateAllele(0).equals(variant2.getAlternateAllele(0)) &&
                    (variant1.getReference().equals(variant2.getReference()))) {
                return Collections.singletonList(Pair.of(0,0));
            } else {
                return Collections.emptyList();
            }
        }

        // Handle the case where one or both of the input VCs are not biallelic.
        final List<Pair<Integer,Integer>> result = new ArrayList<>();

        // First split (and trim) all variant contexts into biallelics.  We are only going ot be interested in the alleles.
        final List<VariantContext> splitVariants1 = simpleSplitIntoBiallelics(variant1);
        final List<VariantContext> splitVariants2 = simpleSplitIntoBiallelics(variant2);

        // Second, match on ref and alt.  If match occurs add it to the output list.
        for (int i = 0; i < splitVariants1.size(); i++) {
            for (int j = 0; j < splitVariants2.size(); j++) {
                final VariantContext splitVariant1 = splitVariants1.get(i);
                final VariantContext splitVariant2 = splitVariants2.get(j);
                if (splitVariant1.getAlternateAllele(0).equals(splitVariant2.getAlternateAllele(0))
                        && splitVariant1.getReference().equals(splitVariant2.getReference())) {
                    result.add(Pair.of(i,j));
                }
            }
        }

        return result;
    }

    /**
     * TODO: Docs, since this requires some explanation.
     * Ignores the INFO field and genotype fields.  These will not be present.  This method is trying to be a bit fast.
     * @param vc
     * @return
     */
    private List<VariantContext> simpleSplitIntoBiallelics(final VariantContext vc) {

        final List<VariantContext> result = new ArrayList<>();

        for (final Allele allele : vc.getAlternateAlleles()) {
            result.add(
                    isBiallelic(vc) ? vc :
                    GATKVariantContextUtils.trimAlleles(
                        new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(), vc.getEnd(),
                            Arrays.asList(vc.getReference(), allele))
                            .make(), true, true)
            );
        }

        return result;
    }

    @Override
    /**
     * {@inheritDoc}
     * {@link VcfFuncotationFactory} can be used with or without Gencode annotations.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        return createFuncotationsOnVariant(variant, referenceContext, featureList);
    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private List<Funcotation> createDefaultFuncotationsOnVariantHelper( final VariantContext variant, final ReferenceContext referenceContext, final Set<Allele> annotatedAltAlleles  ) {

        final List<Funcotation> funcotationList = new ArrayList<>();

        if ( supportedFieldNames.size() != 0 ) {

            final List<Allele> alternateAlleles = variant.getAlternateAlleles();

            for ( final Allele altAllele : alternateAlleles ) {
                if ( !annotatedAltAlleles.contains(altAllele) ) {
                    funcotationList.add(TableFuncotation.create(supportedFieldNamesAndDefaults, altAllele, name, supportedFieldMetadata));
                }
            }
        }

        return funcotationList;
    }

    /**
     * Populates {@link VcfFuncotationFactory#supportedFieldNames} and {@link VcfFuncotationFactory#supportedFieldNamesAndDefaults}.
     */
    private void populateSupportedFieldNamesFromVcfFile() {
        final VCFFileReader reader = new VCFFileReader(sourceFilePath.toFile());
        final VCFHeader header = reader.getFileHeader();

        final List<String> infoLineKeys = new ArrayList<>();
        final Map<String, Boolean> infoFieldFlagMap = new HashMap<>();

        // Get our list of keys and sort them:
        for ( final VCFInfoHeaderLine infoLine : header.getInfoHeaderLines() ) {
            infoLineKeys.add(infoLine.getID());
            infoFieldFlagMap.put(infoLine.getID(), infoLine.getType() == VCFHeaderLineType.Flag);
        }
        infoLineKeys.sort(Comparator.naturalOrder());

        // Add our sorted names to the supported list:
        for ( final String key : infoLineKeys ) {
            if ( infoFieldFlagMap.get(key) ) {
                supportedFieldNamesAndDefaults.put(createFinalFieldName(name, key), "false" );
            }
            else {
                supportedFieldNamesAndDefaults.put(createFinalFieldName(name, key), "" );
            }
            supportedFieldNames.add(createFinalFieldName(name, key));
        }
    }

    @VisibleForTesting
    static String createFinalFieldName(final String funcotationFactoryName, final String fieldName) {
        return funcotationFactoryName + "_" + fieldName;
    }

    private Triple<VariantContext, ReferenceContext, List<Feature>> createCacheKey(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        return Triple.of(variant, referenceContext, featureList);
    }

    @Override
    public void close() {
        super.close();
        logger.info(getName() + " " + getVersion() + " cache hits/total: " + cacheHits + "/" + (cacheMisses + cacheHits));
    }

    //==================================================================================================================
    // Helper Data Types:

    // Modifed from https://docs.oracle.com/javase/7/docs/api/java/util/LinkedHashMap.html#removeEldestEntry(java.util.Map.Entry)
    private class LRUCache<K, V> extends LinkedHashMap<K, V> {
        static final long serialVersionUID = 55337L;
        private static final int MAX_ENTRIES = 20;
        public LRUCache() {
            super(MAX_ENTRIES);
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
            return size() >= MAX_ENTRIES;
        }
    }
}
