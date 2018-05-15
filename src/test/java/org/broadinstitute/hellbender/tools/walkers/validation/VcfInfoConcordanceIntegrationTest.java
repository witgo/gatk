package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class VcfInfoConcordanceIntegrationTest extends CommandLineProgramTest {
    final double epsilon = 1e-3;

    @Test
    public void testSimple() throws Exception {
        final String inputVcf = largeFileTestDir + "VQSR/expected/chr20_tiny_tf_python_gpu2.vcf";
        final File summary = createTempFile("summary", ".txt");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, inputVcf)
                .addArgument(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, inputVcf)
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("eval-info-key", "CNN_2D")
                .addArgument("truth-info-key", "SMALL_2D_TF_MODEL")
                .addArgument("epsilon", "0.01")
                .addArgument(VcfInfoConcordance.SUMMARY_LONG_NAME, summary.toString());
        runCommandLine(argsBuilder);

        InfoConcordanceRecord.InfoConcordanceReader reader = new InfoConcordanceRecord.InfoConcordanceReader(summary);
        InfoConcordanceRecord snpRecord = reader.readRecord();
        InfoConcordanceRecord indelRecord = reader.readRecord();

        Assert.assertEquals(snpRecord.getVariantType(), VariantContext.Type.SNP);
        Assert.assertEquals(indelRecord.getVariantType(), VariantContext.Type.INDEL);

        // numbers verified by manual inspection
        Assert.assertEquals(snpRecord.getMean(), 0.01289, epsilon);
        Assert.assertEquals(snpRecord.getStd(), 0.179343, epsilon);
        Assert.assertEquals(indelRecord.getMean(), 0.0102386, epsilon);
        Assert.assertEquals(indelRecord.getStd(), 0.038635, epsilon);
    }
}
