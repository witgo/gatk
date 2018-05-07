package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class Mutect2FilteringEngineUnitTest extends BaseTest {
    Mutect2FilteringEngine engine;

    @BeforeSuite
    public void init(){
        engine = new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), "badTumor");
    }

    @DataProvider(name = "failingCase")
    public Object[][] makeFailingCaseData() {
        return new Object[][]{
                {Arrays.asList(0.01, 0.01, 0.05, 0.01, 0.3, 0.4), 0.0, 0.0},
        };
    }

    @DataProvider(name = "falsePositiveRateData")
    public Object[][] makeFalsePositiveRateData() {
        return new Object[][]{
                {Arrays.asList(0.01, 0.01, 0.05, 0.01, 0.3, 0.4), 0.1, 0.3},
                // If there are ties e.g. 0.9's below, set the threshold such that the expected false positive rate does not exceed the requested max false positive rate
                // Filtering variants with p > 0.9 gives you the expected false positive rate of 0.076 > 0.5, so the threshold must be 0.01
                {Arrays.asList(0.01, 0.01, 0.9, 0.9, 0.9), 0.05, 0.01},
                {Arrays.asList(0.01, 0.01, 0.01, 0.01, 0.01), 0.05, 1.0}, // The FPR never exceeds the max FPR so do not filter
                {Arrays.asList(0.99, 0.99, 0.99, 1.0, 1.0), 0.05, 0.0}, // Impossible to meet the FPR with this data so filter everything
                {Arrays.asList(0.01, 0.01, 0.05, 0.01, 0.3, 0.4), 0.0, 0.0},
        };
    }

    @Test(dataProvider = "failingCase")
    public void testCalculateThresholdForReadOrientationFilter(final List<Double> posteriors,
                                                               final double maxErrorRate,
                                                               final double expectedThreshold){
        Mutect2FilterStats stats = engine.calculateThresholdForReadOrientationFilter(posteriors, maxErrorRate);
        Assert.assertEquals(stats.getThreshold(), expectedThreshold);
    }

}