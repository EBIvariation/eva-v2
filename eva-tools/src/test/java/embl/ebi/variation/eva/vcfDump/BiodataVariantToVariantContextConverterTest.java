/*
 *
 *  * Copyright 2016 EMBL - European Bioinformatics Institute
 *  *
 *  * Licensed under the Apache License, Version 2.0 (the "License");
 *  * you may not use this file except in compliance with the License.
 *  * You may obtain a copy of the License at
 *  *
 *  *      http://www.apache.org/licenses/LICENSE-2.0
 *  *
 *  * Unless required by applicable law or agreed to in writing, software
 *  * distributed under the License is distributed on an "AS IS" BASIS,
 *  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  * See the License for the specific language governing permissions and
 *  * limitations under the License.
 *
 */

package embl.ebi.variation.eva.vcfdump;

import embl.ebi.variation.eva.vcfdump.cellbasewsclient.CellbaseWSClient;
import embl.ebi.variation.eva.vcfdump.exception.CellbaseSequenceDownloadError;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.After;
import org.junit.BeforeClass;
import org.junit.Test;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.*;
import org.opencb.opencga.lib.common.Config;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

/**
 * Created by pagarcia on 27/06/2016.
 */
public class BiodataVariantToVariantContextConverterTest {

    private static VariantFactory variantFactory;
    private static final String CHR_1 = "1";
    private static final String STUDY_1 = "study_1";
    private static ArrayList<String> s1s6SampleList;
    private static CellbaseWSClient cellbaseWSClientStub;
    private static Map<String, Map<String, String>> noSampleNamesConflictSampleNameCorrections = null;


    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        variantFactory = new VariantVcfFactory();
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        // cellbase client stub that returns some recorded sequences
        cellbaseWSClientStub = new CellbaseWSClient("anSpecie") {
            @Override
            public String getSequence(Region region) throws IOException {
                if (region.getChromosome().equals(CHR_1) && region.getStart() == 1099 && region.getEnd() == 1099) {
                    return "T";
                } else if (region.getChromosome().equals(CHR_1) && region.getStart() == 1089 && region.getEnd() == 1110) {
                    return "AATGTCTGCATAGCACGTTGAC";
                }
                return null;
            }
        };

        // example samples list
        s1s6SampleList = new ArrayList<>();
        for (int i = 1; i <= 6; i++) {
            s1s6SampleList.add("s" + i);
        }

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void singleStudySNV() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1000", "id", "C", "A", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        // TODO: cellbaseClient and region are null because we don't need them. Change API?
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), null, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "C", "A", new HashSet<>(s1s6SampleList));
    }

    @Test
    public void singleStudyInsertionRegionNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "", "G", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        // TODO: test with region
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1100, 1100, "T", "TG", new HashSet<>(s1s6SampleList));
    }

    @Test
    public void singleStudyDeletionRegionNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "A", "", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        // TODO: test with region
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1100, 1101, "TA", "T", new HashSet<>(s1s6SampleList));
    }

    @Test
    public void singleStudyInsertionRegionNotNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "", "G", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        // TODO: test with region
        VariantContext variantContext = variantConverter.transform(variants.get(0), new Region("1", 1090, 1110));
        checkVariantContext(variantContext, CHR_1, 1100, 1100, "T", "TG", new HashSet<>(s1s6SampleList));
    }

    @Test
    public void singleStudyDeletionRegionNotNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "A", "", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        // TODO: test with region
        VariantContext variantContext = variantConverter.transform(variants.get(0),new Region("1", 1090, 1110));
        checkVariantContext(variantContext, CHR_1, 1100, 1101, "TA", "T", new HashSet<>(s1s6SampleList));
    }

    @Test(expected=IllegalArgumentException.class)
    public void singleStudyInsertionMissingCellbase() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "", "G", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), null, noSampleNamesConflictSampleNameCorrections);
        // should thrown IllegalArgumentException:
        variantConverter.transform(variants.get(0), null);
        // TODO: assertions
    }

    @Test
    public void singleStudyMultiAllelicVariant() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1000", "id", "C", "A,T", "100", "PASS", ".", "GT", "0|0", "0|2", "0|1", "1|1", "1|2", "2|2");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(2, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variant1Context = variantConverter.transform(variants.get(0), null);
        VariantContext variant2Context = variantConverter.transform(variants.get(1), null);

        // TODO: assertions
    }

    // TODO: multiallelic indel test: String multiallelicIndelLine = "1\t1000\tid\tN\tNA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";

    @Test
    public void twoStudiesNoConflictingNamesSingleVariant() throws CellbaseSequenceDownloadError {
        // create test variant, with two studies and samples with not conflicting names
        Variant variant = new Variant(CHR_1, 1000, 1000, "T", "G");

        VariantSource source1 = new VariantSource("testFile1", "file_1", "study_1", "testStudy1");
        source1.setSamples(Arrays.asList("SX_1", "SX_2", "SX_3", "SX_4"));
        VariantSourceEntry study1Entry = new VariantSourceEntry("file_1", "study_1");
        addGenotype(study1Entry, "SX_1", "0|0");
        addGenotype(study1Entry, "SX_2", "0|1");
        addGenotype(study1Entry, "SX_3", "0|1");
        addGenotype(study1Entry, "SX_4", "0|0");
        variant.addSourceEntry(study1Entry);

        VariantSource source2 = new VariantSource("testFile2", "file_2", "study_2", "testStudy2");
        source1.setSamples(Arrays.asList("SY_1", "SY_2", "SY_3", "SY_4", "SY_5", "SY_6"));
        VariantSourceEntry study2Entry = new VariantSourceEntry("file_2", "study_2");
        addGenotype(study2Entry, "SY_1", "0|0");
        addGenotype(study2Entry, "SY_2", "1|0");
        addGenotype(study2Entry, "SY_3", "1|1");
        addGenotype(study2Entry, "SY_4", "1|1");
        addGenotype(study2Entry, "SY_5", "0|0");
        addGenotype(study2Entry, "SY_6", "1|0");
        variant.addSourceEntry(study2Entry);

        // transform variant
        BiodataVariantToVariantContextConverter variantConverter = new BiodataVariantToVariantContextConverter(Arrays.asList(source1, source2), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variant, null);

        // check transformed variant
        Set<String> sampleNames = new HashSet<>(source1.getSamples());
        sampleNames.addAll(source2.getSamples());
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "T", "G", sampleNames);
    }

    @Test
    public void twoStudiesConflictingNamesSingleVariant() throws CellbaseSequenceDownloadError {
        // create test variant, with two studies and samples with not conflicting names
        Variant variant = new Variant(CHR_1, 1000, 1000, "T", "G");

        // studies and samples names
        String study1 = "study_1";
        String study2 = "study_2";
        String sampleX1 = "SX_1";
        String sampleX2 = "SX_2";
        String sampleX3 = "SX_3";
        String sampleX4 = "SX_4";
        String sampleX5 = "SX_5";
        String sampleX6 = "SX_6";

        // sample name corrections map (usually created by VariantExporter)
        Map<String, Map<String, String>> sampleNamesCorrections = new HashMap<>();
        Map<String, String> study1SampleNameCorrections = new HashMap<>();
        study1SampleNameCorrections.put(sampleX1, study1 + "_" + sampleX1);
        study1SampleNameCorrections.put(sampleX2, study1 + "_" + sampleX2);
        study1SampleNameCorrections.put(sampleX3, study1 + "_" + sampleX3);
        study1SampleNameCorrections.put(sampleX4, study1 + "_" + sampleX4);
        sampleNamesCorrections.put(study1, study1SampleNameCorrections);
        Map<String, String> study2SampleNameCorrections = new HashMap<>();
        study2SampleNameCorrections.put(sampleX1, study2 + "_" + sampleX1);
        study2SampleNameCorrections.put(sampleX2, study2 + "_" + sampleX2);
        study2SampleNameCorrections.put(sampleX3, study2 + "_" + sampleX3);
        study2SampleNameCorrections.put(sampleX4, study2 + "_" + sampleX4);
        study2SampleNameCorrections.put(sampleX5, study2 + "_" + sampleX5);
        study2SampleNameCorrections.put(sampleX6, study2 + "_" + sampleX6);
        sampleNamesCorrections.put(study2, study2SampleNameCorrections);

        // variant sources
        VariantSource source1 = new VariantSource("testFile1", "file_1", study1, "testStudy1");
        source1.setSamples(Arrays.asList(sampleX1, sampleX2, sampleX3, sampleX4));
        VariantSourceEntry study1Entry = new VariantSourceEntry("file_1", study1);
        addGenotype(study1Entry, sampleX1, "0|0");
        addGenotype(study1Entry, sampleX2, "0|1");
        addGenotype(study1Entry, sampleX3, "0|1");
        addGenotype(study1Entry, sampleX4, "0|0");
        variant.addSourceEntry(study1Entry);
        VariantSource source2 = new VariantSource("testFile2", "file_2", study2, "testStudy2");
        source2.setSamples(Arrays.asList(sampleX1, sampleX2, sampleX3, sampleX4, sampleX5, sampleX6));
        VariantSourceEntry study2Entry = new VariantSourceEntry("file_2", study2);
        addGenotype(study2Entry, sampleX1, "0|0");
        addGenotype(study2Entry, sampleX2, "1|0");
        addGenotype(study2Entry, sampleX3, "1|1");
        addGenotype(study2Entry, sampleX4, "1|1");
        addGenotype(study2Entry, sampleX5, "0|0");
        addGenotype(study2Entry, sampleX6, "1|0");
        variant.addSourceEntry(study2Entry);

        // transform variant
        BiodataVariantToVariantContextConverter variantConverter = new BiodataVariantToVariantContextConverter(Arrays.asList(source1, source2), cellbaseWSClientStub, sampleNamesCorrections);
        VariantContext variantContext = variantConverter.transform(variant, null);

        // check transformed variant
        Set<String> sampleNames = source1.getSamples().stream().map(s -> study1Entry.getStudyId() + "_" + s).collect(Collectors.toSet());
        sampleNames.addAll(source2.getSamples().stream().map(s -> study2Entry.getStudyId() + "_" + s).collect(Collectors.toSet()));
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "T", "G", sampleNames);
    }

    // TODO: test with alleles -1|-1

    private void addGenotype(VariantSourceEntry sourceEntry, String sampleName, String genotype) {
        Map<String, String> sampleData = new HashMap<>();
        sampleData.put("GT", genotype);
        sourceEntry.addSampleData(sampleName, sampleData);
    }


    // TODO: check if those test cases are already covered by the preceding tests, and remove if so (or create new test with the test cases not covered)
//    @Test(expected=IllegalArgumentException.class)
//    public void singleStudyVariantTransformMissingCellbase() throws Exception {
//        String studyId = "studyId";
//        final VariantSource variantSource = createTestVariantSource(studyId);
//
//        List<String> studies = Collections.singletonList(studyId);
//
//        // TODO: split the test into different ones
//        // test multiallelic. these first conversions should NOT fail, as they doesn't need the src
//        String multiallelicLine = "1\t1000\tid\tC\tA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
//        List<Variant> variants = variantFactory.create(variantSource, multiallelicLine);
//        assertEquals(2, variants.size());
//        removeSrc(variants);    // <---- this is the key point of the test
//
//        VariantExporter variantExporter = new VariantExporter(null);
//        Map<String, VariantContext> variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(0), studies, null);
//
//        List<String> alleles = Arrays.asList("C", "A", ".");
//        assertEqualGenotypes(variants.get(0), variantContext.get(studyId), alleles);
//
//        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), studies, null);
//        alleles = Arrays.asList("C", "T", ".");
//        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);
//
//
//        // test multiallelic + indel
//        String multiallelicIndelLine = "1\t1000\tid\tC\tCA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
//        variants = variantFactory.create(variantSource, multiallelicIndelLine);
//        assertEquals(2, variants.size());
//        removeSrc(variants);    // <---- this is the key point of the test
//
//        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), studies, null);
//        alleles = Arrays.asList("C", "T", ".");
//        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);
//
//        // TODO: is the exception just thrown in the last line?
//        // the next exception is the only one that should throw
////        thrown.expect(IllegalArgumentException.class);
//        variantExporter.convertBiodataVariantToVariantContext(variants.get(0), studies, null);
//
//    }

    private void checkVariantContext(VariantContext variantContext, String chromosome, int start, int end, String ref, String alt, Set<String> sampleNames) {
        assertEquals(chromosome, variantContext.getContig());
        assertEquals(start, variantContext.getStart());
        assertEquals(end, variantContext.getEnd());
        assertEquals(Allele.create(ref, true), variantContext.getReference());
        assertEquals(Collections.singletonList(Allele.create(alt, false)), variantContext.getAlternateAlleles());
        assertTrue(variantContext.emptyID());
        // TODO: QUAL
//        assertEquals(100, variantContext.getPhredScaledQual(), 0);
        assertEquals(Collections.singleton(BiodataVariantToVariantContextConverter.FILTER_PASS), variantContext.getFilters());
        assertEquals(0, variantContext.getCommonInfo().getAttributes().size());
        assertEquals(sampleNames, variantContext.getSampleNames());
        // TODO check genotypes. use assertEqualGenotypes?
        // TODO: move outside?
//        checkGenotypes(variantContext);
    }

    private void assertEqualGenotypes(Variant variant, VariantContext variantContext, List<String> alleles) {
        for (Map.Entry<String, Map<String, String>> data : variant.getSourceEntries().values().iterator().next().getSamplesData().entrySet()) {
            Genotype genotype = variantContext.getGenotype(data.getKey());
            String gt = data.getValue().get("GT");
            org.opencb.biodata.models.feature.Genotype biodataGenotype = new org.opencb.biodata.models.feature.Genotype(gt, alleles.get(0), alleles.get(1));
            assertEquals(Allele.create(alleles.get(biodataGenotype.getAllele(0)), biodataGenotype.isAlleleRef(0)),
                    genotype.getAllele(0));
            assertEquals(Allele.create(alleles.get(biodataGenotype.getAllele(1)), biodataGenotype.isAlleleRef(1)),
                    genotype.getAllele(1));
        }
    }

    private VariantSource createTestVariantSource(String studyId) {
        return createTestVariantSource(studyId, s1s6SampleList);
    }

    private VariantSource createTestVariantSource(String studyId, List<String> sampleList) {
        final VariantSource variantSource = new VariantSource("name", "fileId", studyId, "studyName");
        variantSource.setSamples(sampleList);
        return variantSource;
    }


    private void removeSrc(List<Variant> variants) {
        for (Variant variant : variants) {
            for (VariantSourceEntry variantSourceEntry : variant.getSourceEntries().values()) {
                variantSourceEntry.getAttributes().remove("src");
            }
        }
    }

}