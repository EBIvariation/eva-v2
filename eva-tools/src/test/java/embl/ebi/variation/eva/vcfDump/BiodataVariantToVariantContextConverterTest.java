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

    public static final String FILE_ID = "fileId";
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
        cellbaseWSClientStub = new CellbaseWSClient("aSpecies") {
            @Override
            public String getSequence(Region region) throws IOException {
                if (region.getChromosome().equals(CHR_1) && region.getStart() == 1100 && region.getEnd() == 1100) {
                    return "T";
                } else if (region.getChromosome().equals(CHR_1) && region.getStart() == 1000 && region.getEnd() == 1000) {
                    return "C";
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
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), null, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "C", "A", variants.get(0).getSourceEntries(), false);
    }

    @Test
    public void singleStudyInsertionRegionNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "T", "TG", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1100, 1100, "T", "TG", variants.get(0).getSourceEntries(), false);
    }

    @Test
    public void singleStudyDeletionRegionNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "TA", "T", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0), null);
        checkVariantContext(variantContext, CHR_1, 1100, 1101, "TA", "T", variants.get(0).getSourceEntries(), false);
    }

    @Test
    public void singleStudyInsertionRegionNotNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "T", "TG", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0), new Region("1", 1090, 1110));
        checkVariantContext(variantContext, CHR_1, 1100, 1100, "T", "TG", variants.get(0).getSourceEntries(), false);
    }

    @Test
    public void singleStudyDeletionRegionNotNull() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "TA", "T", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variantContext = variantConverter.transform(variants.get(0),new Region("1", 1090, 1110));
        checkVariantContext(variantContext, CHR_1, 1100, 1101, "TA", "T", variants.get(0).getSourceEntries(), false);
    }

    @Test(expected=IllegalStateException.class)
    public void singleStudyInsertionMissingCellbase() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1100", "id", "", "G", "100", "PASS", ".", "GT", "0|0", "0|0", "0|1", "1|1", "1|1", "0|1");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(1, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), null, noSampleNamesConflictSampleNameCorrections);
        // should thrown IllegalStateException:
        variantConverter.transform(variants.get(0), null);
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

        checkVariantContext(variant1Context, CHR_1, 1000, 1000, "C", "A", variants.get(0).getSourceEntries(), false);
        checkVariantContext(variant2Context, CHR_1, 1000, 1000, "C", "T", variants.get(1).getSourceEntries(), false);
    }

    @Test
    public void singleStudyMultiAllelicIndel() throws CellbaseSequenceDownloadError {
        // create SNV variant
        VariantSource variantSource = createTestVariantSource(STUDY_1);
        String variantLine =  String.join("\t", CHR_1, "1000", "id", "C", "CA,T", "100", "PASS", ".", "GT", "0|0", "0|2", "0|1", "1|1", "1|2", "2|2");
        List<Variant> variants = variantFactory.create(variantSource, variantLine);
        assertEquals(2, variants.size());

        // export variant
        BiodataVariantToVariantContextConverter variantConverter =
                new BiodataVariantToVariantContextConverter(Collections.singletonList(variantSource), cellbaseWSClientStub, noSampleNamesConflictSampleNameCorrections);
        VariantContext variant1Context = variantConverter.transform(variants.get(0), null);
        VariantContext variant2Context = variantConverter.transform(variants.get(1), null);

        checkVariantContext(variant1Context, CHR_1, 1000, 1000, "C", "CA", variants.get(0).getSourceEntries(), false);
        checkVariantContext(variant2Context, CHR_1, 1000, 1000, "C", "T", variants.get(1).getSourceEntries(), false);
    }

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
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "T", "G", variant.getSourceEntries(), false);
    }

    @Test
    public void twoStudiesConflictingNamesSingleVariant() throws CellbaseSequenceDownloadError {
        // create test variant, with two studies and samples with not conflicting names
        Variant variant = new Variant(CHR_1, 1000, 1000, "T", "G");

        // studies and samples names
        String study1 = "study_1";
        String study2 = "study_2";
        String file1 = "file_1";
        String file2 = "file_2";
        String sampleX1 = "SX_1";
        String sampleX2 = "SX_2";
        String sampleX3 = "SX_3";
        String sampleX4 = "SX_4";
        String sampleX5 = "SX_5";
        String sampleX6 = "SX_6";

        // sample name corrections map (usually created by VariantExporter)
        Map<String, Map<String, String>> sampleNamesCorrections = new HashMap<>();
        Map<String, String> file1SampleNameCorrections = new HashMap<>();
        file1SampleNameCorrections.put(sampleX1, file1 + "_" + sampleX1);
        file1SampleNameCorrections.put(sampleX2, file1 + "_" + sampleX2);
        file1SampleNameCorrections.put(sampleX3, file1 + "_" + sampleX3);
        file1SampleNameCorrections.put(sampleX4, file1 + "_" + sampleX4);
        sampleNamesCorrections.put(file1, file1SampleNameCorrections);
        Map<String, String> file2SampleNameCorrections = new HashMap<>();
        file2SampleNameCorrections.put(sampleX1, file2 + "_" + sampleX1);
        file2SampleNameCorrections.put(sampleX2, file2 + "_" + sampleX2);
        file2SampleNameCorrections.put(sampleX3, file2 + "_" + sampleX3);
        file2SampleNameCorrections.put(sampleX4, file2 + "_" + sampleX4);
        file2SampleNameCorrections.put(sampleX5, file2 + "_" + sampleX5);
        file2SampleNameCorrections.put(sampleX6, file2 + "_" + sampleX6);
        sampleNamesCorrections.put(file2, file2SampleNameCorrections);

        // variant sources
        VariantSource source1 = new VariantSource("testFile1", file1, study1, "testStudy1");
        source1.setSamples(Arrays.asList(sampleX1, sampleX2, sampleX3, sampleX4));
        VariantSourceEntry source1Entry = new VariantSourceEntry(file1, study1);
        addGenotype(source1Entry, sampleX1, "0|0");
        addGenotype(source1Entry, sampleX2, "0|1");
        addGenotype(source1Entry, sampleX3, "0|1");
        addGenotype(source1Entry, sampleX4, "0|0");
        variant.addSourceEntry(source1Entry);
        VariantSource source2 = new VariantSource("testFile2", file2, study2, "testStudy2");
        source2.setSamples(Arrays.asList(sampleX1, sampleX2, sampleX3, sampleX4, sampleX5, sampleX6));
        VariantSourceEntry source2Entry = new VariantSourceEntry(file2, study2);
        addGenotype(source2Entry, sampleX1, "0|0");
        addGenotype(source2Entry, sampleX2, "1|0");
        addGenotype(source2Entry, sampleX3, "1|1");
        addGenotype(source2Entry, sampleX4, "-1|-1");
        addGenotype(source2Entry, sampleX5, "0|0");
        addGenotype(source2Entry, sampleX6, "1|0");
        variant.addSourceEntry(source2Entry);

        // transform variant
        BiodataVariantToVariantContextConverter variantConverter = new BiodataVariantToVariantContextConverter(Arrays.asList(source1, source2), cellbaseWSClientStub, sampleNamesCorrections);
        VariantContext variantContext = variantConverter.transform(variant, null);

        // check transformed variant
        Set<String> sampleNames = source1.getSamples().stream().map(s -> source1Entry.getFileId() + "_" + s).collect(Collectors.toSet());
        sampleNames.addAll(source2.getSamples().stream().map(s -> source2Entry.getFileId() + "_" + s).collect(Collectors.toSet()));
        checkVariantContext(variantContext, CHR_1, 1000, 1000, "T", "G", variant.getSourceEntries(), true);
    }

    private void addGenotype(VariantSourceEntry sourceEntry, String sampleName, String genotype) {
        Map<String, String> sampleData = new HashMap<>();
        sampleData.put("GT", genotype);
        sourceEntry.addSampleData(sampleName, sampleData);
    }

    private void checkVariantContext(VariantContext variantContext, String chromosome, int start, int end, String ref, String alt,
                                     Map<String, VariantSourceEntry> sourceEntries, boolean sampleNameConflicts) {
        assertEquals(chromosome, variantContext.getContig());
        assertEquals(start, variantContext.getStart());
        assertEquals(end, variantContext.getEnd());
        assertEquals(Allele.create(ref, true), variantContext.getReference());
        assertEquals(Collections.singletonList(Allele.create(alt, false)), variantContext.getAlternateAlleles());
        assertTrue(variantContext.emptyID());
        assertTrue(variantContext.getFilters().isEmpty());
        assertEquals(0, variantContext.getCommonInfo().getAttributes().size());
        checkGenotypes(sourceEntries, variantContext, sampleNameConflicts);
    }

    private void checkGenotypes(Map<String, VariantSourceEntry> sourceEntries, VariantContext variantContext, boolean sampleNameConflicts) {
        // check that variantContext has the same number of samples than the input variant
        int inputVariantsSampleCount = sourceEntries.values().stream().mapToInt(variantSourceEntry -> variantSourceEntry.getSamplesData().size()).sum();
        assertEquals(inputVariantsSampleCount, variantContext.getSampleNames().size());

        for (Map.Entry<String, VariantSourceEntry> sourcesMapEntry : sourceEntries.entrySet()) {
            checkStudyGenotypes(sourcesMapEntry.getValue().getFileId(), sourcesMapEntry.getValue(), variantContext, sampleNameConflicts);
        }
    }

    private void checkStudyGenotypes(String fileId, VariantSourceEntry sourceEntry, VariantContext variantContext, boolean sampleNameConflicts) {
        for (String sample : sourceEntry.getSampleNames()) {
            String inputVariantSampleGenotype = sourceEntry.getSampleData(sample, "GT");
            String sampleNameInOutputVariant;
            if (sampleNameConflicts) {
                sampleNameInOutputVariant = fileId + "_" + sample;
            } else {
                sampleNameInOutputVariant = sample;
            }
            Genotype outputVariantSampleGenotype = variantContext.getGenotype(sampleNameInOutputVariant);
            compareAlleles(inputVariantSampleGenotype.charAt(0), outputVariantSampleGenotype.getAllele(0));
            compareAlleles(inputVariantSampleGenotype.charAt(2), outputVariantSampleGenotype.getAllele(1));
        }
    }

    private void compareAlleles(char inputAlleleIndex, Allele allele) {
        if (inputAlleleIndex == '0') {
            assertTrue(allele.isReference());
        } else {
            assertTrue(allele.isNonReference());
        }
    }

    private VariantSource createTestVariantSource(String studyId) {
        return createTestVariantSource(studyId, s1s6SampleList);
    }

    private VariantSource createTestVariantSource(String studyId, List<String> sampleList) {
        final VariantSource variantSource = new VariantSource("name", FILE_ID, studyId, "studyName");
        variantSource.setSamples(sampleList);
        return variantSource;
    }

}