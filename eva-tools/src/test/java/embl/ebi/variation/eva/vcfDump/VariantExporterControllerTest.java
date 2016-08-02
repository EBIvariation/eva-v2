/*
 * Copyright 2015 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package embl.ebi.variation.eva.vcfdump;

import embl.ebi.variation.eva.vcfdump.cellbasewsclient.CellbaseWSClient;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.ws.rs.core.MultivaluedHashMap;
import javax.ws.rs.core.MultivaluedMap;
import java.io.*;
import java.net.URISyntaxException;
import java.net.UnknownHostException;
import java.util.*;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

/**
 * Created by pagarcia on 10/05/2016.
 */
public class VariantExporterControllerTest {

    private static final String DB_NAME = "VariantExporterTest";
    public static final String OUTPUT_DIR = "/tmp/";

    private static CellbaseWSClient cellBaseClient;
    private static VariantStorageManager variantStorageManager;
    private static VariantDBAdaptor variantDBAdaptor;
    private static VariantDBAdaptor cowVariantDBAdaptor;
    private static final Logger logger = LoggerFactory.getLogger(VariantExporterControllerTest.class);
    private static final MultivaluedMap<String, String> emptyFilter = new MultivaluedHashMap<>();
    private static List<String> testOutputFiles;


    @BeforeClass
    public static void setUpClass() throws IllegalAccessException, ClassNotFoundException, InstantiationException, URISyntaxException, StorageManagerException, IOException, InterruptedException {
        VariantExporterTestDB.cleanDBs();
        VariantExporterTestDB.fillDB();

        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");
        cellBaseClient = new CellbaseWSClient("hsapiens");
        logger.info("Using cellbase: " + cellBaseClient.getUrl() + " version " + cellBaseClient.getVersion());
        variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        variantDBAdaptor = variantStorageManager.getDBAdaptor(DB_NAME, null);
        cowVariantDBAdaptor = variantStorageManager.getDBAdaptor(VariantExporterTestDB.COW_TEST_DB_NAME, null);
        testOutputFiles = new ArrayList<>();
    }

    /**
     * Clears and populates the Mongo collection used during the tests.
     *
     * @throws UnknownHostException
     */
    @AfterClass
    public static void tearDownClass() throws UnknownHostException {
        VariantExporterTestDB.cleanDBs();
        testOutputFiles.forEach(f -> new File(f).delete());
    }

    @Test
    public void testGetQuery() throws Exception {
        List<String> studies = Arrays.asList("s1", "s2");
        List<String> files = Arrays.asList("f3", "f4", "f5");

        VariantExporterController controller = new VariantExporterController("hsapiens", DB_NAME, studies, files, OUTPUT_DIR, emptyFilter);

        // empty query parameters
        MultivaluedMap<String, String> emptyParameters = new MultivaluedHashMap<>();
        QueryOptions query = controller.getQuery(emptyParameters);
        assertEquals(2, query.size());
        assertEquals(studies, query.getAsStringList(VariantDBAdaptor.STUDIES));
        assertEquals(files, query.getAsStringList(VariantDBAdaptor.FILES));

        // some not accepted parameters
        MultivaluedMap<String, String> nonAcceptedParameters = new MultivaluedHashMap<>();
        nonAcceptedParameters.add("wrongParameter1", "1");
        nonAcceptedParameters.add("wrongParameter1", "1b");
        nonAcceptedParameters.add("wrongParameter2", "2");
        query = controller.getQuery(nonAcceptedParameters);
        assertEquals(2, query.size());
        assertEquals(studies, query.getAsStringList(VariantDBAdaptor.STUDIES));
        assertEquals(files, query.getAsStringList(VariantDBAdaptor.FILES));

        // some accepted parameters
        MultivaluedMap<String, String> acceptedParameters = new MultivaluedHashMap<>();
        String region1 = "1:1000-2000";
        String region2 = "1:2500-3000";
        acceptedParameters.add(VariantDBAdaptor.REGION, region1);
        acceptedParameters.add(VariantDBAdaptor.REGION, region2);
        String id = "rs1234";
        acceptedParameters.add(VariantDBAdaptor.ID, id);
        query = controller.getQuery(acceptedParameters);
        assertEquals(4, query.size());
        assertEquals(studies, query.getAsStringList(VariantDBAdaptor.STUDIES));
        assertEquals(files, query.getAsStringList(VariantDBAdaptor.FILES));
        assertEquals(Arrays.asList(region1, region2), query.getAsStringList(VariantDBAdaptor.REGION));
        assertEquals(id, query.getString(VariantDBAdaptor.ID));

        // mixed accepted and non accepted parameters
        MultivaluedMap<String, String> mixedParameters = new MultivaluedHashMap<>();
        mixedParameters.add(VariantDBAdaptor.REGION, region1);
        mixedParameters.add(VariantDBAdaptor.REGION, region2);
        mixedParameters.add(VariantDBAdaptor.ID, id);
        mixedParameters.add("wrongParameter1", "1");
        mixedParameters.add("wrongParameter1", "1b");
        mixedParameters.add("wrongParameter2", "2");
        query = controller.getQuery(mixedParameters);
        assertEquals(4, query.size());
        assertEquals(studies, query.getAsStringList(VariantDBAdaptor.STUDIES));
        assertEquals(files, query.getAsStringList(VariantDBAdaptor.FILES));
        assertEquals(Arrays.asList(region1, region2), query.getAsStringList(VariantDBAdaptor.REGION));
        assertEquals(id, query.getString(VariantDBAdaptor.ID));
    }

    @Test
    public void testVcfHtsExportOneStudy() throws ClassNotFoundException, StorageManagerException, URISyntaxException, InstantiationException, IllegalAccessException, IOException {
        String studyId = "PRJEB6119";
        List<String> studies = Collections.singletonList(studyId);

        VariantExporterController controller = new VariantExporterController("btaurus", VariantExporterTestDB.COW_TEST_DB_NAME, studies, OUTPUT_DIR, emptyFilter);
        controller.run();

        ////////// checks
        String outputFile = controller.getOuputFilePath();
        testOutputFiles.add(outputFile);
        assertEquals(0, controller.getFailedVariants());   // test file should not have failed variants
        QueryOptions query = getQuery(studies);
        VariantDBIterator iterator = cowVariantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFile, iterator);
        checkOrderInOutputFile(outputFile);
    }

    @Test
    public void testVcfHtsExportSeveralStudies() throws Exception {
        String study6119 = "PRJEB6119";
        String study7061 = "PRJEB7061";
        List<String> studies = Arrays.asList(study6119, study7061);

        VariantExporterController controller = new VariantExporterController("btaurus", VariantExporterTestDB.COW_TEST_DB_NAME, studies, OUTPUT_DIR, emptyFilter);
        controller.run();

        ////////// checks
        String outputFile = controller.getOuputFilePath();
        testOutputFiles.add(outputFile);
        assertEquals(0, controller.getFailedVariants());   // test file should not have failed variants
        QueryOptions query = getQuery(Arrays.asList(study6119, study7061));
        VariantDBIterator iterator = cowVariantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFile, iterator);
        checkOrderInOutputFile(outputFile);
    }

    @Test
    public void testConsequenceTypeFilter() throws Exception {
        String studyId = "7";
        List<String> studies = Collections.singletonList(studyId);

        MultivaluedMap<String, String> filter = new MultivaluedHashMap<>();
        filter.putSingle(VariantDBAdaptor.ANNOT_CONSEQUENCE_TYPE, "1627");
        VariantExporterController controller = new VariantExporterController("hsapiens", DB_NAME, studies, OUTPUT_DIR, filter);
        controller.run();

        ////////// checks
        String outputFile = controller.getOuputFilePath();
        testOutputFiles.add(outputFile);
        assertEquals(0, controller.getFailedVariants());   // test file should not have failed variants
        QueryOptions query = controller.getQuery(filter);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFile, iterator);
        checkOrderInOutputFile(outputFile);
    }

    @Test
    public void testConsequenceTypeAndRegionFilter() throws Exception {
        String studyId = "7";
        List<String> studies = Collections.singletonList(studyId);

        MultivaluedMap<String, String> filter = new MultivaluedHashMap<>();
        filter.putSingle(VariantDBAdaptor.REGION, "20:60000-61000");
        filter.putSingle(VariantDBAdaptor.ANNOT_CONSEQUENCE_TYPE, "1627");
        VariantExporterController controller = new VariantExporterController("hsapiens", DB_NAME, studies, OUTPUT_DIR, filter);
        controller.run();

        ////////// checks
        String outputFile = controller.getOuputFilePath();
        testOutputFiles.add(outputFile);
        assertEquals(0, controller.getFailedVariants());   // test file should not have failed variants
        QueryOptions query = controller.getQuery(filter);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFile, iterator);
        checkOrderInOutputFile(outputFile);
    }



    @Test
    public void testFilterUsingIntersectingRegions() throws Exception {
        String studyId = "7";
        List<String> studies = Collections.singletonList(studyId);

        // tell all variables to filter with
        MultivaluedMap<String, String> filter = new MultivaluedHashMap<>();
        filter.put(VariantDBAdaptor.REGION, Arrays.asList("20:61000-66000", "20:63000-69000"));

        VariantExporterController controller = new VariantExporterController("hsapiens", DB_NAME, studies, OUTPUT_DIR, filter);
        controller.run();

        ////////// checks
        String outputFile = controller.getOuputFilePath();
        testOutputFiles.add(outputFile);
        assertEquals(0, controller.getFailedVariants());   // test file should not have failed variants

        QueryOptions query = getQuery(studies);
        query.put(VariantDBAdaptor.REGION, String.join(",", filter.get(VariantDBAdaptor.REGION)));
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);

        assertEqualLinesFilesAndDB(outputFile, iterator);

        checkOrderInOutputFile(outputFile);
    }

    private void checkOrderInOutputFile(String outputFile) {
        assertVcfOrderedByCoordinate(outputFile);
        this.logger.info("Deleting output temp file {}", outputFile);
        boolean delete = new File(outputFile).delete();
        assertTrue(delete);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testMissingStudy() throws Exception {
        List<String> studies = Arrays.asList("7", "9"); // study 9 doesn't exist

        VariantExporterController controller = new VariantExporterController("hsapiens", DB_NAME, studies, OUTPUT_DIR, emptyFilter);

        controller.run();
    }

    @Test(expected = IllegalArgumentException.class)
    public void nullSpeciesThrowsIllegalArgumentException() throws Exception {
        List<String> studies = Collections.singletonList("8");
        new VariantExporterController(null, DB_NAME, studies, OUTPUT_DIR, emptyFilter);
    }

    @Test(expected = IllegalArgumentException.class)
    public void nullDbnameThrowsIllegalArgumentException() throws Exception {
        List<String> studies = Collections.singletonList("8");
        new VariantExporterController("hsapiens", null, studies, OUTPUT_DIR, emptyFilter);
    }

    @Test(expected = IllegalArgumentException.class)
    public void emtpyStudiesThrowsIllegalArgumentException() throws Exception {
        new VariantExporterController("hsapiens", DB_NAME, Collections.EMPTY_LIST, OUTPUT_DIR, emptyFilter);
    }

    @Test(expected = IllegalArgumentException.class)
    public void nullOutputDirThrowsIllegalArgumentException() throws Exception {
        List<String> studies = Collections.singletonList("8");
        new VariantExporterController("hsapiens", DB_NAME, studies, null, emptyFilter);
    }

    private void assertEqualLinesFilesAndDB(String fileName, VariantDBIterator iterator) throws IOException {
        List<Variant> exportedVariants = getVariantsFromOutputFile(fileName);

        // counting variants in the DB
        List<Variant> variantsInDb = getVariantsFromDB(iterator);
        assertEquals(variantsInDb.size(), exportedVariants.size());
    }

    private List<Variant> getVariantsFromDB(Iterator<Variant> iterator) {
        List<Variant> variants = new ArrayList<>();
        while(iterator.hasNext()) {
            variants.add(iterator.next());
        }
        return variants;
    }

    private List<Variant> getVariantsFromOutputFile(String fileName) throws IOException {
        List<Variant> variantIds = new ArrayList<>();
        BufferedReader file = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
        String line;
        while ((line = file.readLine()) != null) {
            if (line.charAt(0) != '#') {
                String[] fields = line.split("\t", 6);
                Variant variant = new Variant(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[1]), fields[3], fields[4]);
                //variant.setEnd(variant.getStart() + variant.getLength() - 1);
                if (variant.getAlternate().substring(0, 1).equals(variant.getReference().substring(0, 1))) {
                    variant.setAlternate(variant.getAlternate().substring(1));
                    variant.setReference(variant.getReference().substring(1));
                }
                variantIds.add(variant);
            }
        }
        file.close();
        return variantIds;
    }

    private void assertVcfOrderedByCoordinate(String fileName) {
        logger.info("Checking that {} is sorted by coordinate", fileName);
        Set<String> finishedContigs = new HashSet<>();
        VCFFileReader vcfReader = new VCFFileReader(new File(fileName), false);
        String lastContig = null;
        int previousStart = -1;

        for (VariantContext variant : vcfReader) {
            // check chromosome
            if (lastContig == null || !variant.getContig().equals(lastContig)) {
                if (lastContig != null) {
                    finishedContigs.add(lastContig);
                }
                lastContig = variant.getContig();
                assertFalse("The variants should by grouped by contig in the vcf output", finishedContigs.contains(lastContig));
                previousStart = -1;
            }
            assertTrue("The vcf is not sorted by coordinate: " + variant.getContig() + ":" + variant.getStart() + ":" +
                    variant.getReference() + "->" + variant.getAlternateAlleles() + "; Previous variant start: " + previousStart,
                    variant.getStart() >= previousStart);
            previousStart = variant.getStart();
        }

    }

    private QueryOptions getQuery(List<String> studies) {
        QueryOptions query = new QueryOptions();
        query.put(VariantDBAdaptor.STUDIES, studies);
        return query;
    }
}