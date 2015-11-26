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
package embl.ebi.variation.eva.vcfDump;

import com.mongodb.DB;
import com.mongodb.MongoClient;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;
import org.opencb.biodata.models.variant.*;
import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URI;
import java.net.UnknownHostException;
import java.util.*;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by jmmut on 2015-10-29.
 * 
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporterTest {
    private static final String DB_NAME = "VariantExporterTest";

    private static final Logger logger = LoggerFactory.getLogger(VariantExporterTest.class);

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    @Test
    public void testVcfHtsExport() throws Exception {
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        QueryOptions query = new QueryOptions();
//        List<String> files = Arrays.asList("5");
        List<String> files = Arrays.asList("5", "6");
        List<String> studies = Collections.singletonList("7");
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);
        String outputDir = "/tmp/";

        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(DB_NAME, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");


        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, "hsapiens");

        VariantExporter variantExporter = new VariantExporter(cellBaseClient);
        List<String> outputFiles = variantExporter.VcfHtsExport(iterator, outputDir, variantSourceDBAdaptor, query);

        ////// checks 
        assertEquals(studies.size(), outputFiles.size());
        assertEquals(0, variantExporter.getFailedVariants());   // test file should not have failed variants

        iterator = variantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFiles, iterator);

        for (String outputFile : outputFiles) {
            boolean delete = new File(outputFile).delete();
            assertTrue(delete);
        }
    }
    
    @Test
    public void testVcfHtsExportSeveralStudies() throws Exception {
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        QueryOptions query = new QueryOptions();
//        List<String> files = Arrays.asList("5");
        List<String> files = Arrays.asList("5", "6");
        List<String> studies = Arrays.asList("7", "8");
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);
        String outputDir = "/tmp/";


        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(DB_NAME, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");


        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, "hsapiens");


        VariantExporter variantExporter = new VariantExporter(cellBaseClient);
        List<String> outputFiles = variantExporter.VcfHtsExport(iterator, outputDir, variantSourceDBAdaptor, query);

        ////////// checks

        assertEquals(studies.size(), outputFiles.size());
        assertEquals(0, variantExporter.getFailedVariants());
        
        // for study 7
        query.put(VariantDBAdaptor.STUDIES, Collections.singletonList("7"));
        iterator = variantDBAdaptor.iterator(query);
        assertEquals(countRows(iterator), countLines(outputFiles.get(0)));
        
        // for study 8
        query.put(VariantDBAdaptor.STUDIES, Collections.singletonList("8"));
        iterator = variantDBAdaptor.iterator(query);
        assertEquals(countRows(iterator), countLines(outputFiles.get(1)));
        

        for (String outputFile : outputFiles) {
            boolean delete = new File(outputFile).delete();
            assertTrue(delete);
        }
    }
    
    @Test
    public void testFilter() throws Exception {
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        QueryOptions query = new QueryOptions();
//        List<String> files = Arrays.asList("5");
        List<String> files = Arrays.asList("5");
        List<String> studies = Arrays.asList("7");
        String outputDir = "/tmp/";

        // tell all variables to filter with
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);
        query.put(VariantDBAdaptor.REGION, "20:61000-69000");
        query.put(VariantDBAdaptor.REFERENCE, "A");
//        query.put(VariantDBAdaptor., "A");

        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(DB_NAME, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");


        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, "hsapiens");


        VariantExporter variantExporter = new VariantExporter(cellBaseClient);
        List<String> outputFiles = variantExporter.VcfHtsExport(iterator, outputDir, variantSourceDBAdaptor, query);
        
        ////////// checks

        assertEquals(studies.size(), outputFiles.size());
        assertEquals(0, variantExporter.getFailedVariants());   // test file should not have failed variants

        iterator = variantDBAdaptor.iterator(query);
        assertEqualLinesFilesAndDB(outputFiles, iterator);

        for (String outputFile : outputFiles) {
            boolean delete = new File(outputFile).delete();
            assertTrue(delete);
        }
    }

    @Test
    public void testMissingStudy() throws Exception {
        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        QueryOptions query = new QueryOptions();
        List<String> files = Arrays.asList("5");
        List<String> studies = Arrays.asList("7", "9"); // study 9 doesn't exist
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);
        String outputDir = "/tmp/";

        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(DB_NAME, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");


        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, "hsapiens");


        VariantExporter variantExporter = new VariantExporter(cellBaseClient);

        thrown.expect(IllegalArgumentException.class);  // comment this line to see the actual exception, making the test fail
        variantExporter.VcfHtsExport(iterator, outputDir, variantSourceDBAdaptor, query);
    }

    @Test
    public void testMissingCellbase() throws Exception {
        final VariantSource variantSource = new VariantSource("name", "fileId", "studyId", "studyName");
        List<String> samples = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            samples.add("s"+i);
        }
        variantSource.setSamples(samples);
        VariantFactory factory = new VariantVcfFactory();
        Map<String, VariantSource> sources = Collections.singletonMap(variantSource.getStudyId(), variantSource);
        QueryOptions options = new QueryOptions("fileId", "fileId");
        String studyId = "studyId";
        options.add("studyId", studyId);
        List<Variant> variants;
        Map<String, VariantContext> variantContext;
        List<String> alleles;

        // test multiallelic. these first conversions should NOT fail, as they doesn't need the src
        String multiallelicLine = "1\t1000\tid\tC\tA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicLine);
        assertEquals(2, variants.size());
        removeSrc(variants);    // <---- this is the key point of the test

        VariantExporter variantExporter = new VariantExporter(null);
        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(0), sources);

        alleles = Arrays.asList("C", "A", ".");
        assertEqualGenotypes(variants.get(0), variantContext.get(studyId), alleles);

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), sources);
        alleles = Arrays.asList("C", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);

        // test multiallelic + indel
        String multiallelicIndelLine = "1\t1000\tid\tC\tCA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicIndelLine);
        assertEquals(2, variants.size());
        removeSrc(variants);    // <---- this is the key point of the test

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), sources);
        alleles = Arrays.asList("C", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);

        // the next conversion is the only one that should throw
        thrown.expect(IllegalArgumentException.class);
        variantExporter.convertBiodataVariantToVariantContext(variants.get(0), sources);

    }

    @Test
    public void testGetVariantContextFromVariant() throws Exception {
        final VariantSource variantSource = new VariantSource("name", "fileId", "studyId", "studyName");
        List<String> samples = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            samples.add("s" + i);
        }
        variantSource.setSamples(samples);
        VariantFactory factory = new VariantVcfFactory();
        Map<String, VariantSource> sources = Collections.singletonMap(variantSource.getStudyId(), variantSource);
        QueryOptions options = new QueryOptions("fileId", "fileId");
        String studyId = "studyId";
        options.add("studyId", studyId);
        List<Variant> variants;
        Map<String, VariantContext> variantContext;
        List<String> alleles;
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");


        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, "hsapiens");


        // test multiallelic
        String multiallelicLine = "1\t1000\tid\tC\tA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicLine);
        assertEquals(2, variants.size());

        VariantExporter variantExporter = new VariantExporter(cellBaseClient);
        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(0), sources);
        
        alleles = Arrays.asList("C", "A", ".");
        assertEqualGenotypes(variants.get(0), variantContext.get(studyId), alleles);

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), sources);
        alleles = Arrays.asList("C", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);


        // test indel
        String indelLine = "1\t1000\tid\tN\tNA\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|0\t0|1";
        variants = factory.create(variantSource, indelLine);

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(0), sources);
        alleles = Arrays.asList("N", "NA");
        assertEqualGenotypes(variants.get(0), variantContext.get(studyId), alleles);


        // test multiallelic + indel
        String multiallelicIndelLine = "1\t1000\tid\tN\tNA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicIndelLine);
        assertEquals(2, variants.size());

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(0), sources);
        alleles = Arrays.asList("N", "NA", ".");
        assertEqualGenotypes(variants.get(0), variantContext.get(studyId), alleles);

        variantContext = variantExporter.convertBiodataVariantToVariantContext(variants.get(1), sources);
        alleles = Arrays.asList("N", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext.get(studyId), alleles);

    }

    /**
     * Clears and populates the Mongo collection used during the tests.
     * @throws java.io.IOException
     * @throws java.lang.InterruptedException
     */
    @BeforeClass
    public static void setUpClass() throws IOException, InterruptedException {
        cleanDBs();
        fillDB();
    }

    /**
     * Clears and populates the Mongo collection used during the tests.
     * 
     * @throws UnknownHostException 
     */
    @AfterClass
    public static void tearDownClass() throws UnknownHostException {
        cleanDBs();
    }

    private static void cleanDBs() throws UnknownHostException {
        MongoClient mongoClient = new MongoClient("localhost");
        List<String> dbs = Arrays.asList(DB_NAME);
        for (String dbName : dbs) {
            DB db = mongoClient.getDB(dbName);
            db.dropDatabase();
        }
        mongoClient.close();
    }

    private static void fillDB() throws IOException, InterruptedException {
        String dump = VariantExporterTest.class.getResource("/dump/").getFile();
        logger.info("restoring DB from " + dump);
        Process exec = Runtime.getRuntime().exec("mongorestore " + dump);
        exec.waitFor();
        String line;
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(exec.getInputStream()));
        while ((line = bufferedReader.readLine()) != null) {
            logger.info("mongorestore output:" + line);
        }
        bufferedReader.close();
        bufferedReader = new BufferedReader(new InputStreamReader(exec.getErrorStream()));
        while ((line = bufferedReader.readLine()) != null) {
            logger.info("mongorestore errorOutput:" + line);
        }
        bufferedReader.close();

        logger.info("mongorestore exit value: " + exec.exitValue());
    }
    
    private void removeSrc(List<Variant> variants) {
        for (Variant variant : variants) {
            for (VariantSourceEntry variantSourceEntry : variant.getSourceEntries().values()) {
                variantSourceEntry.getAttributes().remove("src");
            }
        }
    }

    private void assertEqualLinesFilesAndDB(List<String> fileNames, VariantDBIterator iterator) throws IOException {
        long lines = 0;
        for (String fileName : fileNames) {
            lines += countLines(fileName);
        }

        // counting variants in the DB
        long variantRows = countRows(iterator);

        assertEquals(variantRows, lines);
    }

    private long countRows(Iterator<Variant> iterator) {
        int variantRows = 0;
        while(iterator.hasNext()) {
            iterator.next();
            variantRows++;
        }
        return variantRows;
    }

    private long countLines(String fileName) throws IOException {
        long lines;
        lines = 0;
        BufferedReader file = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
        String line;
        while ((line = file.readLine()) != null) {
            if (line.charAt(0) != '#') {
                lines++;
            }
        }
        file.close();
        return lines;
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

}
