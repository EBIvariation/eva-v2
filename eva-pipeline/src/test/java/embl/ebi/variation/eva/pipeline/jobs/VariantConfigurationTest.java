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
package embl.ebi.variation.eva.pipeline.jobs;

import java.io.*;

import embl.ebi.variation.eva.pipeline.steps.VariantsLoad;
import embl.ebi.variation.eva.pipeline.steps.VariantsStatsCreate;
import embl.ebi.variation.eva.pipeline.steps.VariantsStatsLoad;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.batch.core.*;
import org.springframework.batch.core.launch.JobLauncher;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import java.net.UnknownHostException;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static embl.ebi.variation.eva.pipeline.jobs.JobTestUtils.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static embl.ebi.variation.eva.pipeline.listeners.JobParametersListener.*;

/**
 * Created by jmmut on 2015-10-14.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(classes = {VariantConfiguration.class})
public class VariantConfigurationTest {

    public static final String FILE_20 = "/small20.vcf.gz";
    public static final String FILE_22 = "/small22.vcf.gz";
    public static final String FILE_WRONG_NO_ALT = "/wrong_no_alt.vcf.gz";

    private static final Logger logger = LoggerFactory.getLogger(VariantConfigurationTest.class);

    // iterable doing an enum. Does it worth it?
    private static final String VALID_TRANSFORM = "VariantConfigurationTest_vt";
    private static final String INVALID_TRANSFORM = "VariantConfigurationTest_it";
    private static final String VALID_LOAD = "VariantConfigurationTest_vl";
//    private static final String INVALID_LOAD = "invalidLoad";
    private static final String VALID_CREATE_STATS = "VariantConfigurationTest_vcs";
//    private static final String INVALID_CREATE_STATS = "invalidCreateStats";
    private static final String VALID_LOAD_STATS = "VariantConfigurationTest_vls";
//    private static final String INVALID_LOAD_STATS = "invalidLoadStats";

    @Autowired
    private Job job;

    @Autowired
    private JobLauncher jobLauncher;

    @Test
    public void validTransform() throws JobExecutionException, IOException {
        String input = VariantConfigurationTest.class.getResource(FILE_20).getFile();
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = VALID_TRANSFORM;

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, "/tmp")
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, ".gz")
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, "studyName")
                .addString(STUDY_ID, "1")
                .addString(FILE_ID, "1")
                .addString(OPENCGA_APP_HOME, opencgaHome)
                .addString(VariantsLoad.SKIP_LOAD, "true")
                .addString(VariantsStatsCreate.SKIP_STATS_CREATE, "true")
                .addString(VariantsStatsLoad.SKIP_STATS_LOAD, "true")
                .toJobParameters();

        String outputFilename = getTransformedOutputPath(Paths.get(FILE_20).getFileName(),
                parameters.getString(STORAGE_COMPRESS_EXTENSION), parameters.getString(STORAGE_OUTPUT_DIR));
        logger.info("transformed output will be at: " + outputFilename);
        File file = new File(outputFilename);
        file.delete();
        assertFalse(file.exists());

        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.COMPLETED.getExitCode(), execution.getExitStatus().getExitCode());

        ////////// check transformed file
        long lines = getLines(new GZIPInputStream(new FileInputStream(outputFilename)));
        assertEquals(300, lines);
    }

    /**
     * This test has to fail because the vcf FILE_WRONG_NO_ALT is malformed, in
     * a variant has an empty alternate allele
     */
    @Test
    public void invalidTransform() throws JobExecutionException {
        String input = VariantConfigurationTest.class.getResource(FILE_WRONG_NO_ALT).getFile();
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = INVALID_TRANSFORM;

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, "/tmp")
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, ".gz")
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, "studyName")
                .addString(STUDY_ID, "2")
                .addString(FILE_ID, "2")
                .addString(OPENCGA_APP_HOME, opencgaHome)
                .addString(VariantsLoad.SKIP_LOAD, "true")
                .addString(VariantsStatsCreate.SKIP_STATS_CREATE, "true")
                .addString(VariantsStatsLoad.SKIP_STATS_LOAD, "true")
                .toJobParameters();

        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.FAILED.getExitCode(), execution.getExitStatus().getExitCode());
    }

    @Test
    public void validLoad() throws JobExecutionException, IllegalAccessException, ClassNotFoundException,
            InstantiationException, IOException, StorageManagerException {
        String input = VariantConfigurationTest.class.getResource(FILE_20).getFile();
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = VALID_LOAD;

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, "/tmp")
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, ".gz")
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, "studyName")
                .addString(STUDY_ID, "1")
                .addString(FILE_ID, "1")
                .addString(OPENCGA_APP_HOME, opencgaHome)
                .addString(VariantsStatsCreate.SKIP_STATS_CREATE, "true")
                .addString(VariantsStatsLoad.SKIP_STATS_LOAD, "true")
                .toJobParameters();

        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.COMPLETED.getExitCode(), execution.getExitStatus().getExitCode());

        // check ((documents in DB) == (lines in transformed file))
        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(dbName, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(new QueryOptions());

        String outputFilename = getTransformedOutputPath(Paths.get(FILE_20).getFileName(),
                parameters.getString(STORAGE_COMPRESS_EXTENSION), parameters.getString(STORAGE_OUTPUT_DIR));
        long lines = getLines(new GZIPInputStream(new FileInputStream(outputFilename)));

        assertEquals(countRows(iterator), lines);
    }

    @Test
    public void validCreateStats() throws JobExecutionException, IOException, InterruptedException,
            IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException {

        String input = VariantConfigurationTest.class.getResource(FILE_20).getFile();
        VariantSource source = new VariantSource(input, "1", "1", "studyName");
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = VALID_CREATE_STATS;
        String compressExtension = ".gz";
        String outputDir = "/tmp";
        File statsFile = new File(Paths.get(outputDir).resolve(VariantStorageManager.buildFilename(source)) + ".variants.stats.json.gz");

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, outputDir)
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, compressExtension)
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, source.getStudyName())
                .addString(STUDY_ID, source.getStudyId())
                .addString(FILE_ID, source.getFileId())
                .addString(OPENCGA_APP_HOME, opencgaHome)
                .addString(VariantsStatsLoad.SKIP_STATS_LOAD, "true")
                .toJobParameters();

        statsFile.delete();
        assertFalse(statsFile.exists());  // ensure the stats file doesn't exist from previous executions
        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.COMPLETED.getExitCode(), execution.getExitStatus().getExitCode());
        assertTrue(statsFile.exists());
    }

    @Test
    public void validLoadStats() throws JobExecutionException, IOException, IllegalAccessException,
            ClassNotFoundException, InstantiationException, StorageManagerException {

        String input = VariantConfigurationTest.class.getResource(FILE_20).getFile();
        VariantSource source = new VariantSource(input, "1", "1", "studyName");
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = VALID_LOAD_STATS;
        String compressExtension = ".gz";
        String outputDir = "/tmp";
        File statsFile = new File(Paths.get(outputDir).resolve(VariantStorageManager.buildFilename(source)) + ".variants.stats.json.gz");

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, outputDir)
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, compressExtension)
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, source.getStudyName())
                .addString(STUDY_ID, source.getStudyId())
                .addString(FILE_ID, source.getFileId())
                .addString(OPENCGA_APP_HOME, opencgaHome)
                .toJobParameters();

        statsFile.delete();
        assertFalse(statsFile.exists());  // ensure the stats file doesn't exist from previous executions
        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.COMPLETED.getExitCode(), execution.getExitStatus().getExitCode());
        assertTrue(statsFile.exists());

        // check ((documents in DB) == (lines in transformed file))
        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(dbName, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(new QueryOptions());

        String outputFilename = getTransformedOutputPath(Paths.get(FILE_20).getFileName(),
                parameters.getString(STORAGE_COMPRESS_EXTENSION), parameters.getString(STORAGE_OUTPUT_DIR));
        long lines = getLines(new GZIPInputStream(new FileInputStream(outputFilename)));

        assertEquals(countRows(iterator), lines);

        // check the DB docs have the field "st"
        variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        variantDBAdaptor = variantStorageManager.getDBAdaptor(dbName, null);
        iterator = variantDBAdaptor.iterator(new QueryOptions());

        assertEquals(1, iterator.next().getSourceEntries().values().iterator().next().getCohortStats().size());
    }

    @BeforeClass
    public static void beforeTests() throws UnknownHostException {
        cleanDBs();
    }

    @AfterClass
    public static void afterTests() throws UnknownHostException {
        cleanDBs();
    }

    private static void cleanDBs() throws UnknownHostException {
        JobTestUtils.cleanDBs(
                VALID_TRANSFORM,
                INVALID_TRANSFORM,
                VALID_LOAD,
                VALID_CREATE_STATS,
                VALID_LOAD_STATS);
    }
}
