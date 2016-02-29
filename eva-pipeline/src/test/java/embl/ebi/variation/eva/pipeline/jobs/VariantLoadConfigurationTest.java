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

import embl.ebi.variation.eva.pipeline.steps.VariantsStatsCreate;
import embl.ebi.variation.eva.pipeline.steps.VariantsStatsLoad;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
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

import java.io.*;
import java.net.UnknownHostException;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static embl.ebi.variation.eva.pipeline.jobs.JobTestUtils.*;
import static org.junit.Assert.*;
import static embl.ebi.variation.eva.pipeline.listeners.JobParametersListener.*;

/**
 * Created by jmmut on 2015-10-14.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(classes = {VariantLoadConfiguration.class})
public class VariantLoadConfigurationTest {

    public static final String FILE_20 = "/small20.vcf.gz";
    public static final String FILE_22 = "/small22.vcf.gz";
    public static final String FILE_WRONG_NO_ALT = "/wrong_no_alt.vcf.gz";

    private static final Logger logger = LoggerFactory.getLogger(VariantLoadConfigurationTest.class);

    // iterable doing an enum. Does it worth it?
    private static final String VALID_LOAD = "VariantLoadConfigurationTest_v";
    private static final String INVALID_LOAD = "VariantLoadConfigurationTest_i";

    @Autowired
    VariantLoadConfiguration variantConfiguration;

    @Autowired
    private Job job;

    @Autowired
    private JobLauncher jobLauncher;

    @Test
    public void validLoad() throws JobExecutionException, IllegalAccessException, ClassNotFoundException, 
            InstantiationException, StorageManagerException, IOException {
        String input = VariantLoadConfigurationTest.class.getResource(FILE_20).getFile();
        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";
        String dbName = VALID_LOAD;
        
        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, Paths.get(input).getParent().toString())    // reusing transformed path in resources
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, ".gz")
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, "studyName")
                .addString(STUDY_ID, "7")
                .addString(FILE_ID, "10")
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

    /**
     * This test has to fail because the opencgaHome is not set, so it will fail at loading the storage engine configuration.
     */
    @Test
    public void invalidLoad() throws JobExecutionException {
        String input = VariantLoadConfigurationTest.class.getResource(FILE_20).getFile();
        String outdir = input;
        String dbName = INVALID_LOAD;
//        String opencgaHome = System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga";  // TODO make it fail better

        JobParameters parameters = new JobParametersBuilder()
                .addString(INPUT, input)
                .addString(STORAGE_OUTPUT_DIR, outdir)
                .addString(DB_NAME, dbName)
                .addString(STORAGE_COMPRESS_EXTENSION, ".gz")
                .addString(STORAGE_COMPRESS_GENOTYPES, "true")
                .addString(STORAGE_INCLUDE_SRC, "FIRST_8_COLUMNS")
                .addString(STUDY_AGGREGATED, "NONE")
                .addString(STUDY_TYPE, "COLLECTION")
                .addString(STUDY_NAME, "studyName")
                .addString(STUDY_ID, "1")
                .addString(FILE_ID, "1")
                .addString(OPENCGA_APP_HOME, null)
                .addString(VariantsStatsCreate.SKIP_STATS_CREATE, "true")
                .addString(VariantsStatsLoad.SKIP_STATS_LOAD, "true")
                .toJobParameters();

        System.out.println("parameters in load tests" + parameters.toString());
        JobExecution execution = jobLauncher.run(job, parameters);

        assertEquals(input, execution.getJobParameters().getString(INPUT));
        assertEquals(ExitStatus.FAILED.getExitCode(), execution.getExitStatus().getExitCode());
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
        JobTestUtils.cleanDBs(VALID_LOAD, INVALID_LOAD);
    }

}
