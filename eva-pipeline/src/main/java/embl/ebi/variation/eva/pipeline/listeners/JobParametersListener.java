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
package embl.ebi.variation.eva.pipeline.listeners;

import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantStudy;
import org.opencb.datastore.core.ObjectMap;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobExecutionListener;
import org.springframework.batch.core.JobParameters;
import org.springframework.stereotype.Component;


public abstract class JobParametersListener implements JobExecutionListener {
     
    private static final Logger logger = LoggerFactory.getLogger(JobParametersListener.class);

    public static final String DB_NAME = "db.name";
    public static final String COMPRESS_GENOTYPES = "compress.genotypes";
    public static final String INPUT = "input";
    public static final String PEDIGREE = "pedigree";
    public static final String OUTPUT_DIR = "output.dir";
    public static final String FILE_ID = "file.id";
    public static final String STUDY_ID = "study.id";
    public static final String STUDY_NAME = "study.name";
    public static final String STUDY_TYPE = "study.type";
    public static final String STUDY_AGGREGATED = "study.aggregated";
    public static final String INCLUDE_SRC = "include.src";
    public static final String INCLUDE_STATS = "include.stats";
    public static final String COMPRESS_EXTENSION = "compress.extension";
    public static final String OPENCGA_APP_HOME = "opencga.app.home";

    protected final ObjectMap variantOptions;
    
    public JobParametersListener() {
        variantOptions = new ObjectMap();
    }
    
    @Override
    public void afterJob(JobExecution jobExecution) {
        logger.info("afterJob STATUS + " + jobExecution.getStatus());
        logger.info("afterJob : " + jobExecution);
    }
    
    public ObjectMap getVariantOptions() {
        return variantOptions;
    }
}
