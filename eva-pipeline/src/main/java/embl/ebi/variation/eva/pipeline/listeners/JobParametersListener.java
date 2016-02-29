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

import org.opencb.datastore.core.ObjectMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobExecutionListener;


public abstract class JobParametersListener implements JobExecutionListener {
     
    private static final Logger logger = LoggerFactory.getLogger(JobParametersListener.class);

    public static final String INPUT = "input";
    public static final String DB_NAME = "db.name";
    public static final String STORAGE_COMPRESS_GENOTYPES = "storage.compress.genotypes";
    public static final String STORAGE_COMPRESS_EXTENSION = "storage.compress.extension";
    public static final String STORAGE_INCLUDE_SRC = "storage.include.src";
    public static final String STORAGE_OUTPUT_DIR = "storage.output.dir";
    public static final String PEDIGREE = "pedigree";   // TODO remove ?
    public static final String STATS_OUTPUT_DIR = "stats.output.dir";   // TODO use this in tests and stats steps
    public static final String FILE_ID = "file.id";
    public static final String STUDY_ID = "study.id";
    public static final String STUDY_NAME = "study.name";
    public static final String STUDY_TYPE = "study.type";
    public static final String STUDY_AGGREGATED = "study.aggregated";
    public static final String STUDY_AGGREGATED_MAPPING_FILE = "study.aggregated.mapping.file";
    public static final String STATS_INCLUDE = "stats.include";
    public static final String STATS_OVERWRITE = "stats.overwrite"; // TODO when merged feature/overwrite-stats don't use overwriteStats, use this one
    public static final String STATS_COHORTS = "stats.cohorts";
    public static final String OPENCGA_APP_HOME = "opencga.app.home";
    // TODO when merging feature/annotation, don't forget to add "annot.overwrite" that imitates to opencga's overwriteAnnotations

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
