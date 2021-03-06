/*
 * Copyright 2016 EMBL - European Bioinformatics Institute
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

package embl.ebi.variation.eva.pipeline.steps;

import embl.ebi.variation.eva.pipeline.MongoDBHelper;
import embl.ebi.variation.eva.pipeline.annotation.GzipLazyResource;
import embl.ebi.variation.eva.pipeline.annotation.load.VariantAnnotationLineMapper;
import embl.ebi.variation.eva.pipeline.annotation.load.VariantAnnotationMongoItemWriter;
import embl.ebi.variation.eva.pipeline.jobs.VariantJobArgsConfig;
import embl.ebi.variation.eva.pipeline.listener.SkipCheckingListener;
import org.opencb.biodata.models.variant.annotation.VariantAnnotation;
import org.opencb.datastore.core.ObjectMap;
import org.springframework.batch.core.Step;
import org.springframework.batch.core.configuration.annotation.EnableBatchProcessing;
import org.springframework.batch.core.configuration.annotation.StepBuilderFactory;
import org.springframework.batch.item.ItemWriter;
import org.springframework.batch.item.data.MongoItemWriter;
import org.springframework.batch.item.file.FlatFileItemReader;
import org.springframework.batch.item.file.FlatFileParseException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Import;
import org.springframework.core.io.Resource;
import org.springframework.data.mongodb.core.MongoOperations;

import java.io.IOException;

/**
 * @author Diego Poggioli
 *
 * Step class that:
 * - READ: read a list of VEP {@link VariantAnnotation} from flat file
 * - LOAD: write the {@link VariantAnnotation} into Mongo db
 *
 */

@Configuration
@EnableBatchProcessing
@Import(VariantJobArgsConfig.class)
public class VariantsAnnotLoad {

    @Autowired
    private StepBuilderFactory steps;

    @Autowired
    private ObjectMap pipelineOptions;

    @Bean
    @Qualifier("variantAnnotLoadBatchStep")
    public Step variantAnnotLoadBatchStep() throws IOException {
        return steps.get("variantAnnotLoadBatchStep").<VariantAnnotation, VariantAnnotation> chunk(10)
                .reader(variantAnnotationReader())
                .writer(variantAnnotationWriter())
                .faultTolerant().skipLimit(50).skip(FlatFileParseException.class)
                .listener(skipCheckingListener())
                .build();
    }

    @Bean
    public FlatFileItemReader<VariantAnnotation> variantAnnotationReader() throws IOException {
        Resource resource = new GzipLazyResource(pipelineOptions.getString("vepOutput"));
        FlatFileItemReader<VariantAnnotation> reader = new FlatFileItemReader<>();
        reader.setResource(resource);
        reader.setLineMapper(new VariantAnnotationLineMapper());
        return reader;
    }

    @Bean
    public ItemWriter<VariantAnnotation> variantAnnotationWriter(){
        MongoOperations mongoOperations = MongoDBHelper.getMongoOperationsFromPipelineOptions(pipelineOptions);
        MongoItemWriter<VariantAnnotation> writer = new VariantAnnotationMongoItemWriter(mongoOperations);
        writer.setCollection(pipelineOptions.getString("dbCollectionVariantsName"));
        writer.setTemplate(mongoOperations);
        return writer;
    }

    @Bean
    public SkipCheckingListener skipCheckingListener(){
        return new SkipCheckingListener();
    }

}
