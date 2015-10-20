/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package embl.ebi.variation.eva;

import org.springframework.batch.core.Job;
import org.springframework.batch.core.Step;
import org.springframework.batch.core.StepContribution;
import org.springframework.batch.core.configuration.annotation.EnableBatchProcessing;
import org.springframework.batch.core.configuration.annotation.JobBuilderFactory;
import org.springframework.batch.core.configuration.annotation.StepBuilderFactory;
import org.springframework.batch.core.job.builder.JobBuilder;
import org.springframework.batch.core.job.builder.SimpleJobBuilder;
import org.springframework.batch.core.launch.support.RunIdIncrementer;
import org.springframework.batch.core.repository.JobRepository;
import org.springframework.batch.core.scope.context.ChunkContext;
import org.springframework.batch.core.step.builder.StepBuilder;
import org.springframework.batch.core.step.builder.TaskletStepBuilder;
import org.springframework.batch.core.step.tasklet.Tasklet;
import org.springframework.batch.repeat.RepeatStatus;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.PropertySource;
import org.springframework.context.support.PropertySourcesPlaceholderConfigurer;
import org.springframework.core.env.Environment;

/**
 *
 * @author Cristina Yenyxe Gonzalez Garcia <cyenyxe@ebi.ac.uk>
 */
@Configuration
@EnableBatchProcessing
@PropertySource("classpath:job-test.properties")
public class SimplestConfig {
    
    public static final String jobName = "simpleJob";

    @Autowired JobRepository jobRepository;
    @Autowired JobBuilderFactory jobBuilderFactory;
    @Autowired StepBuilderFactory stepBuilderFactory;
    
    @Autowired Environment env;
    
    @Autowired SimplestProperties properties;

    @Bean
    public Job simpleJob() {
        JobBuilder jobBuilder = jobBuilderFactory.get(jobName)
                .repository(jobRepository)
                .incrementer(new RunIdIncrementer());

        SimpleJobBuilder simpleJobBuilder = jobBuilder.start(step1());
        return simpleJobBuilder.build();
    }

    Step step1() {
        StepBuilder b = stepBuilderFactory.get("step1");
        TaskletStepBuilder t = b.tasklet(new Tasklet() {

            @Override
            public RepeatStatus execute(StepContribution sc, ChunkContext cc) throws Exception {
                System.out.println(env.getProperty("compressGenotypes"));
                System.out.println(properties.compressGenotypes);
                System.out.println("input = " + properties.input);
                System.out.println("Step 1 done!!");
                return RepeatStatus.FINISHED;
            }
        });
        
        return t.build();
    }
    
    @Bean
    public SimplestProperties simplestProperties(){
        return new SimplestProperties();
    }
    
    // To resolve ${} in @Value
    @Bean
    public static PropertySourcesPlaceholderConfigurer propertyConfigInDev() {
        return new PropertySourcesPlaceholderConfigurer();
    }
    
}
