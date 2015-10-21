/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package embl.ebi.variation.eva;

import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.runner.RunWith;
import org.springframework.batch.core.*;
import org.springframework.batch.core.launch.JobLauncher;
import org.springframework.batch.test.JobLauncherTestUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.SpringApplication;
import org.springframework.context.annotation.Bean;
import org.springframework.core.env.Environment;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Cristina Yenyxe Gonzalez Garcia <cyenyxe@ebi.ac.uk>
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration(classes = SimplestConfig.class)
//@ContextConfiguration(classes = VariantConfiguration.class)
public class SimplestTest {

    @Autowired Job job;
    @Autowired JobLauncher jobLauncher;
    @Autowired Environment env;
    @Autowired SimplestProperties properties;
    @Autowired SimplestConfig simplestConfig;

    @Test public void simpleTest() throws Exception {
        Map<String, JobParameter> params = new HashMap<>();
        String testInput = "testInput";
        params.put("--input", new JobParameter(testInput));

        System.out.println("en simple test");
        properties.input = testInput;
        JobExecution execution = jobLauncher.run(job, new JobParameters());

        assertTrue(env.getProperty("compressGenotypes", Boolean.class));
        assertTrue(properties.compressGenotypes);
        assertEquals(testInput, simplestConfig.properties.input);
        assertEquals(ExitStatus.COMPLETED, execution.getExitStatus());
    }

//    @Test
//    public void simpleMainTest() {
//
//        String testInput = "testInput";
//        String[] args = {
//                "--spring.batch.job.names=" + SimplestConfig.jobName,
//                "--input=" + testInput,
//        };
//
//        SpringApplication.run(Application.class, args);
//
//
//
//    }
    
    
//        @Autowired
//        JobLauncherTestUtils jobLauncherTestUtils;
//
//	@Test
//        public void simpleTest() throws Exception {
//            jobLauncherTestUtils.launchJob();
//            JobExecution execution = jobLauncherTestUtils.getJobLauncher().run(jobLauncherTestUtils.getJob(), new JobParameters());
//            assertEquals(ExitStatus.COMPLETED, execution.getExitStatus());
//	}
//	
//        @Bean
//        public JobLauncherTestUtils jobLauncherTestUtils() {
//            return new JobLauncherTestUtils();
//        }
}