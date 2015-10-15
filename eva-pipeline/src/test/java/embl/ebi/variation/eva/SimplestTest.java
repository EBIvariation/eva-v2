/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package embl.ebi.variation.eva;

import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.runner.RunWith;
import org.springframework.batch.core.ExitStatus;
import org.springframework.batch.core.Job;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobParameters;
import org.springframework.batch.core.launch.JobLauncher;
import org.springframework.batch.test.JobLauncherTestUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Bean;
import org.springframework.core.env.Environment;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

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
	
	@Test public void simpleTest() throws Exception {
            JobExecution execution = jobLauncher.run(job, new JobParameters());
            assertTrue(env.getProperty("compressGenotypes", Boolean.class));
            assertEquals(ExitStatus.COMPLETED, execution.getExitStatus());
	}
	
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