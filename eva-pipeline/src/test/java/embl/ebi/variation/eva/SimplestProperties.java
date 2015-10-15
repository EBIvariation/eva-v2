/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package embl.ebi.variation.eva;

import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantStudy;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.PropertySource;
import org.springframework.context.support.PropertySourcesPlaceholderConfigurer;

/**
 *
 * @author Cristina Yenyxe Gonzalez Garcia <cyenyxe@ebi.ac.uk>
 */
public class SimplestProperties {
    
    @Value("${input}")              public String input;
    @Value("${outputDir}")          public String outputDir;
    @Value("${pedigree}")           public String pedigree;
    @Value("${dbName}")             public String dbName;
    @Value("${storageEngine}")      public String storageEngine;
    @Value("${compressGenotypes}")  public boolean compressGenotypes;
    @Value("${compressExtension}")  public String compressExtension;
    @Value("${includeSrc}")         public String includeSrc;
//    @Value("${credentials}")        public String credentials;
//    @Value("${opencga.configfile}") public String configFile;
    @Value("${opencga.app.home}")   public String appHome;
    @Value("${fileId}")             public String fileId;
    @Value("${aggregated}")         public VariantSource.Aggregation aggregated;
    @Value("${studyType}")          public VariantStudy.StudyType studyType;
    @Value("${studyName}")          public String studyName;
    @Value("${studyId}")            public String studyId;

    // transform

    // load
    @Value("${loadThreads}")        public Integer loadThreads;
    // Integer bulkSize, batchSize?

    //stats
    @Value("${calculateStats}")     public boolean calculateStats;
    @Value("${overwriteStats}")     public boolean overwriteStats;

    // annotation
    @Value("${annotate}")           public boolean annotate;

    // job repository DB
    @Value("${jobRepositoryDriverClassName}")   public String jobRepositoryDriverClassName;
    @Value("${jobRepositoryUrl}")               public String jobRepositoryUrl;
    @Value("${jobRepositoryUsername}")          public String jobRepositoryUsername;
    @Value("${jobRepositoryPassword}")          public String jobRepositoryPassword;
    
}
