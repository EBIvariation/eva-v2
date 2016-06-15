/*
 * Copyright 2015-2016 EMBL - European Bioinformatics Institute
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
package embl.ebi.variation.eva.vcfdump;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.ws.rs.core.MultivaluedHashMap;
import org.opencb.opencga.lib.common.Config;
import org.slf4j.LoggerFactory;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

/**
 * The variant exporter tool allows to dump a valid VCF from a query against 
 * the EVA database.
 * 
 * Mandatory arguments are: species, database name, studies and files
 * Optional arguments are: output directory
 * 
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 * @author Cristina Yenyxe Gonzalez Garcia &lt;cyenyxe@ebi.ac.uk&gt;
 */
@SpringBootApplication
public class VariantExportBootApplication implements CommandLineRunner {

    private static final org.slf4j.Logger logger = LoggerFactory.getLogger(VariantExportBootApplication.class);
    VariantExportCommand command;
    JCommander commander;

    public VariantExportBootApplication() {
        command = new VariantExportCommand();
        commander = new JCommander(command);
    }

    public void validateArguments(String[] args) throws ParameterException {
        commander.parse(args);
    }

    public void help() {
        commander.usage();
    }

    @Override
    public void run(String[] args) throws Exception {
        try {
            validateArguments(args);
        } catch (ParameterException e) {
            logger.error("Invalid argument: {}", e.getMessage());
            help();
            System.exit(1);
        }

        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");
        
        try {
                    new VariantExporterController(
                            command.species, 
                            command.database, 
                            command.studies, 
                            command.files, 
                            command.outdir,
                            new MultivaluedHashMap<>()).run();
        } catch (Exception e) {
            logger.error("Unsuccessful VCF export: {}", e.getMessage());
            logger.debug("Exception details: ", e);
            System.exit(1);
        }
    }

    public static void main(String[] args) {
        SpringApplication.run(VariantExportBootApplication.class, args);
    }

}
