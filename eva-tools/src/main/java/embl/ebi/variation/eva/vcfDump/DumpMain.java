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
package embl.ebi.variation.eva.vcfDump;

import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;

import javax.ws.rs.core.MultivaluedHashMap;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.List;

/**
 * Created by jmmut on 2015-10-28.
 *
 * Example skeleton usage of the VariantExporter. An actual CLI is intended to be developed here.
 * This class doesn't work out of the box, as it needs an indexed vcf file in mongo with some hardcoded values:
 * - fileId: "5"
 * - studyId: "7"
 * - dbName: "batch";
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class DumpMain {

    public static void main(String args[]) throws IllegalAccessException, ClassNotFoundException,
            InstantiationException, StorageManagerException, IOException, URISyntaxException {

        Config.setOpenCGAHome(System.getenv("OPENCGA_HOME") != null ? System.getenv("OPENCGA_HOME") : "/opt/opencga");

        List<String> files;         // = Arrays.asList("5");
        List<String> studies;       // = Arrays.asList("7");
        String dbName;              // = "batch";
        String outputDir;           // = "./";
        String species;

        if (args.length == 5) {
            species = args[0];
            dbName = args[1];
            studies = Arrays.asList(args[2].split(","));
            files = Arrays.asList(args[3].split(","));
            outputDir = args[4];
        } else {
            System.out.println("usage: java -jar <jar> <species> <dbName> <studies CommaSeparatedValues> <files CSV> <output directory>");
            System.out.println("example: java -jar eva-tools-0.1.jar hsapiens batch 7 5,6 ./");
            return;
        }

        List<String> fileNames = new VariantExporterController(
                species, dbName, studies, files, outputDir, new MultivaluedHashMap<String, String>()).run();
    }
}
