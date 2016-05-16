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
package embl.ebi.variation.eva.vcfdump;

import com.beust.jcommander.IValueValidator;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import java.nio.file.Paths;
import java.util.List;

/**
 *
 * @author Cristina Yenyxe Gonzalez Garcia &lt;cyenyxe@ebi.ac.uk&gt;
 */
public class VariantExportCommand {

    @Parameter(names = "--species", required = true, description = "Species the data are associated to")
    String species;

    @Parameter(names = "--database", required = true, description = "Name of the database to extract data from")
    String database;
    
    @Parameter(names = "--outdir", description = "Output directory", validateValueWith = PathValidator.class)
    String outdir = System.getProperty("user.dir");
    
    @Parameter(names = "--studies", required = true, description = "Comma-separated list of studies to query")
    List<String> studies;

    @Parameter(names = "--files", required = true, description = "Comma-separated list of files to query")
    List<String> files;
    
    
    public static class PathValidator implements IValueValidator {
        
        @Override
        public void validate(String name, Object value) throws ParameterException {
            try {
                Paths.get(value.toString());
            } catch (Exception e) {
                throw new ParameterException(e);
            }
        }
    }
}
