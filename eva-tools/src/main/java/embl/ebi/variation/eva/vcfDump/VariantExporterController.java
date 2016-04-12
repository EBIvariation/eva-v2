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
package embl.ebi.variation.eva.vcfDump;

import org.opencb.biodata.models.variant.Variant;
import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Iterator;
import java.util.List;

/**
 * Created by jmmut on 2016-01-20.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporterController {

    private final CellBaseClient cellBaseClient;
    private final String outputDir;
    private final VariantSourceDBAdaptor variantSourceDBAdaptor;
    private final VariantDBAdaptor variantDBAdaptor;

    public VariantExporterController(String species, String dbName, List<String> studies, List<String> files, String outputDir)
            throws IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException, URISyntaxException {
        this.outputDir = outputDir;
        cellBaseClient = getCellBaseClient(species);
        variantDBAdaptor = getVariantDBAdaptor(dbName);
        variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
    }

    public CellBaseClient getCellBaseClient(String species) throws URISyntaxException {
        String url = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        String version = (String) Config.getStorageProperties().get("CELLBASE.VERSION");
        CellBaseClient cellBaseClient = new CellBaseClient(new URI(url), version, species);
        return cellBaseClient;
    }

    public VariantDBAdaptor getVariantDBAdaptor(String dbName) throws StorageManagerException, IllegalAccessException, ClassNotFoundException, InstantiationException {
        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        return variantStorageManager.getDBAdaptor(dbName, null);
    }

    public List<String> run(List<String> studies, List<String> files) throws URISyntaxException, IOException {
        return new VariantExporter(cellBaseClient, variantDBAdaptor, files, studies).VcfHtsExport(outputDir, variantSourceDBAdaptor);
    }

}
