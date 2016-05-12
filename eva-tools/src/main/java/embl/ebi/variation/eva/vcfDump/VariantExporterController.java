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

import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;

import javax.ws.rs.core.MultivaluedMap;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.List;

/**
 * Created by jmmut on 2016-01-20.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporterController {

    private final CellBaseClient cellBaseClient;
    private final List<String> studies;
    private final List<String> files;
    private final String outputDir;
    private final VariantSourceDBAdaptor variantSourceDBAdaptor;
    private final VariantDBAdaptor variantDBAdaptor;
    private final QueryOptions query;

    public VariantExporterController(String species, String dbName, List<String> studies, List<String> files,
                                     String outputDir, MultivaluedMap<String, String> queryParameters)
            throws IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException, URISyntaxException {
        this.studies = studies;
        this.files = files;
        this.outputDir = outputDir;
        cellBaseClient = getCellBaseClient(species);
        variantDBAdaptor = getVariantDBAdaptor(dbName);
        query = getQuery(queryParameters);
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

    public QueryOptions getQuery(MultivaluedMap<String, String> queryParameters) {
        QueryOptions query = new QueryOptions();
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);

        for (String acceptedValue : VariantDBAdaptor.QueryParams.acceptedValues) {
            if (queryParameters.containsKey(acceptedValue)) {
                List<String> values = queryParameters.get(acceptedValue);
                String csv = values.get(0);
                for (int i = 1; i < values.size(); i++) {
                    csv += "," + values.get(i);
                }
                query.add(acceptedValue, csv);
            }
        }
        return query;
    }

    public List<String> run() throws URISyntaxException, IOException {
        return new VariantExporter(cellBaseClient, variantDBAdaptor, query).VcfHtsExport(outputDir, variantSourceDBAdaptor);
    }

}
