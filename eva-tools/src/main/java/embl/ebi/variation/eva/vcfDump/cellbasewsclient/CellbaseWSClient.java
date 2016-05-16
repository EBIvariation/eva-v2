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
package embl.ebi.variation.eva.vcfdump.cellbasewsclient;

import org.opencb.biodata.models.feature.Region;
import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.cellbase.core.common.GenomeSequenceFeature;
import org.opencb.datastore.core.QueryResponse;
import org.opencb.datastore.core.QueryResult;
import org.opencb.opencga.lib.common.Config;
import org.springframework.core.ParameterizedTypeReference;
import org.springframework.http.HttpMethod;
import org.springframework.http.ResponseEntity;
import org.springframework.web.client.RestTemplate;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Created by pagarcia on 12/05/2016.
 * This class encapsulates CellBaseClient adding some operations to its API
 */
public class CellbaseWSClient {
    private final String species;
    private final CellBaseClient cellbaseClient;
    private final String cellbaseRestURL;
    private final String cellbaseRestVersion;


    public CellbaseWSClient (String species) throws URISyntaxException {
        this.species = species;
        this.cellbaseRestURL = (String) Config.getStorageProperties().get("CELLBASE.REST.URL");
        this.cellbaseRestVersion = (String) Config.getStorageProperties().get("CELLBASE.VERSION");
        this.cellbaseClient = getClient(species);
    }

    private CellBaseClient getClient(String species) throws URISyntaxException {
        CellBaseClient cellBaseClient = new CellBaseClient(new URI(cellbaseRestURL), cellbaseRestVersion, species);
        return cellBaseClient;
    }

    public String getSequence(Region region) throws IOException {
        String regionSequence = null;
        List<Region> regions = Collections.singletonList(region);
        QueryResponse<QueryResult<GenomeSequenceFeature>> sequence = cellbaseClient.getSequence(
                CellBaseClient.Category.genomic, CellBaseClient.SubCategory.region, regions, null);

        List<GenomeSequenceFeature> response = sequence.getResponse().get(0).getResult();
        if (response.size() == 1) {
            regionSequence = response.get(0).getSequence();
        }
        return regionSequence;
    }

    public Set<String> getChromosomes() {
        // call cellbase chromosomes WS
        RestTemplate restTemplate = new RestTemplate();
        ParameterizedTypeReference<QueryResponse<QueryResult<CellbaseChromosomesWSOutput>>> responseType =
                new ParameterizedTypeReference<QueryResponse<QueryResult<CellbaseChromosomesWSOutput>>>() {};
        ResponseEntity<QueryResponse<QueryResult<CellbaseChromosomesWSOutput>>> wsOutput =
                restTemplate.exchange(cellbaseRestURL + "/" + cellbaseRestVersion + "/" + species + "/genomic/chromosome/all",
                        HttpMethod.GET, null, responseType);

        // parse WS output and return all chromosome names
        QueryResponse<QueryResult<CellbaseChromosomesWSOutput>> response = wsOutput.getBody();
        QueryResult<CellbaseChromosomesWSOutput> result = response.getResponse().get(0);
        CellbaseChromosomesWSOutput results = result.getResult().get(0);
        return results.getAllChromosomeNames();
    }

    public String getUrl() {
        return cellbaseRestURL;
    }

    public String getVersion() {
        return cellbaseRestVersion;
    }
}
