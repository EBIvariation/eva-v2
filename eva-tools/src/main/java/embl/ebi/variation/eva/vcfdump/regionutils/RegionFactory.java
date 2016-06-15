/*
 *
 *  * Copyright 2016 EMBL - European Bioinformatics Institute
 *  *
 *  * Licensed under the Apache License, Version 2.0 (the "License");
 *  * you may not use this file except in compliance with the License.
 *  * You may obtain a copy of the License at
 *  *
 *  *      http://www.apache.org/licenses/LICENSE-2.0
 *  *
 *  * Unless required by applicable law or agreed to in writing, software
 *  * distributed under the License is distributed on an "AS IS" BASIS,
 *  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  * See the License for the specific language governing permissions and
 *  * limitations under the License.
 *
 */

package embl.ebi.variation.eva.vcfdump.regionutils;

import com.mongodb.BasicDBObject;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by pagarcia on 27/05/2016.
 */
public class RegionFactory {

    private static final Logger logger = LoggerFactory.getLogger(RegionFactory.class);

    private int windowSize;
    private final VariantDBAdaptor variantAdaptor;
    private final QueryOptions query;
    private List<Region> regionsInFilter;

    public RegionFactory(int windowSize, VariantDBAdaptor variantAdaptor, QueryOptions query) {
        this.windowSize = windowSize;
        this.variantAdaptor = variantAdaptor;
        this.query = query;
    }

    public List<Region> getRegionsForChromosome(String chromosome) {
        String regionFilter = query.getString(VariantDBAdaptor.REGION);
        if (regionFilter == null || regionFilter.isEmpty()) {
            int minStart = getMinStart(chromosome);
            if (minStart == -1) {
                return Collections.EMPTY_LIST;
            } else {
                int maxStart = getMaxStart(chromosome);
                logger.debug("Chromosome {} maxStart: {}", chromosome, maxStart);
                logger.debug("Chromosome {} minStart: {}", chromosome, minStart);
                return divideChromosomeInChunks(chromosome, minStart, maxStart);
            }
        } else {
            List<Region> chromosomeRegionsFromQuery =
                    getRegionsFromQuery(regionFilter).stream().filter(r -> r.getChromosome().equals(chromosome)).collect(new IntersectingRegionsMerger());

            String commaSeparatedRegionList = chromosomeRegionsFromQuery.stream().map(Region::toString).collect(Collectors.joining(", "));
            logger.debug("Chromosome {} regions from query: {}", chromosome, commaSeparatedRegionList);

            return divideRegionListInChunks(chromosomeRegionsFromQuery);

        }
    }

    private List<Region> divideChromosomeInChunks(String chromosome, int chromosomeMinStart, int chromosomeMaxStart) {
        List<Region> regions = divideRegionInChunks(chromosome, chromosomeMinStart, chromosomeMaxStart);
        logger.debug("Number of regions in chromosome{}: {}", chromosome, regions.size());
        if (!regions.isEmpty()) {
            logger.debug("First region: {}", regions.get(0));
            logger.debug("Last region: {}", regions.get(regions.size() - 1));
        }

        return regions;
    }

    private int getMinStart(String chromosome) {
        QueryOptions minQuery = addChromosomeSortAndLimitToQuery(chromosome, true);
        return getVariantStart(minQuery);
    }

    private int getMaxStart(String chromosome) {
        QueryOptions maxQuery = addChromosomeSortAndLimitToQuery(chromosome, false);
        return getVariantStart(maxQuery);
    }

    private QueryOptions addChromosomeSortAndLimitToQuery(String chromosome, boolean ascending) {
        QueryOptions chromosomeSortedByStartQuery = new QueryOptions(query);
        chromosomeSortedByStartQuery.put(VariantDBAdaptor.CHROMOSOME, chromosome);

        BasicDBObject sortDBObject = new BasicDBObject();
        int orderOperator = ascending ? 1 : -1;
        sortDBObject.put("chr", orderOperator);
        sortDBObject.put("start", orderOperator);
        chromosomeSortedByStartQuery.put("sort", sortDBObject);

        chromosomeSortedByStartQuery.put("limit", 1);

        return chromosomeSortedByStartQuery;
    }

    private int getVariantStart(QueryOptions query) {
        int start = -1;
        VariantDBIterator variantDBIterator = variantAdaptor.iterator(query);
        if (variantDBIterator.hasNext()) {
            Variant variant = variantDBIterator.next();
            start = variant.getStart();
        }

        return start;
    }

    private List<Region> getRegionsFromQuery(String regionFilter) {
        if (regionsInFilter == null) {
            regionsInFilter = Region.parseRegions(regionFilter);
        }
        return regionsInFilter;
    }

    private List<Region> divideRegionListInChunks(List<Region> regionsFromQuery) {
        List<Region> regions = new ArrayList<>();
        for (Region region : regionsFromQuery) {
            regions.addAll(divideRegionInChunks(region.getChromosome(), region.getStart(), region.getEnd()));
        }
        return regions;
    }

    public List<Region> divideRegionInChunks(String chromosome, int minStart, int maxStart) {
        if (minStart == -1) {
            return Collections.EMPTY_LIST;
        } else {
            List<Region> regions = new ArrayList<>();
            int nextStart = minStart;
            while (nextStart <= maxStart) {
                int end = Math.min(nextStart + windowSize, maxStart + 1);
                regions.add(new Region(chromosome, nextStart, end - 1));
                nextStart = end;
            }
            return regions;
        }
    }

}
