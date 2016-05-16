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
package embl.ebi.variation.eva.vcfdump.regionutils;

import org.opencb.biodata.models.feature.Region;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by pagarcia on 03/05/2016.
 */
public class RegionDivider {

    private int windowSize;

    public RegionDivider(int windowSize) {
        this.windowSize = windowSize;
    }

    public List<Region> divideRegionInChunks(Region region) {
        return divideRegionInChunks(region.getChromosome(), region.getStart(), region.getEnd());
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
