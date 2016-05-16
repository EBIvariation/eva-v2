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

import java.util.Comparator;

/**
 * Created by pagarcia on 06/05/2016.
 */
public class RegionComparator implements Comparator<Region> {
    @Override
    public int compare(Region r1, Region r2) {
        if (!r1.getChromosome().equals(r2.getChromosome())) {
            return r1.getChromosome().compareTo(r2.getChromosome());
        } else {
            if (r1.getStart() != r2.getStart()) {
                return r1.getStart() - r2.getStart();
            } else {
                return r1.getEnd() - r2.getEnd();
            }
        }
    }
}
