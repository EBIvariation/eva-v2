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

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

/**
 * Created by pagarcia on 09/05/2016.
 */
public class IntersectingRegionsMerger implements Collector<Region, Set<Region>, List<Region>> {
    @Override
    public Supplier<Set<Region>> supplier() {
        return () -> new TreeSet<>(new RegionComparator());
    }

    @Override
    public BiConsumer<Set<Region>, Region> accumulator() {
        return Set::add;
    }

    @Override
    public BinaryOperator<Set<Region>> combiner() {
        return (left, right) -> {
            left.addAll(right);
            return left;
        };
    }

    @Override
    public Function<Set<Region>, List<Region>> finisher() {
        return this::mergeRegionSetIntoList;
    }

    private List<Region> mergeRegionSetIntoList(Set<Region> regionSet) {
        List<Region> mergedList = new ArrayList<>();
        Region previousRegion = null;
        for (Region region: regionSet) {
            if (intersect(region, previousRegion)) {
                previousRegion = merge(region, previousRegion);
            } else {
                if (previousRegion != null) {
                    mergedList.add(previousRegion);
                }
                previousRegion = region;
            }
        }
        if (previousRegion != null) {
            mergedList.add(previousRegion);
        }
        return mergedList;
    }

    private static boolean intersect(Region region, Region previousRegion) {
        boolean intersect = false;
        if (previousRegion != null) {
            // this method assumes that previousRegion.start <= region.start
            if (region.getChromosome().equals(previousRegion.getChromosome()) && region.getStart() <= previousRegion.getEnd()) {
                intersect = true;
            }
        }
        return intersect;
    }

    private static Region merge(Region region, Region previousRegion) {
        // this method assumes that previousRegion.start <= region.start
        previousRegion.setEnd(Math.max(region.getEnd(), previousRegion.getEnd()));
        return previousRegion;
    }

    @Override
    public Set<Characteristics> characteristics() {
        return EnumSet.of(Characteristics.CONCURRENT);
    }
}
