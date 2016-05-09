package embl.ebi.variation.eva.vcfDump.regionUtils;

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
        mergedList.add(previousRegion);
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
