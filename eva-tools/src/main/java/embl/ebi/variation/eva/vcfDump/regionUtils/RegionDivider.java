package embl.ebi.variation.eva.vcfDump.regionUtils;

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
