package embl.ebi.variation.eva.vcfDump.regionUtils;

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
