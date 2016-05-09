package embl.ebi.variation.eva.vcfDump.regionUtils;

import embl.ebi.variation.eva.vcfDump.regionUtils.RegionDivider;
import org.junit.BeforeClass;
import org.junit.Test;
import org.opencb.biodata.models.feature.Region;

import java.util.List;

import static org.junit.Assert.*;

/**
 * Created by pagarcia on 03/05/2016.
 */
public class RegionDividerTest {

    private static RegionDivider regionMaker;

    @BeforeClass
    public static void setUpClass() throws Exception {
        regionMaker = new RegionDivider(1000);
    }

    @Test
    public void testGetRegions() {
        List<Region> regions = regionMaker.divideRegionInChunks("1", 500, 1499);
        assertTrue(regions.size() == 1);
        assertTrue(regions.contains(new Region("1", 500, 1499)));

        regions = regionMaker.divideRegionInChunks("1", 500, 1000);
        assertTrue(regions.size() == 1);
        assertTrue(regions.contains(new Region("1", 500, 1000)));

        regions = regionMaker.divideRegionInChunks("1", 2500, 2500);
        assertTrue(regions.size() == 1);
        assertTrue(regions.contains(new Region("1", 2500, 2500)));

        regions = regionMaker.divideRegionInChunks("1", 500, 2499);
        assertTrue(regions.size() == 2);
        assertTrue(regions.contains(new Region("1", 500, 1499)));
        assertTrue(regions.contains(new Region("1", 1500, 2499)));

        regions = regionMaker.divideRegionInChunks("1", 500, 2000);
        assertTrue(regions.size() == 2);
        assertTrue(regions.contains(new Region("1", 500, 1499)));
        assertTrue(regions.contains(new Region("1", 1500, 2000)));

        regions = regionMaker.divideRegionInChunks("1", 500, 2500);
        assertTrue(regions.size() == 3);
        assertTrue(regions.contains(new Region("1", 500, 1499)));
        assertTrue(regions.contains(new Region("1", 1500, 2499)));
        assertTrue(regions.contains(new Region("1", 2500, 2500)));

        regions = regionMaker.divideRegionInChunks("1", -1, 2500);
        assertTrue(regions.size() == 0);

        regions = regionMaker.divideRegionInChunks("1", 3000, 2500);
        assertTrue(regions.size() == 0);
    }


}