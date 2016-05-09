package embl.ebi.variation.eva.vcfDump.regionUtils;

import org.junit.Test;
import org.opencb.biodata.models.feature.Region;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Created by pagarcia on 09/05/2016.
 */
public class IntersectingRegionsMergerTest {

    @Test
    public void testNoIntersectingRegionList() {
        List<Region> nonIntersectingRegions = Arrays.asList(new Region("1", 100, 200), new Region("1", 300, 500), new Region("2", 150, 250));
        List<Region> mergedRegions = nonIntersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(nonIntersectingRegions, mergedRegions);
    }

    @Test
    public void testSortedOutput() {
        List<Region> unsortedRegions = Arrays.asList(new Region("1", 300, 500), new Region("2", 150, 250), new Region("1", 100, 200));
        List<Region> mergedRegions = unsortedRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 200), new Region("1", 300, 500), new Region("2", 150, 250)));
    }

    @Test
    public void testMergeIntersectingRegions() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 100, 200), new Region("1", 150, 500), new Region("2", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 500), new Region("2", 150, 250)));
    }

    @Test
    public void testMergeUnsortedIntersectingRegions() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 150, 500), new Region("1", 100, 200), new Region("2", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 500), new Region("2", 150, 250)));
    }

    @Test
    public void testMergeRegionContainedInOther() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 100, 300), new Region("1", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 300)));
    }

    @Test
    public void testMoreThanTwoIntersectingRegions() {
        List<Region> intersectingRegions =
                Arrays.asList(new Region("1", 100, 200), new Region("1", 150, 500), new Region("1", 400, 600), new Region("1", 125, 550));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 600)));
    }

    @Test
    public void testMergeTwoIdenticalRegions() {
        List<Region> intersectingRegions =
                Arrays.asList(new Region("1", 100, 200), new Region("1", 100, 200));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(mergedRegions, Arrays.asList(new Region("1", 100, 200)));
    }

}