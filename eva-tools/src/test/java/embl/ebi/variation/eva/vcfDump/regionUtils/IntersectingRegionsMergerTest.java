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
package embl.ebi.variation.eva.vcfdump.regionutils;

import org.junit.Test;
import org.opencb.biodata.models.feature.Region;

import java.util.Arrays;
import java.util.Collections;
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
        assertEquals(Arrays.asList(new Region("1", 100, 200), new Region("1", 300, 500), new Region("2", 150, 250)), mergedRegions);
    }

    @Test
    public void testMergeIntersectingRegions() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 100, 200), new Region("1", 150, 500), new Region("2", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(Arrays.asList(new Region("1", 100, 500), new Region("2", 150, 250)), mergedRegions);
    }

    @Test
    public void testMergeUnsortedIntersectingRegions() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 150, 500), new Region("1", 100, 200), new Region("2", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(Arrays.asList(new Region("1", 100, 500), new Region("2", 150, 250)), mergedRegions);
    }

    @Test
    public void testMergeRegionContainedInOther() {
        List<Region> intersectingRegions = Arrays.asList(new Region("1", 100, 300), new Region("1", 150, 250));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(Collections.singletonList(new Region("1", 100, 300)), mergedRegions);
    }

    @Test
    public void testMoreThanTwoIntersectingRegions() {
        List<Region> intersectingRegions =
                Arrays.asList(new Region("1", 100, 200), new Region("1", 150, 500), new Region("1", 400, 600), new Region("1", 125, 550));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(Collections.singletonList(new Region("1", 100, 600)), mergedRegions);
    }

    @Test
    public void testMergeTwoIdenticalRegions() {
        List<Region> intersectingRegions =
                Arrays.asList(new Region("1", 100, 200), new Region("1", 100, 200));
        List<Region> mergedRegions = intersectingRegions.stream().collect(new IntersectingRegionsMerger());
        assertEquals(Collections.singletonList(new Region("1", 100, 200)), mergedRegions);
    }

}