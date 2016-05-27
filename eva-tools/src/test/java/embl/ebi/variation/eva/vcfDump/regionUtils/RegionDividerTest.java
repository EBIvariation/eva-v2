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

import org.junit.BeforeClass;
import org.junit.Test;
import org.opencb.biodata.models.feature.Region;

import java.util.Arrays;
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

    @Test
    public void testDivideListOfRegions() {
        List<Region> inputRegions = Arrays.asList(new Region("1", 500, 2000), new Region("2", 10000, 12750));
        List<Region> dividedRegions = regionMaker.divideListOfRegions(inputRegions);
        assertTrue(dividedRegions.size() == 5);
        assertTrue(dividedRegions.contains(new Region("1", 500, 1499)));
        assertTrue(dividedRegions.contains(new Region("1", 1500, 2000)));
        assertTrue(dividedRegions.contains(new Region("2", 10000, 10999)));
        assertTrue(dividedRegions.contains(new Region("2", 11000, 11999)));
        assertTrue(dividedRegions.contains(new Region("2", 12000, 12750)));
    }

}