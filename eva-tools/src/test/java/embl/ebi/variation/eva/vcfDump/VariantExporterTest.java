package embl.ebi.variation.eva.vcfDump;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.Test;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantFactory;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantVcfFactory;
import org.opencb.biodata.models.variant.stats.VariantSourceStats;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.datastore.core.QueryResult;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static org.junit.Assert.*;

/**
 * Created by jmmut on 2015-10-29.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporterTest {

    @Test
    public void testVcfHtsExport() throws Exception {

        Config.setOpenCGAHome("/opt/opencga");

        QueryOptions query = new QueryOptions();
        QueryOptions options = new QueryOptions();
//        List<String> files = Arrays.asList("5");
        List<String> files = Arrays.asList("5", "6");
        List<String> studies = Arrays.asList("7");
        String dbName = "batch";
        String fileName = "exported.vcf.gz";


        OutputStream outputStream = new GZIPOutputStream(new FileOutputStream(fileName));
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);


        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(dbName, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();

        int failedVariants = VariantExporter.VcfHtsExport(iterator, outputStream, variantSourceDBAdaptor, options);

        ////////// checks

        // test file should not have failed variants
        assertEquals(0, failedVariants);

        // counting lines (without comments)
        BufferedReader file = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
        long lines = 0;
        String line;
        while ((line = file.readLine()) != null) {
            if (line.charAt(0) != '#') {
                lines++;
            }
        }
        file.close();

        // counting variants in the DB
        iterator = variantDBAdaptor.iterator(query);
        int variantRows = 0;
        while(iterator.hasNext()) {
            iterator.next();
            variantRows++;
        }

        assertEquals(variantRows, lines);
    }

    @Test
    public void testGetVariantContextFromVariant() throws Exception {
        final VariantSource variantSource = new VariantSource("name", "fileId", "studyId", "studyName");
        List<String> samples = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            samples.add("s"+i);
        }
        variantSource.setSamples(samples);
        VariantFactory factory = new VariantVcfFactory();
        VariantSourceDBAdaptor sourceDBAdaptor = new MockVariantSourceDBAdaptor(variantSource);
        QueryOptions options = new QueryOptions("fileId", "fileId");
        options.add("studyId", "studyId");
        List<Variant> variants;
        VariantContext variantContext;
        List<String> alleles;

        // test multiallelic
        String multiallelicLine = "1\t1000\tid\tC\tA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicLine);
        assertEquals(2, variants.size());

        variantContext = VariantExporter.convertBiodataVariantToVariantContext(variants.get(0), sourceDBAdaptor, options);
        alleles = Arrays.asList("C", "A", ".");
        assertEqualGenotypes(variants.get(0), variantContext, alleles);

        variantContext = VariantExporter.convertBiodataVariantToVariantContext(variants.get(1), sourceDBAdaptor, options);
        alleles = Arrays.asList("C", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext, alleles);


        // test indel
        String indelLine = "1\t1000\tid\tC\tCA\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|0\t0|1";
        variants = factory.create(variantSource, indelLine);

        variantContext = VariantExporter.convertBiodataVariantToVariantContext(variants.get(0), sourceDBAdaptor, options);
        alleles = Arrays.asList("C", "CA");
        assertEqualGenotypes(variants.get(0), variantContext, alleles);


        // test multiallelic + indel
        String multiallelicIndelLine = "1\t1000\tid\tC\tCA,T\t100\tPASS\t.\tGT\t0|0\t0|0\t0|1\t1|1\t1|2\t0|1";
        variants = factory.create(variantSource, multiallelicIndelLine);
        assertEquals(2, variants.size());

        variantContext = VariantExporter.convertBiodataVariantToVariantContext(variants.get(0), sourceDBAdaptor, options);
        alleles = Arrays.asList("C", "CA", ".");
        assertEqualGenotypes(variants.get(0), variantContext, alleles);

        variantContext = VariantExporter.convertBiodataVariantToVariantContext(variants.get(1), sourceDBAdaptor, options);
        alleles = Arrays.asList("C", "T", ".");
        assertEqualGenotypes(variants.get(1), variantContext, alleles);

    }

    private void assertEqualGenotypes(Variant variant, VariantContext variantContext, List<String> alleles) {
        for (Map.Entry<String, Map<String, String>> data : variant.getSourceEntries().values().iterator().next().getSamplesData().entrySet()) {
            Genotype genotype = variantContext.getGenotype(data.getKey());
            String gt = data.getValue().get("GT");
            org.opencb.biodata.models.feature.Genotype biodataGenotype = new org.opencb.biodata.models.feature.Genotype(gt, alleles.get(0), alleles.get(1));
            assertEquals(Allele.create(alleles.get(biodataGenotype.getAllele(0)), biodataGenotype.isAlleleRef(0)),
                    genotype.getAllele(0));
            assertEquals(Allele.create(alleles.get(biodataGenotype.getAllele(1)), biodataGenotype.isAlleleRef(1)),
                    genotype.getAllele(1));
        }
    }

    public class MockVariantSourceDBAdaptor implements VariantSourceDBAdaptor {
        private VariantSource variantSource;

        public MockVariantSourceDBAdaptor(VariantSource variantSource) {
            this.variantSource = variantSource;
        }

        @Override
        public QueryResult countSources() {
            return null;
        }

        @Override
        public QueryResult getAllSources(QueryOptions options) {
            QueryResult<VariantSource> queryResult = new QueryResult<>();
            queryResult.setResult(Collections.singletonList(variantSource));
            return queryResult;
        }

        @Override
        public QueryResult getAllSourcesByStudyId(String studyId, QueryOptions options) {
            return null;
        }

        @Override
        public QueryResult getAllSourcesByStudyIds(List<String> studyIds, QueryOptions options) {
            return null;
        }

        @Override
        public QueryResult getSamplesBySource(String fileId, QueryOptions options) {
            return null;
        }

        @Override
        public QueryResult getSamplesBySources(List<String> fileIds, QueryOptions options) {
            return null;
        }

        @Override
        public QueryResult getSourceDownloadUrlByName(String filename) {
            return null;
        }

        @Override
        public List<QueryResult> getSourceDownloadUrlByName(List<String> filenames) {
            return null;
        }

        @Override
        public QueryResult getSourceDownloadUrlById(String fileId, String studyId) {
            return null;
        }

        @Override
        public QueryResult updateSourceStats(VariantSourceStats variantSourceStats, QueryOptions queryOptions) {
            return null;
        }

        @Override
        public boolean close() {
            return false;
        }
    }
}
