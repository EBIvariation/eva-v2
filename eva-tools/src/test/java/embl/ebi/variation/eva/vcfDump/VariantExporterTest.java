package embl.ebi.variation.eva.vcfDump;

import org.junit.Test;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;

import java.io.*;
import java.util.Arrays;
import java.util.List;
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
}
