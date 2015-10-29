package embl.ebi.variation.eva.vcfDump;

import org.junit.Test;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;

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

        VariantExporter variantExporter = new VariantExporter();


        QueryOptions query = new QueryOptions();
        QueryOptions options = new QueryOptions();
        List<String> files = Arrays.asList("5");
        List<String> studies = Arrays.asList("7");
        String dbName = "batch";
        String fileName = "exported.vcf";


        OutputStream outputStream = new FileOutputStream(fileName);
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);


        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        VariantDBAdaptor variantDBAdaptor = variantStorageManager.getDBAdaptor(dbName, null);
        VariantDBIterator iterator = variantDBAdaptor.iterator(query);
        VariantSourceDBAdaptor variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();

        variantExporter.VcfHtsExport(iterator, outputStream, variantSourceDBAdaptor, options);
    }
}
