package embl.ebi.variation.eva.vcfDump;

import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.lib.common.Config;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

/**
 * Created by jmmut on 2015-10-28.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class DumpMain {

    public static void main(String args[]) throws IllegalAccessException, ClassNotFoundException,
            InstantiationException, StorageManagerException, IOException {

        Config.setOpenCGAHome("/opt/opencga");

        QueryOptions query = new QueryOptions();
        QueryOptions options = new QueryOptions();
        List<String> files = Arrays.asList("5");
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

        VariantExporter.VcfHtsExport(iterator, outputStream, variantSourceDBAdaptor, options);


    }
}
