package embl.ebi.variation.eva.vcfDump;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantSourceEntry;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.datastore.core.QueryResult;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.opencb.opencga.storage.mongodb.variant.DBObjectToVariantSourceConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.stream.Stream;

/**
 * Created by jmmut on 2015-10-28.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporter {

    private static final Logger logger = LoggerFactory.getLogger(VariantExporter.class);
    /**
     *
     * @param iterator
     * @param outputStream
     * @param sourceDBAdaptor
     * @param options TODO fill
     * @return num variants not written due to errors
     * @throws Exception
     */
    public static int VcfHtsExport(Iterator<Variant> iterator, OutputStream outputStream,
                                   VariantSourceDBAdaptor sourceDBAdaptor, QueryOptions options) throws IOException {
        final VCFHeader header = getVcfHeader(sourceDBAdaptor, options);
//        header.addMetaDataLine(new VCFFilterHeaderLine("PASS", "Valid variant"));
//        header.addMetaDataLine(new VCFFilterHeaderLine(".", "No FILTER info"));

        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();

        // setup writer
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        VariantContextWriter writer = builder
                .setOutputStream(outputStream)
                .setReferenceDictionary(sequenceDictionary)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        writer.writeHeader(header);

        // actual loop
        int failedVariants = 0;

        while (iterator.hasNext()) {
            Variant variant = iterator.next();
            try {
                VariantContext variantContext = convertBiodataVariantToVariantContext(variant);
                if (variantContext != null) {
                    writer.add(variantContext);
                }
            } catch (Exception e) {
                e.printStackTrace();
                failedVariants++;
            }
        }

        if (failedVariants > 0) {
            logger.warn(failedVariants + " variants were not written due to errors");
        }

        writer.close();
        return failedVariants;
    }

    private static VCFHeader getVcfHeader(VariantSourceDBAdaptor sourceDBAdaptor, QueryOptions options) throws IOException {

        List<VariantSource> sources = sourceDBAdaptor.getAllSources(options).getResult();
        List<String> headers = new ArrayList<>();
        for (VariantSource source : sources) {

            Object headerObject = source.getMetadata().get(DBObjectToVariantSourceConverter.HEADER_FIELD);
            if (headerObject instanceof String) {
                headers.add((String)headerObject);
            }
        }
        //        get header from studyConfiguration

//        List<String> returnedSamples = new ArrayList<>();
//        if (options != null) {
//            returnedSamples = options.getAsStringList(VariantDBAdaptor.VariantQueryParams.RETURNED_SAMPLES.key());
//        }
        if (headers.size() < 1) {
            throw new IllegalArgumentException("file headers not available with options " + options.toString());
        }
        String fileHeader = headers.get(0);
/*
        int lastLineIndex = fileHeader.lastIndexOf("#CHROM");
        if (lastLineIndex >= 0) {
            String substring = fileHeader.substring(0, lastLineIndex);
            if (returnedSamples.isEmpty()) {
                BiMap<Integer, String> samplesPosition = StudyConfiguration.getSamplesPosition(studyConfiguration).inverse();
                returnedSamples = new ArrayList<>(samplesPosition.size());
                for (int i = 0; i < samplesPosition.size(); i++) {
                    returnedSamples.add(samplesPosition.get(i));
                }
            }
            String samples = String.join("\t", returnedSamples);
            logger.debug("export will be done on samples: [{}]", samples);

            fileHeader = substring + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples;
        }
//        return VariantFileMetadataToVCFHeaderConverter.parseVcfHeader(fileHeader);

        // */
        VCFCodec vcfCodec = new VCFCodec();
        LineIterator source = vcfCodec.makeSourceFromStream(new ByteArrayInputStream(fileHeader.getBytes()));
        FeatureCodecHeader featureCodecHeader = vcfCodec.readHeader(source);
        return (VCFHeader) featureCodecHeader.getHeaderValue();
    }


    /**
     * converts org.opencb.biodata.models.variant.Variant into a htsjdk.variant.variantcontext.VariantContext
     * some assumptions:
     * * splitted multiallelic variants will produce only one variantContexts. Merging is done
     * * If some normalization have been done to the variant, the source entries may have an attribute ORI like: "POS:REF:ALT_0(,ALT_N)*:ALT_IDX"
     * @param variant
     * @return
     */
    public static VariantContext convertBiodataVariantToVariantContext(Variant variant) {//, StudyConfiguration studyConfiguration) {
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        int start = variant.getStart();
        int end = variant.getEnd();
        String reference = variant.getReference();
        String alternate = variant.getAlternate();
        String filter = "PASS";
        List<String> allelesArray = Arrays.asList(reference, alternate);  // TODO jmmut: multiallelic
        ArrayList<Genotype> genotypes = new ArrayList<>();

        for (Map.Entry<String, VariantSourceEntry> source : variant.getSourceEntries().entrySet()) {
            for (Map.Entry<String, Map<String, String>> samplesData : source.getValue().getSamplesData().entrySet()) {
                // reminder of samplesData meaning: Map(sampleName -> Map(dataType -> value))
                String sampleName = samplesData.getKey();
                String gt = samplesData.getValue().get("GT");
//                System.out.println("gt = " + gt);
                if (gt != null) {
                    org.opencb.biodata.models.feature.Genotype genotype = new org.opencb.biodata.models.feature.Genotype(gt, reference, alternate);
                    List<Allele> alleles = new ArrayList<>();
//                    System.out.println("\tgenotype = " + genotype.getReference().toString());
//                    System.out.println("\tgenotype = " + genotype.getAlternate().toString());
                    for (int gtIdx : genotype.getAllelesIdx()) {
//                        System.out.println("\t\t>>"+originalAlleles);
//                        System.out.println("\t\t>>"+gtIdx);
//                        System.out.println("\t\t>>"+originalAlleles.get(gtIdx));

                        /*
                        if (gtIdx < originalAlleles.size() && gtIdx >= 0) {
                            alleles.add(Allele.create(originalAlleles.get(gtIdx), gtIdx == 0));    // allele is reference if the alleleIndex is 0
//                            alleles.add(Allele.create(allele, gtIdx == 0));    // allele is reference if the alleleIndex is 0
                        } else {
                            alleles.add(Allele.create(".", false)); // genotype of a secondary alternate, or an actual missing
                        }
                        */

                        if (gtIdx < allelesArray.size() && gtIdx >= 0) {
                            alleles.add(Allele.create(allelesArray.get(gtIdx), gtIdx == 0));    // allele is reference if the alleleIndex is 0
//                            alleles.add(Allele.create(allele, gtIdx == 0));    // allele is reference if the alleleIndex is 0
                        } else {
                            alleles.add(Allele.create(".", false)); // genotype of a secondary alternate, or an actual missing
                        }

                    }
                    genotypes.add(new GenotypeBuilder().name(sampleName).alleles(alleles).phased(genotype.isPhased()).make());
                }
            }
        }

        variantContextBuilder.chr(variant.getChromosome())
                .start(start)
                .stop(end)
                .id(String.join(",", variant.getIds()))
                .alleles(allelesArray)
                .genotypes()
                .filter(filter)
                .genotypes(genotypes);


        return variantContextBuilder.make();




//        if (reference.isEmpty()) {
//            try {
//                QueryResponse<QueryResult<GenomeSequenceFeature>> resultQueryResponse = cellbaseClient.getSequence(
//                        CellBaseClient.Category.genomic,
//                        CellBaseClient.SubCategory.region,
//                        Arrays.asList(Region.parseRegion(variant.getChromosome() + ":" + start + "-" + start)),
//                        new QueryOptions());
//                indelSequence = resultQueryResponse.getResponse().get(0).getResult().get(0).getSequence();
//                reference = indelSequence;
//                alternate = indelSequence + alternate;
//                end = start + reference.length() - 1;
////                if ((end - start) != reference.length()) {
////                    end = start + reference.length() - 1;
////                }
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//        if (alternate.isEmpty()) {
//            try {
//                start -= reference.length();
//                QueryResponse<QueryResult<GenomeSequenceFeature>> resultQueryResponse = cellbaseClient.getSequence(
//                        CellBaseClient.Category.genomic,
//                        CellBaseClient.SubCategory.region,
//                        Arrays.asList(Region.parseRegion(variant.getChromosome() + ":" + start + "-" + start)),
//                        new QueryOptions());
//                indelSequence = resultQueryResponse.getResponse().get(0).getResult().get(0).getSequence();
//                reference = indelSequence + reference;
//                alternate = indelSequence;
//                if ((end - start) != reference.length()) {
//                    end = start + reference.length() - 1;
//                }
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }


        ////////////////////
/*
        List<String> allelesArray = Arrays.asList(reference, alternate);  // TODO jmmut: multiallelic
        ArrayList<Genotype> genotypes = new ArrayList<>();
        Integer originalPosition = null;
        List<String> originalAlleles = null;
        for (StudyEntry studyEntry : variant.getStudies()) {
            String[] ori = getOri(studyEntry);
            Integer auxOriginalPosition = getOriginalPosition(ori);
            if (originalPosition != null && auxOriginalPosition != null && !originalPosition.equals(auxOriginalPosition)) {
                throw new IllegalStateException("Two or more VariantSourceEntries have different origin. Unable to merge");
            }
            originalPosition = auxOriginalPosition;
            originalAlleles = getOriginalAlleles(ori);
            if (originalAlleles == null) {
                originalAlleles = allelesArray;
            }

            //Only print those variants in which the alternate is the first alternate from the multiallelic alternatives
            if (originalAlleles.size() > 2 && !"0".equals(getOriginalAlleleIndex(ori))) {
                logger.debug("Skip multiallelic variant! " + variant);
                return null;
            }


            String sourceFilter = studyEntry.getAttribute("FILTER");
            if (sourceFilter != null && !filter.equals(sourceFilter)) {
                filter = ".";   // write PASS iff all sources agree that the filter is "PASS" or assumed if not present, otherwise write "."
            }

            for (String sampleName : studyEntry.getOrderedSamplesName()) {
                Map<String, String> sampleData = studyEntry.getSampleData(sampleName);
                String gt = sampleData.get("GT");
//                System.out.println("gt = " + gt);
                if (gt != null) {
                    org.opencb.biodata.models.feature.Genotype genotype = new org.opencb.biodata.models.feature.Genotype(gt, reference, alternate);
                    List<Allele> alleles = new ArrayList<>();
//                    System.out.println("\tgenotype = " + genotype.getReference().toString());
//                    System.out.println("\tgenotype = " + genotype.getAlternate().toString());
                    for (int gtIdx : genotype.getAllelesIdx()) {
//                        System.out.println("\t\t>>"+originalAlleles);
//                        System.out.println("\t\t>>"+gtIdx);
//                        System.out.println("\t\t>>"+originalAlleles.get(gtIdx));
//                        String allele = originalAlleles.get(gtIdx);
//                        if (allele == null || allele.isEmpty()) {
//                            allele = "A";
//                        }
                        if (gtIdx < originalAlleles.size() && gtIdx >= 0) {
                            alleles.add(Allele.create(originalAlleles.get(gtIdx), gtIdx == 0));    // allele is reference if the alleleIndex is 0
//                            alleles.add(Allele.create(allele, gtIdx == 0));    // allele is reference if the alleleIndex is 0
                        } else {
                            alleles.add(Allele.create(".", false)); // genotype of a secondary alternate, or an actual missing
                        }
                    }
                    genotypes.add(new GenotypeBuilder().name(sampleName).alleles(alleles).phased(genotype.isPhased()).make());
                }
            }
        }

        List<String> ids = variant.getIds();

        variantContextBuilder.start(originalPosition == null ? start : originalPosition)
                .stop((originalPosition == null ? start : originalPosition) + (originalAlleles == null? allelesArray : originalAlleles).get(0).length() - 1)
                .chr(variant.getChromosome())
                .alleles(originalAlleles == null? allelesArray : originalAlleles)
                .filter(filter)
                .genotypes(genotypes);
//                .attributes(variant.)// TODO jmmut: join attributes from different source entries? what to do on a collision?


        if (ids != null) {
            Optional<String> reduce = variant.getIds().stream().reduce((left, right) -> left + "," + right);
            if (reduce.isPresent()) {
                variantContextBuilder.id(reduce.get());
            }
        }

        return variantContextBuilder.make();
        */
    }

}
