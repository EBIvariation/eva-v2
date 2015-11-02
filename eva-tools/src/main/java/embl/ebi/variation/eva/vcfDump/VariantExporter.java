package embl.ebi.variation.eva.vcfDump;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.opencb.biodata.models.variant.*;
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
     * @param outputStream where to write the output vcf. you can use a GZIPOutputStrem for compressed files
     * @param sourceDBAdaptor to retrieve all the VariantSources in any VariantSourceEntry.
     * @param options not implemented yet, use only for studyId and fileId
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
                VariantContext variantContext = convertBiodataVariantToVariantContext(variant, sourceDBAdaptor, options);
                if (variantContext != null) {
                    writer.add(variantContext);
                }
            } catch (Exception e) {
                logger.info("failed variant: ", e);
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

        //TODO: allow specify which samples to return
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

    private static class VariantFields {
        public int start;
        public int end;
        public String reference;
        public String alternate;
    }

    /**
     * converts org.opencb.biodata.models.variant.Variant into a htsjdk.variant.variantcontext.VariantContext
     * behaviour:
     * * split multiallelic variant will remain split.
     * * in case a normalized INDEL has empty alleles, the original alleles in the vcf line will be used.
     *
     * steps:
     * * get main variant data: position, alleles, filter...
     * * if there are empty alleles, get them from the vcf line
     * * get the genotypes
     * * add all (position, alleles, genotypes...) to a VariantContext.
     *
     * @param variant
     * @return
     */
    public static VariantContext convertBiodataVariantToVariantContext(Variant variant, VariantSourceDBAdaptor sourceDBAdaptor, QueryOptions sourceOptions)
            throws IOException {

        VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        Integer start = variant.getStart();
        Integer end = variant.getEnd();
        String reference = variant.getReference();
        String alternate = variant.getAlternate();
        String filter = "PASS";
        List<String> allelesArray = Arrays.asList(reference, alternate);
        ArrayList<Genotype> genotypes = new ArrayList<>();

        for (Map.Entry<String, VariantSourceEntry> source : variant.getSourceEntries().entrySet()) {

            // if there are indels, we cannot use the normalized alleles, (hts forbids empty alleles) so we have to take them from the original vcf line
            int emtpyAlleles = 0;
            for (String a : allelesArray) {
                emtpyAlleles += a.isEmpty()? 1: 0;
            }
            if (emtpyAlleles != 0) {
                String src = source.getValue().getAttribute("src");
                if (src != null) {
                    sourceOptions.add("studyId", source.getValue().getStudyId());
                    sourceOptions.add("fileId", source.getValue().getFileId());
                    QueryResult<VariantSource> allSources = sourceDBAdaptor.getAllSources(sourceOptions);   // TODO cache this
                    if (allSources.getResult().size() < 1) {
                        throw new IllegalArgumentException("VariantSource not available with options " + sourceOptions.toString());
                    }
                    VariantFields variantFields = getVariantFields(variant, allSources.getResult().get(0), src);

                    // overwrite the initial-guess position and alleles
                    allelesArray = new ArrayList<>();
                    allelesArray.add(variantFields.reference);
                    allelesArray.add(variantFields.alternate);
                    start = variantFields.start;
                    end = variantFields.end;
                }
            }

            // add the genotypes
            for (Map.Entry<String, Map<String, String>> samplesData : source.getValue().getSamplesData().entrySet()) {
                // reminder of samplesData meaning: Map(sampleName -> Map(dataType -> value))
                String sampleName = samplesData.getKey();
                String gt = samplesData.getValue().get("GT");

                if (gt != null) {
                    org.opencb.biodata.models.feature.Genotype genotype = new org.opencb.biodata.models.feature.Genotype(gt, reference, alternate);
                    List<Allele> alleles = new ArrayList<>();
                    for (int gtIdx : genotype.getAllelesIdx()) {
                        if (gtIdx < allelesArray.size() && gtIdx >= 0) {
                            alleles.add(Allele.create(allelesArray.get(gtIdx), gtIdx == 0));    // allele is reference if the alleleIndex is 0
                        } else {
                            alleles.add(Allele.create(".", false)); // genotype of a secondary alternate, or an actual missing
                        }
                    }
                    genotypes.add(new GenotypeBuilder().name(sampleName).alleles(alleles).phased(genotype.isPhased()).make());
                }
            }
        }

        variantContextBuilder
                .chr(variant.getChromosome())
                .start(start)
                .stop(end)
                .id(String.join(";", variant.getIds()))
                .alleles(allelesArray)
                .genotypes()
                .filter(filter)
                .genotypes(genotypes);

        return variantContextBuilder.make();
    }

    /**
     * In case there is an INDEL (multiallelic or not), we have to retrieve the alleles from the original vcf line.
     * @param variant
     * @param variantSource
     * @param srcLine  @return
     */
    private static VariantFields getVariantFields(Variant variant, VariantSource variantSource, String srcLine) {
        String[] split = srcLine.split("\t", 6);
        StringBuilder newLineBuilder = new StringBuilder();
        for (int i = 0; i < split.length - 1; i++) {
            newLineBuilder.append(split[i]).append("\t");
        }
        newLineBuilder.append(".\t.\t.");   // ignoring qual, filter, info, format and genotypes. We just want normalization
        VariantFactory variantFactory = new VariantVcfFactory();
        List<Variant> variants = variantFactory.create(variantSource, newLineBuilder.toString());

        int alleleNumber = 0;
        for (Variant v : variants) {
            // reference should be always equal; this may be shortened
            if (v.getReference().equals(variant.getReference()) && v.getAlternate().equals(variant.getAlternate())) {
                break;
            }
            alleleNumber++;
        }

        String[] alts = split[4].split(",");
        if (alleleNumber >= alts.length) {
            throw new IllegalArgumentException(String.format(
                    "Variant \"%s_%s_%s_%s\" has empty alleles and no original line",
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
        }

        VariantFields variantFields = new VariantFields();
        variantFields.reference = split[3];
        variantFields.alternate = alts[alleleNumber];
        variantFields.start = Integer.parseInt(split[1]);
        variantFields.end = variantFields.start + variantFields.reference.length()-1;

        logger.debug("Using original alleles from vcf line in \"{}_{}_{}_{}\". Original ref and alts: \"{}:{}\". Output: \"{}:{}\"",
                variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate(),
                variantFields.reference, String.join(",", alts),
                variantFields.reference, variantFields.alternate);
        return variantFields;
    }

}
