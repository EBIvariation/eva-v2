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
package embl.ebi.variation.eva.vcfDump;

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
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.opencb.opencga.storage.mongodb.variant.DBObjectToVariantSourceConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by jmmut on 2015-10-28.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporter {

    private static final Logger logger = LoggerFactory.getLogger(VariantExporter.class);
    private static int failedVariants = 0;

    /**
     *
     * @param iterator where to get the variants from
     * @param outputDir directory to write the output vcf(s).
     * @param sourceDBAdaptor to retrieve all the VariantSources in any VariantSourceEntry.
     * @param options not implemented yet, use only for studyId and fileId
     * @return num variants not written due to errors
     */
    public static List<String> VcfHtsExport(Iterator<Variant> iterator, String outputDir,
                                            VariantSourceDBAdaptor sourceDBAdaptor, QueryOptions options) throws IOException {


        List<String> studyIds = options.getAsStringList(VariantDBAdaptor.STUDIES);
        List<VariantSource> sourcesList = sourceDBAdaptor.getAllSourcesByStudyIds(studyIds, options).getResult();

        Map<String, VariantSource> sources = new TreeMap<>();

        for (VariantSource variantSource : sourcesList) {
            sources.put(variantSource.getStudyId(), variantSource);
        }
        
        Map<String, VCFHeader> headers = getVcfHeaders(sources);
//        header.addMetaDataLine(new VCFFilterHeaderLine("PASS", "Valid variant"));
//        header.addMetaDataLine(new VCFFilterHeaderLine(".", "No FILTER info"));

        String suffix = ".exported.vcf.gz";
        List<String> files = new ArrayList<>();
        Map<String, VariantContextWriter> writers = new TreeMap<>();

        // setup writers
        for (String studyId : studyIds) {
            if (!headers.containsKey(studyId)) {
                throw new IllegalArgumentException("aborting VCF export: missing header for study " + studyId);
            }
            VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
            File outFile = Paths.get(outputDir).resolve(studyId + suffix).toFile();
            files.add(outFile.getPath());
            VariantContextWriter writer = builder
                    .setOutputFile(outFile)
                    .setReferenceDictionary(headers.get(studyId).getSequenceDictionary())
                    .unsetOption(Options.INDEX_ON_THE_FLY)
                    .build();
            writers.put(studyId, writer);
            writer.writeHeader(headers.get(studyId));
        }

        // actual loop
        int failedVariants = 0;
        StringBuilder filesConcat = new StringBuilder();
        for (String file : files) {
            filesConcat.append(file).append(" ");
        }
        logger.info("Exporting to files: [ " + filesConcat.toString() + "]");

        while (iterator.hasNext()) {
            Variant variant = iterator.next();
            try {
                Map<String, VariantContext> variantContexts = convertBiodataVariantToVariantContext(variant, sources);
                for (Map.Entry<String, VariantContext> variantContextEntry : variantContexts.entrySet()) {
                    if (writers.containsKey(variantContextEntry.getKey())) {
                        writers.get(variantContextEntry.getKey()).add(variantContextEntry.getValue());
                    }
                }
            } catch (Exception e) {
                logger.info("failed variant: ", e);
                failedVariants++;
            }
        }

        if (failedVariants > 0) {
            logger.warn(failedVariants + " variants were not written due to errors");
        }
        VariantExporter.failedVariants += failedVariants;

        for (VariantContextWriter variantContextWriter : writers.values()) {
            variantContextWriter.close();
        }

        return files;
    }

    private static Map<String, VCFHeader> getVcfHeaders(Map<String, VariantSource> sources) throws IOException {

        Map<String, VCFHeader> headers = new TreeMap<>();
        
        for (VariantSource source : sources.values()) {
            Object headerObject = source.getMetadata().get(DBObjectToVariantSourceConverter.HEADER_FIELD);

            if (headerObject instanceof String) {
                VCFCodec vcfCodec = new VCFCodec();
                ByteArrayInputStream bufferedInputStream = new ByteArrayInputStream(((String) headerObject).getBytes());
                LineIterator sourceFromStream = vcfCodec.makeSourceFromStream(bufferedInputStream);
                FeatureCodecHeader featureCodecHeader = vcfCodec.readHeader(sourceFromStream);
                headers.put(source.getStudyId(), (VCFHeader) featureCodecHeader.getHeaderValue());
            }
        }
        
        for (String studyId : sources.keySet()) {
            if (!headers.containsKey(studyId)) {
                throw new IllegalArgumentException("file headers not available for study " + studyId);
            }
        }

        return headers;
        //TODO: allow specify which samples to return
/*
//        List<String> returnedSamples = new ArrayList<>();
//        if (options != null) {
//            returnedSamples = options.getAsStringList(VariantDBAdaptor.VariantQueryParams.RETURNED_SAMPLES.key());
//        }
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

    }

    public static int getFailedVariants() {
        return failedVariants;
    }

    private static class VariantFields {
        public int start;
        public int end;
        public String reference;
        public String alternate;
    }

    /**
     * converts org.opencb.biodata.models.variant.Variant into one or more htsjdk.variant.variantcontext.VariantContext
     * behaviour:
     * * one VariantContext per study
     * * split multiallelic variant will remain split.
     * * in case a normalized INDEL has empty alleles, the original alleles in the vcf line will be used.
     *
     * steps:
     * * foreach variantSourceEntry, collect genotypes in its study, only if the study was requested
     * * get main variant data: position, alleles, filter...
     * * if there are empty alleles, get them from the vcf line
     * * get the genotypes
     * * add all (position, alleles, genotypes...) to a VariantContext for each study.
     *
     * @param variant
     * @return
     */
    public static Map<String, VariantContext> convertBiodataVariantToVariantContext(
            Variant variant, Map<String, VariantSource> sources) 
            throws IOException {

        List<String> studyIds = new ArrayList<>(sources.keySet());
        Map<String, VariantContext> variantContextMap = new TreeMap<>();
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        Integer start = variant.getStart();
        Integer end = variant.getEnd();
        String reference = variant.getReference();
        String alternate = variant.getAlternate();
        String filter = "PASS";
        List<String> allelesArray = Arrays.asList(reference, alternate);
        ArrayList<Genotype> genotypes = new ArrayList<>();
        Map<String, List<Genotype>> genotypesPerStudy = new TreeMap<>();
        
        for (Map.Entry<String, VariantSourceEntry> source : variant.getSourceEntries().entrySet()) {

            String studyId = source.getValue().getStudyId();

            if (studyIds.contains(studyId)) {   // skipping studies not asked
                
                // if we added this outside the loop, if the study is not present in this variant, the writer would add 
                // a whole line of "./."
                if (!genotypesPerStudy.containsKey(studyId)) {  
                    genotypesPerStudy.put(studyId, new ArrayList<Genotype>());
                }
                
                // if there are indels, we cannot use the normalized alleles, (hts forbids empty alleles) so we have to take them from the original vcf line
                int emtpyAlleles = 0;
                for (String a : allelesArray) {
                    emtpyAlleles += a.isEmpty() ? 1 : 0;
                }
                if (emtpyAlleles != 0) {
                    String src = source.getValue().getAttribute("src");
                    if (src != null) {
                        VariantSource variantSource = sources.get(studyId);
                        if (variantSource == null) {
                            throw new IllegalArgumentException("VariantSource not available for study " + studyId);
                        }
                        VariantFields variantFields = getVariantFields(variant, variantSource, src);

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
                        genotypesPerStudy.get(studyId).add(
                                new GenotypeBuilder().name(sampleName).alleles(alleles).phased(genotype.isPhased()).make());
                    }
                }
            }
        }


        for (Map.Entry<String, List<Genotype>> studyEntry : genotypesPerStudy.entrySet()) {
            VariantContext make = variantContextBuilder
                    .chr(variant.getChromosome())
                    .start(start)
                    .stop(end)
//                .id(String.join(";", variant.getIds()))   // in multiallelic, this results in duplicated ids, across several rows
                    .noID()
                    .alleles(allelesArray)
                    .filter(filter)
                    .genotypes(studyEntry.getValue()).make();
            variantContextMap.put(studyEntry.getKey(), make);
        }
        return variantContextMap;
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

        StringBuilder stringBuilder = new StringBuilder(alts[0]);
        for (int i = 1; i < alts.length; i++) {
            stringBuilder.append(",").append(alts[i]);
        }

        logger.debug("Using original alleles from vcf line in \"{}_{}_{}_{}\". Original ref and alts: \"{}:{}\". Output: \"{}:{}\"",
                variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate(),
                variantFields.reference, stringBuilder.toString(),
                variantFields.reference, variantFields.alternate);
        return variantFields;
    }

}
