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
import org.apache.commons.lang.StringUtils;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.*;
import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.cellbase.core.common.GenomeSequenceFeature;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.datastore.core.QueryResponse;
import org.opencb.datastore.core.QueryResult;
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

    private CellBaseClient cellbaseClient;
    /**
     * Read only. Keeps track of the total failed variants across several dumps. To accumulate, use the same instance of
     * VariantExporter to dump several VCFs. If you just want to count on one dump, use a `new VariantExporter` each time.
     */
    private int failedVariants = 0;

    /**
     * if the variants will have empty alleles (such as normalized deletions: "A" to "") CellBase is mandatory to
     * retrieve the context, (for instance "GA" to "G").
     *
     * This context may be different from the alleles in the original VCF.
     *
     * If there won't be empty alleles, CellBase is not needed, and the parameter may be null.
     *
     * @param cellbaseClient for empty alleles. nullable.
     */
    public VariantExporter(CellBaseClient cellbaseClient) {
        this.cellbaseClient = cellbaseClient;
    }

    /**
     * Main method of this class. It generates one VCF file per requested study,
     * writing both the header meta-data and the variants in the body.
     * 
     * @param iterator where to get the variants from
     * @param outputDir directory to write the output vcf(s) to
     * @param sourceDBAdaptor to retrieve all the VariantSources in any VariantSourceEntry
     * @param options not implemented yet, use only for studyId and fileId
     * @return list of files written, one per study.
     */
    public List<String> VcfHtsExport(Iterator<Variant> iterator, String outputDir,
                                     VariantSourceDBAdaptor sourceDBAdaptor, QueryOptions options) throws IOException {

        // 3 steps to get the headers of all the studyIds: ask the sourceAdaptor for VariantSources, check we got everything, build the headers.
        List<String> studyIds = options.getAsStringList(VariantDBAdaptor.STUDIES);

        // 1) retrieve the sources
        Map<String, VariantSource> sources = new TreeMap<>();
        List<VariantSource> sourcesList = sourceDBAdaptor.getAllSourcesByStudyIds(studyIds, options).getResult();
        for (VariantSource variantSource : sourcesList) {
            sources.put(variantSource.getStudyId(), variantSource);
        }

        // 2) check that sourceDBAdaptor got all the studyIds
        for (String studyId : studyIds) {
            if (!sources.containsKey(studyId)) {
                throw new IllegalArgumentException("Aborting VCF export: missing header for study " + studyId);
            }
        }

        // 3) check and get the headers, one for each source
        Map<String, VCFHeader> headers = getVcfHeaders(sources);

        // from here we grant that `headers` have all the headers requested in `studyIds`

        String suffix = ".exported.vcf.gz";
        List<String> files = new ArrayList<>();
        Map<String, VariantContextWriter> writers = new TreeMap<>();

        // setup writers
        for (String studyId : studyIds) {
            VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
            File outFile = Paths.get(outputDir).resolve(studyId + suffix).toFile();
            files.add(outFile.getPath());
            SAMSequenceDictionary sequenceDictionary;
            try {
                sequenceDictionary = headers.get(studyId).getSequenceDictionary();
            } catch (Exception e) {
                sequenceDictionary = null;
            }
            VariantContextWriter writer = builder
                    .setOutputFile(outFile)
                    .setReferenceDictionary(sequenceDictionary)
                    .unsetOption(Options.INDEX_ON_THE_FLY)
                    .build();
            writers.put(studyId, writer);
            writer.writeHeader(headers.get(studyId));
        }

        logger.info("Exporting to files: [" + StringUtils.join(files, " ") + "]");

        // actual loop
        int failedVariants = 0;
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
                logger.info(String.format("Variant dump failed: \"%s:%d:%s>%s\"", variant.getChromosome(),
                                variant.getStart(), variant.getReference(), variant.getAlternate()),
                        e);
                failedVariants++;
            }
        }

        if (failedVariants > 0) {
            logger.warn(failedVariants + " variants were not written due to errors");
        }
        this.failedVariants += failedVariants;

        for (VariantContextWriter variantContextWriter : writers.values()) {
            variantContextWriter.close();
        }

        return files;
    }

    /**
     * postconditions:
     * - returns one header per study (one header for each key in `sources`).
     * @throws IOException
     */
    private Map<String, VCFHeader> getVcfHeaders(Map<String, VariantSource> sources) throws IOException {
        Map<String, VCFHeader> headers = new TreeMap<>();
        
        for (VariantSource source : sources.values()) {
            Object headerObject = source.getMetadata().get(DBObjectToVariantSourceConverter.HEADER_FIELD);

            if (headerObject instanceof String) {
                VCFCodec vcfCodec = new VCFCodec();
                ByteArrayInputStream bufferedInputStream = new ByteArrayInputStream(((String) headerObject).getBytes());
                LineIterator sourceFromStream = vcfCodec.makeSourceFromStream(bufferedInputStream);
                FeatureCodecHeader featureCodecHeader = vcfCodec.readHeader(sourceFromStream);
                VCFHeader headerValue = (VCFHeader) featureCodecHeader.getHeaderValue();
                headerValue.addMetaDataLine(new VCFFilterHeaderLine("PASS", "Valid variant"));
                headers.put(source.getStudyId(), headerValue);
            } else {
                throw new IllegalArgumentException("File headers not available for study " + source.getStudyId());
            }
        }

        return headers;
        //TODO: allow specify which samples to return
/*
//
//        header.addMetaDataLine(new VCFFilterHeaderLine(".", "No FILTER info"));

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

    public int getFailedVariants() {
        return failedVariants;
    }

    class VariantFields {
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
     * * in case a normalized INDEL has empty alleles, a query to cellbase will be done, to have the previous base as context
     *
     * steps:
     * * foreach variantSourceEntry, collect genotypes in its study, only if the study was requested
     * * get main variant data: position, alleles, filter...
     * * if there are empty alleles, get them from cellbase
     * * get the genotypes
     * * add all (position, alleles, genotypes...) to a VariantContext for each study.
     *
     * @param variant
     * @return
     */
    public Map<String, VariantContext> convertBiodataVariantToVariantContext(
            Variant variant, Map<String, VariantSource> sources) throws IOException {
        int missingGenotypes = 0;
        Set<String> studyIds = sources.keySet();
        Map<String, VariantContext> variantContextMap = new TreeMap<>();
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        Integer start = variant.getStart();
        Integer end = variant.getEnd();
        String reference = variant.getReference();
        String alternate = variant.getAlternate();
        String filter = "PASS";
        List<String> allelesArray = Arrays.asList(reference, alternate);
        Map<String, List<Genotype>> genotypesPerStudy = new TreeMap<>();

        for (VariantSourceEntry source : variant.getSourceEntries().values()) {

            String studyId = source.getStudyId();

            if (studyIds.contains(studyId)) {   // skipping studies not asked

                // if we added this outside the loop, if the study is not present in this variant, the writer would add
                // a whole line of "./."
                if (!genotypesPerStudy.containsKey(studyId)) {
                    genotypesPerStudy.put(studyId, new ArrayList<Genotype>());
                }

                // if there are indels, we cannot use the normalized alleles, (hts forbids empty alleles) so we have to take them from cellbase
                boolean emptyAlleles = false;
                for (String a : allelesArray) {
                    if (a.isEmpty()) {
                        emptyAlleles = true;
                        break;
                    }
                }

                if (emptyAlleles) {
                    if (cellbaseClient != null) {
                        List<Region> regions = Collections.singletonList(new Region(variant.getChromosome(), start - 1, start-1));

                        QueryResponse<QueryResult<GenomeSequenceFeature>> sequence = cellbaseClient.getSequence(
                                CellBaseClient.Category.genomic, CellBaseClient.SubCategory.region, regions, null);

                        List<GenomeSequenceFeature> response = sequence.getResponse().get(0).getResult();
                        if (response.size() == 1) {
                            start--;
                            allelesArray = new ArrayList<>();
                            allelesArray.add(response.get(0).getSequence() + reference);
                            allelesArray.add(response.get(0).getSequence() + alternate);
                            end = start + allelesArray.get(0).length()-1;   // -1 because end is inclusive. [start, end] instead of [start, end)
                        }
                    } else {
                        throw new IllegalArgumentException(String.format(
                                        "CellBase was not provided, needed to fill empty alleles at study %s, in variant %s:%d:%s>%s", studyId,
                                        variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
                    }
                }

                // add the genotypes
                for (Map.Entry<String, Map<String, String>> samplesData : source.getSamplesData().entrySet()) {
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
                    } else {
                        missingGenotypes++;
                    }
                }
            }
        }

        if (missingGenotypes > 0) {
            logger.info("Variant %s:%d:%s>%s lacked the GT field in %d genotypes (they will be printed as \"./.\").", 
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate());
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
     * In case there is an INDEL (multiallelic or not), we have to retrieve the alleles from the original VCF line.
     * 
     * @param variant the variant as available in the variant store
     * @param variantSource VCF file meta-data as available in the variant store
     * @param srcLine the original VCF line as available in the variant store
     * @return the variant coordinates and alleles to write as output
     */
    private VariantFields getVariantFields(Variant variant, VariantSource variantSource, String srcLine) {
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
                    "Variant %s:%d:%s>%s has empty alleles and no original VCF line stored",
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
        }

        VariantFields variantFields = new VariantFields();
        variantFields.reference = split[3];
        variantFields.alternate = alts[alleleNumber];
        variantFields.start = Integer.parseInt(split[1]);
        variantFields.end = variantFields.start + variantFields.reference.length()-1;

        logger.debug("Using alleles from original VCF line in {}:{}:{}>{}. Original REF and ALT: \"{}>{}\". Output: \"{}>{}\"",
                variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate(),
                variantFields.reference, StringUtils.join(alts, ","),
                variantFields.reference, variantFields.alternate);
        return variantFields;
    }

}
