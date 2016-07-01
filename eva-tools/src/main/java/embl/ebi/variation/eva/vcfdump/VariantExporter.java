/*
 * Copyright 2015-2016 EMBL - European Bioinformatics Institute
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
package embl.ebi.variation.eva.vcfdump;

import embl.ebi.variation.eva.vcfdump.cellbasewsclient.CellbaseWSClient;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantSourceEntry;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.opencb.opencga.storage.mongodb.variant.DBObjectToVariantSourceConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.*;

/**
 * Created by jmmut on 2015-10-28.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporter {

    private static final Logger logger = LoggerFactory.getLogger(VariantExporter.class);

    private CellbaseWSClient cellbaseClient;
    /**
     * Read only. Keeps track of the total failed variants across several dumps. To accumulate, use the same instance of
     * VariantExporter to dump several VCFs. If you just want to count on one dump, use a `new VariantExporter` each time.
     */
    private int failedVariants;
    private String regionSequence;

    /**
     * if the variants will have empty alleles (such as normalized deletions: "A" to "") CellBase is mandatory to
     * retrieve the context, (for instance "GA" to "G").
     *
     * This context may be different from the alleles in the original VCF.
     *
     * If there won't be empty alleles, CellBase is not needed, and the parameter may be null.
     *  @param cellbaseClient for empty alleles. nullable.
     */
    public VariantExporter(CellbaseWSClient cellbaseClient) {
        this.cellbaseClient = cellbaseClient;
    }

    public Map<String, List<VariantContext>> export(VariantDBIterator iterator, Region region, List<String> studyIds) {
        // TODO: To be replaced by List<VariantContext>, because the variants wont be splitted by study
        Map<String, List<VariantContext>> variantsToExportByStudy = new HashMap<>();
        failedVariants = 0;

        // region sequence contains the last exported region: we set it to null to get the new region sequence from cellbase if needed
        regionSequence = null;

        while (iterator.hasNext()) {
            Variant variant = iterator.next();
            if (region.contains(variant.getChromosome(), variant.getStart())) {
                try {
                    // TODO: this call should be replaced by a call to BiodataVariantToVariantContextConverter
                    Map<String, VariantContext> variantContexts = convertBiodataVariantToVariantContext(variant, studyIds, region);
                    // TODO: we will get just one variantContest from BiodataVariantToVariantContextConverter
                    for (Map.Entry<String, VariantContext> variantContext : variantContexts.entrySet()) {
                        variantsToExportByStudy.putIfAbsent(variantContext.getKey(), new ArrayList<>());
                        variantsToExportByStudy.get(variantContext.getKey()).add(variantContext.getValue());
                    }
                } catch (Exception e) {
                    logger.warn("Variant {}:{}:{}>{} dump failed: {}", variant.getChromosome(), variant.getStart(), variant.getReference(),
                            variant.getAlternate(), e.getMessage());
                    failedVariants++;
                }
            }
        }
        return variantsToExportByStudy;
    }

    public Map<String, VariantSource> getSources(VariantSourceDBAdaptor sourceDBAdaptor, List<String> studyIds) {
        Map<String, VariantSource> sources = new TreeMap<>();
        List<VariantSource> sourcesList = sourceDBAdaptor.getAllSourcesByStudyIds(studyIds, new QueryOptions()).getResult();
        for (VariantSource variantSource : sourcesList) {
            sources.put(variantSource.getStudyId(), variantSource);
        }
        return sources;
    }

    /**
     * postconditions:
     * - returns one header per study (one header for each key in `sources`).
     * @throws IOException
     */
    public Map<String, VCFHeader> getVcfHeaders(Map<String, VariantSource> sources) throws IOException {
        Map<String, VCFHeader> headers = new TreeMap<>();

        for (VariantSource source : sources.values()) {
            Object headerObject = source.getMetadata().get(DBObjectToVariantSourceConverter.HEADER_FIELD);

            if (headerObject instanceof String) {
                VCFCodec vcfCodec = new VCFCodec();
                ByteArrayInputStream bufferedInputStream = new ByteArrayInputStream(((String) headerObject).getBytes());
                LineIterator sourceFromStream = vcfCodec.makeSourceFromStream(bufferedInputStream);
                FeatureCodecHeader featureCodecHeader = vcfCodec.readHeader(sourceFromStream);
                VCFHeader headerValue = (VCFHeader) featureCodecHeader.getHeaderValue();
                headerValue.addMetaDataLine(new VCFFilterHeaderLine("PASS", "All filters passed"));
                headers.put(source.getStudyId(), headerValue);
            } else {
                throw new IllegalArgumentException("File headers not available for study " + source.getStudyId());
            }
        }

        return headers;
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

    public VCFHeader getMergedVCFHeader(Map<String, VariantSource> sources) throws IOException {
        Map<String, VCFHeader> headers = getVcfHeaders(sources);

        Set<VCFHeaderLine> mergedHeaderLines = VCFUtils.smartMergeHeaders(headers.values(), true);
        // TODO: update sample names if there is conflict: the code is in BiodataVariantToVariantContextConverter and should be moved here
        return new VCFHeader(mergedHeaderLines);
    }

    public int getFailedVariants() {
        return failedVariants;
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
    @Deprecated
    public Map<String, VariantContext> convertBiodataVariantToVariantContext(
            Variant variant, List<String> studyIds, Region region) throws IOException {
        int missingGenotypes = 0;
        Map<String, VariantContext> variantContextMap = new TreeMap<>();
        VariantContextBuilder variantContextBuilder = new VariantContextBuilder();

        String reference = variant.getReference();
        String alternate = variant.getAlternate();
        Integer start = variant.getStart();
        Integer end = start + reference.length() -1;
        String filter = "PASS";
        String[] allelesArray = {reference, alternate};
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
                    // TODO: the context nucleotide is being retrieved several times for the same variant
                    // TODO: study this loop, extracting all common code to avoid being executed more than once if not necessary
                    // TODO: to be fixed in the merge studies feature
                    String contextNucleotide;
                    if (region != null) {
                        contextNucleotide = getContextNucleotideFromCellbaseCachingRegions(variant, start, region, studyId);
                    } else {
                        contextNucleotide = getContextNucleotideFromCellbase(variant, start, studyId);
                    }
                    // TODO: maybe this can be more efficient
                    allelesArray[0] = contextNucleotide + reference;
                    allelesArray[1] =  contextNucleotide + alternate;
                    end = start + allelesArray[0].length() - 1;
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
                            if (gtIdx < allelesArray.length && gtIdx >= 0) {
                                alleles.add(Allele.create(allelesArray[gtIdx], gtIdx == 0));    // allele is reference if the alleleIndex is 0
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

    private String getContextNucleotideFromCellbase(Variant variant, Integer start, String studyId) throws IOException {
        if (cellbaseClient != null) {
            return cellbaseClient.getSequence(new Region(variant.getChromosome(), start - 1, start-1));
        } else {
            throw new IllegalArgumentException(String.format(
                    "CellBase was not provided, needed to fill empty alleles at study %s, in variant %s:%d:%s>%s", studyId,
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
        }
    }

    private String getContextNucleotideFromCellbaseCachingRegions(Variant variant, int start, Region region, String studyId) throws IOException {
        if (cellbaseClient != null) {
            try {
                if (regionSequence == null) {
                    // if an indel start is the first nucleotide of the region, we will need the previous nucleotide, so we are adding
                    // the preceding nucleotide to the region (region.getStart()-1)
                    regionSequence = cellbaseClient.getSequence(new Region(variant.getChromosome(), region.getStart()-1, region.getEnd()));
                }
                String nucleotide = getNucleotideFromRegionSequence(start, region.getStart(), regionSequence);
                return nucleotide;
            } catch (Exception e) {
                logger.error("Getting nucleotide for variant {} in region {} using start {}: {}", variant, region, start, e.getMessage());
                throw e;
            }
        } else {
            throw new IllegalArgumentException(String.format(
                    "CellBase was not provided, needed to fill empty alleles at study %s, in variant %s:%d:%s>%s", studyId,
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
        }
    }

    private String getNucleotideFromRegionSequence(int start, int regionStart, String regionSequence) {
        int relativePosition = start - regionStart;
        return regionSequence.substring(relativePosition, relativePosition + 1);
    }


}
