/*
 * Copyright 2016 EMBL - European Bioinformatics Institute
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
import embl.ebi.variation.eva.vcfdump.regionutils.RegionFactory;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.datastore.core.QueryOptions;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StorageManagerFactory;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;
import org.opencb.opencga.storage.core.variant.adaptors.VariantSourceDBAdaptor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.ws.rs.core.MultivaluedMap;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by jmmut on 2016-01-20.
 *
 * @author Jose Miguel Mut Lopez &lt;jmmut@ebi.ac.uk&gt;
 */
public class VariantExporterController {

    private static final Logger logger = LoggerFactory.getLogger(VariantExporterController.class);
    private static final int WINDOW_SIZE = 20000;

    private final CellbaseWSClient cellBaseClient;
    private final List<String> studies;
    private final List<String> files;
    private final String outputDir;
    private final VariantSourceDBAdaptor variantSourceDBAdaptor;
    private final VariantDBAdaptor variantDBAdaptor;
    private final QueryOptions query;
   // private final RegionDivider regionDivider;
    private final RegionFactory regionFactory;
    private final VariantExporter exporter;
    private Map<String, String> outputFilePaths;
    private int failedVariants;

    public VariantExporterController(String species, String dbName, List<String> studies, List<String> files,
                                     String outputDir, MultivaluedMap<String, String> queryParameters)
            throws IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException, URISyntaxException {
        this.studies = studies;
        this.files = files;
        assert outputDir != null;
        this.outputDir = outputDir;
        cellBaseClient = new CellbaseWSClient(species);
        variantDBAdaptor = getVariantDBAdaptor(dbName);
        query = getQuery(queryParameters);
        variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        //this.regionDivider = new RegionDivider(WINDOW_SIZE);
        regionFactory = new RegionFactory(WINDOW_SIZE, variantDBAdaptor, query);
        exporter = new VariantExporter(cellBaseClient);
        failedVariants = 0;
    }

    public VariantDBAdaptor getVariantDBAdaptor(String dbName) throws StorageManagerException, IllegalAccessException, ClassNotFoundException, InstantiationException {
        assert dbName != null;
        VariantStorageManager variantStorageManager = StorageManagerFactory.getVariantStorageManager();
        return variantStorageManager.getDBAdaptor(dbName, null);
    }

    public QueryOptions getQuery(MultivaluedMap<String, String> queryParameters) {
        QueryOptions query = new QueryOptions();
        query.put(VariantDBAdaptor.FILES, files);
        query.put(VariantDBAdaptor.STUDIES, studies);

        queryParameters.forEach((parameterName, parameterValues) -> {
            if (VariantDBAdaptor.QueryParams.acceptedValues.contains(parameterName)) {
                query.add(parameterName, String.join(",", parameterValues));
            }
        });

        return query;
    }

    public void run() {
        Map<String, VCFHeader> headers = getVcfHeaders();
        if (headers != null) {
            Map<String, VariantContextWriter> writers = getWriters(headers);
            writers.forEach((study, writer) -> writer.writeHeader(headers.get(study)));

            // get all chromosomes in the query or organism, and export the variants for each chromosome
            Set<String> chromosomes = getChromosomes(headers, studies, query);
            for (String chromosome : chromosomes) {
                exportChromosomeVariants(writers, chromosome);
            }

            writers.values().forEach(VariantContextWriter::close);
        }
    }

    private Map<String, VCFHeader> getVcfHeaders() {
        // get VCF header(s) and write them to output file(s)
        logger.info("Generating VCF headers ...");
        Map<String, VariantSource> sources = exporter.getSources(variantSourceDBAdaptor, studies);
        Map<String, VCFHeader> headers = null;
        try {
            headers = exporter.getVcfHeaders(sources);
            checkHeaders(headers);
        } catch (IOException e) {
            logger.error("Error getting vcf headers: {}", e.getMessage());
        }
        return headers;
    }

    private void exportChromosomeVariants(Map<String, VariantContextWriter> writers, String chromosome) {
        List<Region> allRegionsInChromosome = regionFactory.getRegionsForChromosome(chromosome);
        for (Region region : allRegionsInChromosome) {
            VariantDBIterator regionVariantsIterator = variantDBAdaptor.iterator(getRegionQuery(region));
            Map<String, List<VariantContext>> exportedVariants = exporter.export(regionVariantsIterator, region, studies);
            failedVariants += exporter.getFailedVariants();
            writeRegionVariants(writers, exportedVariants);
        }
    }

    private void writeRegionVariants(Map<String, VariantContextWriter> writers, Map<String, List<VariantContext>> exportedVariants) {
        exportedVariants.forEach((study, variants) -> {
            if (writers.containsKey(study)) {
                Collections.sort(variants, (v1, v2)-> v1.getStart() - v2.getStart());
                VariantContextWriter writer = writers.get(study);
                variants.forEach(writer::add);
            }
        });
    }

    private QueryOptions getRegionQuery(Region region) {
        QueryOptions regionQuery = new QueryOptions(query);
        regionQuery.put(VariantDBAdaptor.REGION, region.toString());
        return regionQuery;
    }

    private void checkHeaders(Map<String, VCFHeader> headers) {
        for (String studyId : studies) {
            if (!headers.containsKey(studyId)) {
                throw new IllegalArgumentException("Aborting VCF export: missing header for study " + studyId);
            }
        }
    }

    private Map<String, VariantContextWriter> getWriters(Map<String, VCFHeader> headers) {
        String suffix = ".exported.vcf.gz";
        Map<String, VariantContextWriter> writers = new TreeMap<>();
        outputFilePaths = new HashMap<>();

        // setup writers
        for (String studyId : studies) {

            Path outputPath = Paths.get(outputDir).resolve(studyId + suffix);
            outputFilePaths.put(studyId, outputPath.toString());
            // get sequence dictionary from header
            SAMSequenceDictionary sequenceDictionary;
            try {
                sequenceDictionary = headers.get(studyId).getSequenceDictionary();
            } catch (SAMException e) {
                logger.warn("Incorrect sequence / contig meta-data found in study {}: {}", studyId, e.getMessage());
                logger.warn("It won't be included in study {} output VCF header", studyId);
                sequenceDictionary = null;
            }

            VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
            VariantContextWriter writer = builder.setOutputFile(outputPath.toFile())
                    .setReferenceDictionary(sequenceDictionary)
                    .unsetOption(Options.INDEX_ON_THE_FLY)
                    .build();
            writers.put(studyId, writer);
        }

        //logger.info("Exporting to files: [" + String.join(" ", files) + "]");

        return writers;
    }

    private Set<String> getChromosomes(Map<String, VCFHeader> headers, List<String> studyIds, QueryOptions options) {
        Set<String> chromosomes;

        List<String> regions = options.getAsStringList(VariantDBAdaptor.REGION);
        if (regions.size() > 0) {
            chromosomes = getChromosomesFromRegionFilter(regions);
        } else {
            chromosomes = cellBaseClient.getChromosomes();
            if (chromosomes == null || chromosomes.isEmpty()) {
                chromosomes = getChromosomesFromVCFHeader(headers, studyIds);
            }
        }
        logger.debug("Chromosomes: {}", String.join(", ", chromosomes));
        return chromosomes;
    }

    private Set<String> getChromosomesFromRegionFilter(List<String> regions) {
        return regions.stream().map(r -> r.split(":")[0]).collect(Collectors.toSet());
    }

    private Set<String> getChromosomesFromVCFHeader(Map<String, VCFHeader> headers, List<String> studyIds) {
        Set<String> chromosomes = new HashSet<>();
        // setup writers
        for (String studyId : studyIds) {
            SAMSequenceDictionary sequenceDictionary = headers.get(studyId).getSequenceDictionary();
            chromosomes.addAll(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet()));
        }

        return chromosomes;
    }

    public Map<String, String> getOuputFiles() {
        return outputFilePaths;
    }

    public int getFailedVariants() {
        return failedVariants;
    }
}
