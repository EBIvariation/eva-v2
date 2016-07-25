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
import java.time.LocalDateTime;
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
    private final String species;
    private final List<String> studies;
    private final List<String> files;
    private final String outputDir;
    private final VariantSourceDBAdaptor variantSourceDBAdaptor;
    private final VariantDBAdaptor variantDBAdaptor;
    private final QueryOptions query;
    private final RegionFactory regionFactory;
    private final VariantExporter exporter;
    private Path outputFilePath;
    private int failedVariants;

    public VariantExporterController(String species, String dbName, List<String> studies, String outputDir,
                                     MultivaluedMap<String, String> queryParameters)
            throws IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException, URISyntaxException {
        this(species, dbName, studies, null, outputDir, queryParameters);
    }

    public VariantExporterController(String species, String dbName, List<String> studies, List<String> files,
                                     String outputDir, MultivaluedMap<String, String> queryParameters)
            throws IllegalAccessException, ClassNotFoundException, InstantiationException, StorageManagerException, URISyntaxException {
        checkParams(species, studies, outputDir, dbName);
        this.species = species;
        this.studies = studies;
        this.files = files;
        this.outputDir = outputDir;
        cellBaseClient = new CellbaseWSClient(species);
        variantDBAdaptor = getVariantDBAdaptor(dbName);
        query = getQuery(queryParameters);
        variantSourceDBAdaptor = variantDBAdaptor.getVariantSourceDBAdaptor();
        regionFactory = new RegionFactory(WINDOW_SIZE, variantDBAdaptor, query);
        exporter = new VariantExporter(cellBaseClient);
        failedVariants = 0;
    }

    private void checkParams(String species, List<String> studies, String outputDir, String dbName) {
        if (species == null || species.isEmpty()) {
            throw new IllegalArgumentException("'species' is required");
        } else if (studies == null || studies.isEmpty()) {
            throw new IllegalArgumentException("'studies' is required");
        } else if (outputDir == null || outputDir.isEmpty()) {
            throw new IllegalArgumentException("'outputDir' is required");
        } else if (dbName == null || dbName.isEmpty()) {
            throw new IllegalArgumentException("'dbName' is required");
        }
    }

    public VariantDBAdaptor getVariantDBAdaptor(String dbName) throws StorageManagerException, IllegalAccessException, ClassNotFoundException, InstantiationException {
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
        VCFHeader header = getOutputVcfHeader();
        VariantContextWriter writer = getWriter(header);
        writer.writeHeader(header);

        // get all chromosomes in the query or organism, and export the variants for each chromosome
        Set<String> chromosomes = getChromosomes(header, studies, query);
        for (String chromosome : chromosomes) {
            exportChromosomeVariants(writer, chromosome);
        }

        writer.close();
    }

    private VCFHeader getOutputVcfHeader() {
        // get VCF header(s) and write them to output file(s)
        logger.info("Generating VCF headers ...");
        Map<String, VariantSource> sources = exporter.getSources(variantSourceDBAdaptor, studies);
        VCFHeader header = null;
        try {
            header = exporter.getMergedVCFHeader(sources);
        } catch (IOException e) {
            logger.error("Error getting vcf headers: {}", e.getMessage());
        }
        return header;
    }

    private void exportChromosomeVariants(VariantContextWriter writer, String chromosome) {
        List<Region> allRegionsInChromosome = regionFactory.getRegionsForChromosome(chromosome);
        for (Region region : allRegionsInChromosome) {
            VariantDBIterator regionVariantsIterator = variantDBAdaptor.iterator(getRegionQuery(region));
            List<VariantContext> exportedVariants = exporter.export(regionVariantsIterator, region);
            failedVariants += exporter.getFailedVariants();
            exportedVariants.forEach(writer::add);
        }
    }

    private QueryOptions getRegionQuery(Region region) {
        QueryOptions regionQuery = new QueryOptions(query);
        regionQuery.put(VariantDBAdaptor.REGION, region.toString());
        return regionQuery;
    }

    private VariantContextWriter getWriter(VCFHeader vcfHeader) {
        LocalDateTime now = LocalDateTime.now();
        String fileName = species + "_exported_" + now + ".vcf.gz";
        outputFilePath = Paths.get(outputDir).resolve(fileName);

        // get sequence dictionary from header
        SAMSequenceDictionary sequenceDictionary;
        try {
            sequenceDictionary = vcfHeader.getSequenceDictionary();
        } catch (SAMException e) {
            logger.warn("Incorrect sequence / contig meta-data: ", e.getMessage());
            logger.warn("It won't be included in output VCF header");
            sequenceDictionary = null;
        }

        // setup writer
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        VariantContextWriter writer = builder.setOutputFile(outputFilePath.toFile())
                .setReferenceDictionary(sequenceDictionary)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        return writer;
    }

    private Set<String> getChromosomes(VCFHeader headers, List<String> studyIds, QueryOptions options) {
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

    // TODO: this method is not going to work well because the sequence dictionarys are not stored correctly in the files collection
    private Set<String> getChromosomesFromVCFHeader(VCFHeader header, List<String> studyIds) {
        Set<String> chromosomes = new HashSet<>();
        // setup writers
        for (String studyId : studyIds) {
            SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
            chromosomes.addAll(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet()));
        }

        return chromosomes;
    }

    public String getOuputFilePath() {
        return outputFilePath.toString();
    }

    public int getFailedVariants() {
        return failedVariants;
    }
}
