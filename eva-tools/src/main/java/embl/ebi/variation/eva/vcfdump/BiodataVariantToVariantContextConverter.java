/*
 *
 *  * Copyright 2016 EMBL - European Bioinformatics Institute
 *  *
 *  * Licensed under the Apache License, Version 2.0 (the "License");
 *  * you may not use this file except in compliance with the License.
 *  * You may obtain a copy of the License at
 *  *
 *  *      http://www.apache.org/licenses/LICENSE-2.0
 *  *
 *  * Unless required by applicable law or agreed to in writing, software
 *  * distributed under the License is distributed on an "AS IS" BASIS,
 *  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  * See the License for the specific language governing permissions and
 *  * limitations under the License.
 *
 */

package embl.ebi.variation.eva.vcfdump;

import embl.ebi.variation.eva.vcfdump.cellbasewsclient.CellbaseWSClient;
import embl.ebi.variation.eva.vcfdump.exception.CellbaseSequenceDownloadError;
import htsjdk.variant.variantcontext.*;
import org.opencb.biodata.models.feature.Region;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.biodata.models.variant.VariantSourceEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by pagarcia on 27/06/2016.
 */
public class BiodataVariantToVariantContextConverter {

    // TODO: remove if not used
    private static final Logger logger = LoggerFactory.getLogger(BiodataVariantToVariantContextConverter.class);


    private final VariantContextBuilder variantContextBuilder;
    private final List<VariantSource> sources;
    private CellbaseWSClient cellbaseClient;
    private String regionSequence;
    private Map<String, Map<String,String>> studiesSampleNamesEquivalences;

    // TODO: use a different value?
    public static final double UNKNOWN_QUAL = -10;

    // TODO: check if there is a PASS constant in any of the libraries
    public static final String FILTER_PASS = "PASS";

    public BiodataVariantToVariantContextConverter(List<VariantSource> sources, CellbaseWSClient cellbaseWSClient,
                                                   Map<String, Map<String,String>> studiesSampleNamesEquivalences)
    {
        this.sources = sources;
        this.cellbaseClient = cellbaseWSClient;
        this.studiesSampleNamesEquivalences = studiesSampleNamesEquivalences;
        variantContextBuilder = new VariantContextBuilder();
    }

    public VariantContext transform(Variant variant, Region region) throws CellbaseSequenceDownloadError {
        // TODO: alleles array: use array or another structure?
        String[] allelesArray = getAllelesArray(variant, region);

//        double qual = getQual(variant);

        Set<Genotype> genotypes = getGenotypes(variant, allelesArray);

        VariantContext variantContext = variantContextBuilder
                .chr(variant.getChromosome())
                // TODO: check start and end for indels
                .start(variant.getStart())
                .stop(variant.getEnd())
//                .id(String.join(";", variant.getIds()))   // in multiallelic, this results in duplicated ids, across several rows
                .noID()
                .alleles(allelesArray)
                // TODO: if qual is not set, which value is written to the VCF?
//                .log10PError(qual)
                // TODO: Info fields
                .filter(FILTER_PASS)
                .genotypes(genotypes).make();
        return variantContext;
    }

    private String[] getAllelesArray(Variant variant, Region region) throws CellbaseSequenceDownloadError {
        String[] allelesArray = {variant.getReference(), variant.getAlternate()};
        // if there are indels, we cannot use the normalized alleles, (hts forbids empty alleles) so we have to take them from cellbase
        boolean emptyAlleles = false;
        for (String a : allelesArray) {
            if (a.isEmpty()) {
                emptyAlleles = true;
                break;
            }
        }

        if (emptyAlleles) {
            String contextNucleotide;
            if (region != null) {
                contextNucleotide = getContextNucleotideFromCellbaseCachingRegions(variant, variant.getStart(), region);
            } else {
                contextNucleotide = getContextNucleotideFromCellbase(variant, variant.getStart());
            }
            // TODO: can be this made more efficient?
            allelesArray[0] = contextNucleotide + variant.getReference();
            allelesArray[1] = contextNucleotide + variant.getAlternate();
            variant.setEnd(variant.getStart() + allelesArray[0].length() - 1);
        }

        return allelesArray;
    }

    private String getContextNucleotideFromCellbase(Variant variant, Integer start) throws CellbaseSequenceDownloadError {
        if (cellbaseClient != null) {
            String contextNucleotide;
            int contextNucleotidePosition = start - 1;
            try {
                contextNucleotide = cellbaseClient.getSequence(new Region(variant.getChromosome(), contextNucleotidePosition, contextNucleotidePosition));
                return contextNucleotide;
            } catch (Exception e) {
                throw new CellbaseSequenceDownloadError("Error getting from Cellbase sequence for Region " + variant.getChromosome() + ":" +
                        contextNucleotidePosition + "-" + contextNucleotidePosition, e);
            }
        } else {
            throw new IllegalArgumentException(String.format(
                    "CellBase was not provided, needed to fill empty alleles in variant %s:%d:%s>%s", variant.getChromosome(),
                    variant.getStart(), variant.getReference(), variant.getAlternate()));
        }
    }

    private String getContextNucleotideFromCellbaseCachingRegions(Variant variant, int start, Region region) throws CellbaseSequenceDownloadError {
        if (cellbaseClient != null) {
            if (regionSequence == null) {
                // if an indel start is the first nucleotide of the region, we will need the previous nucleotide, so we are adding
                // the preceding nucleotide to the region (region.getStart()-1)
                int regionStart = region.getStart() - 1;
                int regionEnd = region.getEnd();
                try {
                    regionSequence = cellbaseClient.getSequence(new Region(variant.getChromosome(), regionStart, regionEnd));
                } catch (Exception e) {
                    throw new CellbaseSequenceDownloadError("Error getting from Cellbase sequence for Region " + variant.getChromosome() +
                            ":" + regionStart + "-" + regionEnd, e);
                }
            }
            String nucleotide = getNucleotideFromRegionSequence(start, region.getStart(), regionSequence);
            return nucleotide;
        } else {
            throw new IllegalArgumentException(String.format(
                    "CellBase was not provided, needed to fill empty alleles in variant %s:%d:%s>%s",
                    variant.getChromosome(), variant.getStart(), variant.getReference(), variant.getAlternate()));
        }
    }

    private String getNucleotideFromRegionSequence(int start, int regionStart, String regionSequence) {
        int relativePosition = start - regionStart;
        return regionSequence.substring(relativePosition, relativePosition + 1);
    }

    private double getQual(Variant variant) {
        double qual;
//        if (studies.size() == 1) {
//            // TODO: we need the file id
//            qual = Double.parseDouble(variant.getSourceEntries().get(studies.get(0)).getAttribute("QUAL"));
//        } else {
//            // TODO: there is a qual for study, how can we merge them?
//            qual = UNKNOWN_QUAL;
//        }
        // TODO: constant for unknown qual (.) en samtools library?
        qual = UNKNOWN_QUAL;
        return qual;
    }

    private Set<Genotype> getGenotypes(Variant variant, String[] allelesArray) {
        Set<Genotype> genotypes = new HashSet<>();

        Allele[] variantAlleles =
                {Allele.create(allelesArray[0], true), Allele.create(allelesArray[1]), Allele.create(Allele.NO_CALL, false)};

        for (VariantSource source : sources) {
            List<VariantSourceEntry> variantStudyEntries =
                    variant.getSourceEntries().values().stream().filter(s -> s.getStudyId().equals(source.getStudyId())).collect(Collectors.toList());
            for (VariantSourceEntry variantStudyEntry : variantStudyEntries) {
                genotypes = getStudyGenotypes(genotypes, variantAlleles, variantStudyEntry);
            }
            // TODO: fill with missing values (./.), or if we don't add the genotypes the writer will write the missing values? Check in VariantExporterController
//                GenotypeBuilder builder = new GenotypeBuilder().name(sampleName).alleles(null);

        }
        return genotypes;
    }

    private Set<Genotype> getStudyGenotypes(Set<Genotype> genotypes, Allele[] variantAlleles, VariantSourceEntry variantStudyEntry) {
        for (Map.Entry<String, Map<String, String>> sampleEntry : variantStudyEntry.getSamplesData().entrySet()) {
            // TODO: constant for GT?
            String gt = sampleEntry.getValue().get("GT");
            org.opencb.biodata.models.feature.Genotype genotype = new org.opencb.biodata.models.feature.Genotype(gt, variantAlleles[0].getBaseString(), variantAlleles[1].getBaseString());
            List<Allele> genotypeAlleles = new ArrayList<>();
            for (int index : genotype.getAllelesIdx()) {
                genotypeAlleles.add(variantAlleles[index]);
            }

            genotypes.add(new GenotypeBuilder().name(getFixedSampleName(variantStudyEntry.getStudyId(), sampleEntry.getKey())).alleles(genotypeAlleles).phased(genotype.isPhased()).make());
        }
        return genotypes;
    }

    private String getFixedSampleName(String studyId, String sampleName) {
        // this method returns the "studyId appended" sample name if there are sample name conflicts
        if (studiesSampleNamesEquivalences != null) {
            return studiesSampleNamesEquivalences.get(studyId).get(sampleName);
        } else {
            return sampleName;
        }
    }

}
