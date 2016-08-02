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

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by pagarcia on 27/06/2016.
 */
public class BiodataVariantToVariantContextConverter {

    public static final String GENOTYPE_KEY = "GT";
    private final VariantContextBuilder variantContextBuilder;
    private List<VariantSource> sources;
    private CellbaseWSClient cellbaseClient;
    private String regionSequence;
    private Map<String, Map<String,String>> filesSampleNamesEquivalences;
    private static final int NO_CALL_ALLELE_INDEX = 2;

    public BiodataVariantToVariantContextConverter(CellbaseWSClient cellbaseWSClient) {
        this(null, cellbaseWSClient, null);
    }

    public BiodataVariantToVariantContextConverter(List<VariantSource> sources, CellbaseWSClient cellbaseWSClient,
                                                   Map<String, Map<String,String>> filesSampleNamesEquivalences)
    {
        this.sources = sources;
        this.cellbaseClient = cellbaseWSClient;
        this.filesSampleNamesEquivalences = filesSampleNamesEquivalences;
        variantContextBuilder = new VariantContextBuilder();
    }

    public VariantContext transform(Variant variant, Region region) throws CellbaseSequenceDownloadError {
        String[] allelesArray = getAllelesArray(variant, region);

        Set<Genotype> genotypes = getGenotypes(variant, allelesArray);

        VariantContext variantContext = variantContextBuilder
                .chr(variant.getChromosome())
                .start(variant.getStart())
                .stop(variant.getEnd())
                .noID()
                .alleles(allelesArray)
                .unfiltered()
                .genotypes(genotypes).make();
        return variantContext;
    }

    private String[] getAllelesArray(Variant variant, Region region) throws CellbaseSequenceDownloadError {
        String[] allelesArray;
        // if there are indels, we cannot use the normalized alleles, (hts forbids empty alleles) so we have to take them from cellbase
        if (variant.getReference().isEmpty() || variant.getAlternate().isEmpty()) {
            String contextNucleotide = getContextNucleotideFromCellbase(variant, variant.getStart(), region);
            // update variant ref, alt, start and end adding the context nucleotide
            variant.setReference(contextNucleotide + variant.getReference());
            variant.setAlternate(contextNucleotide + variant.getAlternate());
            variant.setStart(variant.getStart() - 1);
            variant.setEnd(variant.getStart() + variant.getReference().length() - 1);
        }
        allelesArray = new String[] {variant.getReference(), variant.getAlternate()};

        return allelesArray;
    }

    private String getContextNucleotideFromCellbase(Variant variant, Integer start, Region region) throws CellbaseSequenceDownloadError {
        if (cellbaseClient != null) {
            String contextNucleotide;
            int contextNucleotidePosition = start - 1;
            if (region != null) {
                contextNucleotide = getContextNucleotideFromCellbaseCachingRegions(variant.getChromosome(), contextNucleotidePosition, region);
            } else {
                try {
                    contextNucleotide = cellbaseClient.getSequence(new Region(variant.getChromosome(), contextNucleotidePosition, contextNucleotidePosition));
                } catch (Exception e) {
                    throw new CellbaseSequenceDownloadError("Error getting from Cellbase sequence for Region " + variant.getChromosome() + ":" +
                            contextNucleotidePosition + "-" + contextNucleotidePosition, e);
                }
            }
            return contextNucleotide;
        } else {
            throw new IllegalStateException(String.format(
                    "CellBase was not provided, needed to fill empty alleles in variant %s:%d:%s>%s", variant.getChromosome(),
                    variant.getStart(), variant.getReference(), variant.getAlternate()));
        }
    }

    private String getContextNucleotideFromCellbaseCachingRegions(String chromosome, int contextNucleotidePosition, Region region) throws CellbaseSequenceDownloadError {
        int regionStart = region.getStart() - 1;
        if (regionSequence == null) {
            // if an indel start is the first nucleotide of the region, we will need the previous nucleotide, so we are adding
            // the preceding nucleotide to the region (region.getStart()-1)
            int regionEnd = region.getEnd();
            try {
                regionSequence = cellbaseClient.getSequence(new Region(chromosome, regionStart, regionEnd));
            } catch (Exception e) {
                throw new CellbaseSequenceDownloadError("Error getting from Cellbase sequence for Region " + chromosome + ":"
                        + regionStart + "-" + regionEnd, e);
            }
        }
        String nucleotide = getNucleotideFromRegionSequence(contextNucleotidePosition, regionStart, regionSequence);
        return nucleotide;
    }

    private String getNucleotideFromRegionSequence(int nucleotidePosition, int regionStart, String regionSequence) {
        int relativePosition = nucleotidePosition - regionStart;
        return regionSequence.substring(relativePosition, relativePosition + 1);
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
        }
        return genotypes;
    }

    private Set<Genotype> getStudyGenotypes(Set<Genotype> genotypes, Allele[] variantAlleles, VariantSourceEntry variantStudyEntry) {
        for (Map.Entry<String, Map<String, String>> sampleEntry : variantStudyEntry.getSamplesData().entrySet()) {
            String sampleGenotypeString = sampleEntry.getValue().get(GENOTYPE_KEY);
            Genotype sampleGenotype =
                    parseSampleGenotype(variantAlleles, variantStudyEntry.getFileId(), sampleEntry.getKey(), sampleGenotypeString);
            genotypes.add(sampleGenotype);
        }
        return genotypes;
    }

    private Genotype parseSampleGenotype(Allele[] variantAlleles, String fileId, String sampleName, String sampleGenotypeString) {
        // use opencb biodata-models Genotype class for parsing the genotype string and get the list of genotype allele indexes
        org.opencb.biodata.models.feature.Genotype genotype =
                new org.opencb.biodata.models.feature.Genotype(sampleGenotypeString, variantAlleles[0].getBaseString(), variantAlleles[1].getBaseString());
        List<Allele> genotypeAlleles = new ArrayList<>(2);
        for (int index : genotype.getAllelesIdx()) {
            // every allele not 0 or 1 will be considered no call
            if (index == -1 || index > NO_CALL_ALLELE_INDEX) {
                index = NO_CALL_ALLELE_INDEX;
            }
            genotypeAlleles.add(variantAlleles[index]);
        }

        GenotypeBuilder builder = new GenotypeBuilder()
                .name(getFixedSampleName(fileId, sampleName))
                .alleles(genotypeAlleles);

        return builder.make();
    }

    private String getFixedSampleName(String fileId, String sampleName) {
        // this method returns the "studyId appended" sample name if there are sample name conflicts
        if (filesSampleNamesEquivalences != null) {
            return filesSampleNamesEquivalences.get(fileId).get(sampleName);
        } else {
            return sampleName;
        }
    }

    public void cleanCachedRegionSequence() {
        regionSequence = null;
    }

    public void setSources(List<VariantSource> sources) {
        this.sources = sources;
    }

    public void setFilesSampleNamesEquivalences(Map<String, Map<String, String>> filesSampleNamesEquivalences) {
        this.filesSampleNamesEquivalences = filesSampleNamesEquivalences;
    }
}
