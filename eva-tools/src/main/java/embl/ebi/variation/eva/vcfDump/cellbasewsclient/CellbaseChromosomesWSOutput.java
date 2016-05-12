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
package embl.ebi.variation.eva.vcfdump.cellbasewsclient;

import org.opencb.cellbase.core.common.core.Chromosome;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by pagarcia on 18/04/2016.
 */
public class CellbaseChromosomesWSOutput {
    private String species;
    private List<Chromosome> chromosomes;
    private List<Chromosome> supercontigs;

    public CellbaseChromosomesWSOutput() {
    }

    public CellbaseChromosomesWSOutput(String species, List<Chromosome> chromosomes, List<Chromosome> supercontigs) {
        this.species = species;
        this.chromosomes = chromosomes;
        this.supercontigs = supercontigs;
    }

    public String getSpecies() {
        return this.species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public List<Chromosome> getChromosomes() {
        return this.chromosomes;
    }

    public void setChromosomes(List<Chromosome> chromosomes) {
        this.chromosomes = chromosomes;
    }

    public List<Chromosome> getSupercontigs() {
        return supercontigs;
    }

    public void setSupercontigs(List<Chromosome> supercontigs) {
        this.supercontigs = supercontigs;
    }

    public Set<String> getAllChromosomeNames() {
        List<Chromosome> allChromosomes = new ArrayList<Chromosome>(chromosomes);
        allChromosomes.addAll(supercontigs);
        return allChromosomes.stream().map(c -> c.getName()).collect(Collectors.toSet());
    }
}
