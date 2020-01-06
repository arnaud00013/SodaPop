#ifndef CELL_H
#define CELL_H

#include "Gene.h"
#include "global.h"

class Gene;

/*SodaPop
Copyright (C) 2019 Louis Gauthier

    SodaPop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SodaPop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
 */

class Cell {

    typedef double(Cell::*funcPtr)(void) const;
public:

    static int ff_;
    static bool useDist_;
    static bool fromS_;
    static Gene selected_gene;
    static std::map<int, std::map<std::string, double>> map_codon_usage_maps_;//Key = species ID; Value =  Key-Value map of Codon usage frequency (will be converted to Relative Synonymous Codon Usage (RSCU) and normalized by maximum synonymous RSCU at initialization);

    Cell();
    Cell(std::ifstream &);
    Cell(std::ifstream & ,const std::string & );

    virtual ~Cell(){};

    int total_mutations(const int & );
    void FillGene_L();
    void linkGenes();

    void UpdateRates();

    int ID() const {return ID_;}
    uint32_t parent() const {return parent_;}
    double mrate() const {return c_mrate_;}
    int gene_count() const {return genomeVec_.size();}
    int genome_size() const {return geneBlocks_.back();}
    std::string barcode() const {return barcode_;}
    double fitness() const;
    double abs_fitness() const;
    double get_relative_cuf(const std::string &);

    void change_ID(int a) {ID_ = a;}
    void setParent(uint32_t a) {parent_ = a;}
    void ch_barcode(std::string s) {barcode_ = s;}

    void ch_Fitness(double f){fitness_ = f;}
    void ch_abs_Fitness(double f){abs_fitness_ = f;}

    // Fitness functions
    void selectFitness();
    double flux() const;
    double toxicity() const;
    double metabolicOutput() const;
    double multiplicative() const;
    double neutral() const;
    double noMut() const;
    double fold() const;
    double growthRate() const;
    double multiplicative_without_genes_fit_mean() const;

    void ranmut_Gene();
    void ranmut_Gene(std::ofstream&, int);
    void change_exprlevel();
    double normalizeFit(double);
    void dump(std::ofstream&, int) const;
    void dumpShort(std::ofstream&) const;
    void dumpParent(std::ofstream&) const;

    /**** HGT ****/
    bool any_mobile_gene_present() const;
    bool select_random_gene_gain();
    int add_gene(const double &);
    int remove_rand_gene(const double &);
    void print_summary_Gene_arr_();
    void print_summary_Gene_L_();

    double get_PevFe() const {return pev_fe_;}
    void set_PevFe(double PevFe) {pev_fe_ = PevFe;}

    double getSelCoeffCurrentMutation() const;
    void setSelCoeffCurrentMutation(double selCoeffCurrentMutation);
    void dumpCellGeneContent(std::ofstream&,int) const;
    /**** HGT ****/

    void PrintCell(int) const;

    int Na() const {return Total_Na_;}
    int Ns() const {return Total_Ns_;}
    void UpdateNsNa();

protected:
    // organism barcode
    std::string barcode_;

    // organism ID
    int ID_;

    // parent ID
    int parent_;

    // initial mutation rate
    double o_mrate_;

    // current mutation rate
    double c_mrate_;

    // organismal fitness (relative to fittest)
    double fitness_;

    // ABSOLUTE fitness
    double abs_fitness_;

    /**** HGT ****/
    //Pangenome evolution event (gain or loss of genes) fitness effect
    double pev_fe_;
    /**** HGT ****/

    //used to save the selection coefficient of a mutation event
    double sel_coeff_current_mutation_;

    //Vector of Mobile Genes index in genomeVec_; **NOTE THAT IF YOU WANT TO CONSIDER GENE ORDER / SYNTENY, each time genomeVec_ order would be changed, mobilomeVec_ would have to be updated in consequence!
    std::vector <int> mobilomeVec_;

    //Array of genes (mobile + non-mobile)
    std::vector <Gene> genomeVec_;
    
    
    //Cumulative sum of gene lengths (i.e. genome size)
    VectInt geneBlocks_;

    int Total_Ns_;
    int Total_Na_;

    // function pointer to select fitness function
    funcPtr fit;

};
#endif
