#ifndef POPULATION_H
#define POPULATION_H

#include "Cell.h"

class Population {

public:
	Population();
	Population(int);
	Population(std::ifstream & ,const std::string &, int, Init_Pop);

	void initMonoclonal(std::ifstream & ,const std::string &, int);
	void initPolyclonal(std::ifstream & ,const std::string &, int);

	static void initLandscape(int, std::vector<std::string>,std::string,std::string);

	void saveSnapshot(std::ofstream&, std::string, int, Encoding_Type);
	void writeSnapshotHeader(std::ofstream&, Encoding_Type);
	void writePop(std::ofstream&, Encoding_Type);

	void zip_PEV_Logs(std::string, std::string, Encoding_Type);
	void simul_pev_before_cell_division(std::ofstream& cells_gene_content_log,std::ofstream& gain_log,std::ofstream& loss_log, double& nb_gain_to_sim, double& nb_loss_to_sim, int& gain_event_ctr, int& loss_event_ctr,int& total_nb_event_pev_to_sim, double& ratio_gain_current_gen, double& lambda_plus, double& lambda_minus, double& r_prime, double& s_prime, double& a_for_s_x, double& b_for_s_x, bool& simul_pangenomes_evolution, bool& track_pangenomes_evolution, int& GENERATION_CTR, int& DT);
	void divide(int, int, std::ofstream&, bool,std::ofstream& pev_log, bool& track_pangenomes_evolution, double& lambda_plus, double& lambda_minus, double& r_prime, double& s_prime, int& GENERATION_CTR, int& DT);

	int getSize() const {return size_;}
 
        //std::unordered_set<int> Population::getUniqueSpeciesId() const;

	int getMutationCount() const {return mutationCounter_;}

	double getSumFitness() const {return sumFitness_;}

	int getGeneration() const {return generation_;}

	void incrementGeneration() {generation_++;}

	std::vector<Cell> getCells() const {return cells_;}

	void fill_n(int, const Cell&);

	void incrementMutationCount(int c) {mutationCounter_ += c;}

	void incrementSize(int c) {size_ += c;}

	void setSize(int s) {size_ = s;}

	void resetSumFitness() {sumFitness_ = 0;}

	double addSumFitness(double);

	static int numberOfGenes;

	static Input_Type simType;

	static bool noMut;

protected:
	int size_;
	double sumFitness_;
	int mutationCounter_;
	int generation_;
	std::vector<Cell> cells_;

};

#endif
