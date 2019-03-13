// Gene.cpp
#include "Cell.h"
Gene Cell::selected_gene = Gene();

Cell::Cell():
    barcode_(getBarcode()),
    ID_(0),
    parent_(0),
	accum_pev_fe(0),
    o_mrate_(0),
    c_mrate_(0),
    fitness_(0)
    {
        Gene_L_.reserve(GENECOUNTMAX);
        Gene_arr_.reserve(GENECOUNTMAX);
    }

// Construct from cell file
Cell::Cell(std::fstream & cell_in) {
	accum_pev_fe = 0;
    char buffer[140];
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);
    ch_barcode(getBarcode());
    setParent(0);
    std::string line;
    // default genesPath, can be changed by user in command-line
    std::string genesPath = "files/genes/";
    while (!cell_in.eof()) {
        getline(cell_in, line);
        std::string word;
        std::istringstream iss(line, std::istringstream:: in );
        iss >> word;
        if (word == "genes_path") {
            iss >> word;
            genesPath = word.c_str();
        }
        if (word == "org_id") {
            iss >> word;
            ID_ = atoi(word.c_str());
        } else if (word == "mrate") {
            iss >> word;
            o_mrate_ = atof(word.c_str());
            c_mrate_ = atof(word.c_str());
        }
        else if (word == "fitness") {
            iss >> word;
            fitness_ = atof(word.c_str());
        } else if (word == "G") {
            //reading gene files; 
            //concentration and stability from gene file take precedence
            iss >> word;
            //open gene file
            sprintf(buffer, "%s%s.gene", genesPath.c_str(), word.c_str());
            std::fstream gene_data(buffer);
            if (!gene_data.is_open()) {
                std::cerr << "File could not be open: " << buffer << std::endl;
                exit(2);
            }
            Gene A(gene_data,this);
            Gene_arr_.push_back(A);

            //Check if gene is correctly inserted
            std::vector <Gene>::iterator i = Gene_arr_.end();
            i--;
            // if gene fitness is null, assign a randomly fit value
            if((*i).f()==0){
                (*i).ch_f(0.95+randomNumber()*0.02);
            }
            gene_data.close();
        }
    }
}

// Constructs from a unit cell stored in binary 
Cell::Cell(std::fstream & IN,
    const std::string & genesPath) {
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);

    char buffer[140];
    int cell_id, cell_index, gene_size;
    double m, f;

    IN.read((char*)(&cell_index), sizeof(int));
    IN.read((char*)(&cell_id), sizeof(int));
    ID_ = cell_id;
    parent_ = 0;
    accum_pev_fe = 0;

    //read barcode
    int l;
    IN.read((char*) &l, sizeof(int));
    //construct vector container with nl elements
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    ch_barcode(std::string().assign(buf.data(), buf.size()));

    IN.read((char*)(&f), sizeof(double));
    fitness_ = f;
    IN.read((char*)(&m), sizeof(double));
    c_mrate_ = m;
    o_mrate_ = m;
    IN.read((char*)(&gene_size), sizeof(int));
    //read gene info
    for (int j = 0; j < gene_size; j++) {
        double e, c, dg, f, eff;
        int gene_nid, Ns, Na;
        std::string DNAsequence;

        IN.read((char*)(&gene_nid), sizeof(int));
        IN.read((char*)(&e), sizeof(double));
        IN.read((char*)(&c), sizeof(double));
        IN.read((char*)(&eff), sizeof(double));
        IN.read((char*)(&dg), sizeof(double));
        IN.read((char*)(&f), sizeof(double));
        IN.read((char*)(&Ns), sizeof(int));
        IN.read((char*)(&Na), sizeof(int));

        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        //construct vector container with nl elements
        std::vector<char> buf(nl);
        IN.read(& buf[0], nl);
        DNAsequence.assign(buf.data(), buf.size());

        sprintf(buffer, "%s%d.gene", genesPath.c_str(), gene_nid);
        std::fstream gene_data(buffer);
        if (!gene_data.is_open()) {
            std::cerr << "ERROR: Cannot open gene file " << buffer << std::endl;
            exit(2);
        }
        Gene G(gene_data,this);
        //update gene information
        G.ch_conc(c);
        dg = exp(-dg / kT);
        G.ch_dg(dg);
        G.ch_eff(eff);
        G.ch_f(f);
        G.Update_Sequences(DNAsequence);
        G.ch_Na(Na);
        G.ch_Ns(Ns);
        Gene_arr_.push_back(G);
    }
}

void Cell::linkGenes()
{
    for(auto gene_it = this->Gene_arr_.begin(); gene_it != this->Gene_arr_.end(); gene_it++) {
        gene_it->setCell(this);
    }
}

//Select a random gene from the cell to be transferred to another cell by saving a copy of it in the static member selected_gene of the Gene class
void Cell::select_random_gene() {
	// 1. Create a temporary vector of indices corresponding to the actual gene objects
	std::vector<int> indices(Gene_arr_.size());
	std::iota(indices.begin(), indices.end(), 0);
	// 2. Shuffle the vector of indices using the already instantiated rng
	std::shuffle(indices.begin(), indices.end(), g_rng);
	// 3. Take the last element as ID of the gene to be selected
	int ID_random_gene = indices.back();
	// 4. Return the random gene
	Cell::selected_gene =  Gene(*(this->Gene_arr_.begin() + ID_random_gene),this);
}

int Cell::remove_rand_gene() {
	// 1. Create a temporary vector of indices corresponding to the actual gene objects
	std::vector<int> indices(Gene_arr_.size());
	std::iota(indices.begin(), indices.end(), 0);
	// 2. Shuffle the vector of indices using the already instantiated rng
	std::shuffle(indices.begin(), indices.end(), g_rng);
	// 3. Take the last element as ID of the gene to be removed
	int indice_removed_gene = indices.back();
	int ID_removed_gene = (Gene_arr_.begin() + indice_removed_gene)->num();
	// 4. Erase the gene from the vector. Note that the vector is automatically resized
	Gene_arr_.erase(Gene_arr_.begin() + indice_removed_gene);
	//Update Gene_L_
	this->Gene_L_.clear();
	this->Gene_L_.reserve(this->gene_count()-1);
	this->FillGene_L();

	return ID_removed_gene;
}

//Add the selected gene saved in the static memeber selected_gene in the present cell
int Cell::add_gene() {
	Cell::selected_gene.setCell(this);
	Gene_arr_.push_back(Cell::selected_gene);
	//std::cout<<"Gain event : Cell"<<this->ID()<<" new gene number is "<<n_G.num()<<" and has length : "<<n_G.length()<<std::endl;
	//Update Gene_L_
	this->Gene_L_.clear();
	this->Gene_L_.reserve(this->gene_count()+1);
	this->FillGene_L();

	return Cell::selected_gene.num();
}

void Cell::print_summary_Gene_arr_() {
	for (auto current_gene_it = this->Gene_arr_.begin(); current_gene_it!=this->Gene_arr_.end();++current_gene_it){
		 std::cout<<"current gene of cell"<< this->ID()<<" is gene"<<current_gene_it->num()<<" and has length "<<current_gene_it->length()<<std::endl;
	}
}

void Cell::print_summary_Gene_L_() {
	std::cout<<"cumulative lengths of cell"<< this->ID()<<" is :"<<std::endl;
	for (auto it_cumul_length = this->Gene_L_.begin(); it_cumul_length!=this->Gene_L_.end();++it_cumul_length){
		std::cout<<*(it_cumul_length)<<std::endl;
	}

}

double Cell::getCellCumulSumFitEffectMutCurrentGen() const {
	return Cell_cumul_sum_fit_effect_mut_current_gen;
}

void Cell::setCellCumulSumFitEffectMutCurrentGen(
		double cellCumulSumFitEffectMutCurrentGen) {
	Cell_cumul_sum_fit_effect_mut_current_gen =
			cellCumulSumFitEffectMutCurrentGen;
}

// // copy constructor
// Cell::Cell(const Cell& C)
// {
//     barcode_ = C.barcode_;
//     ID_ = C.ID_;
//     parent_ = C.parent_;
//     o_mrate_ = C.o_mrate_;
//     c_mrate_ = C.c_mrate_;
//     fitness_ = C.fitness_;
//     Gene_L_ = C.Gene_L_;
//     std::vector < Gene > ::iterator i;
//     for (auto cell_it = C.Gene_arr_.begin(); cell_it != C.Gene_arr_.end(); cell_it++) {
//         Gene A(*cell_it,this);
//         Gene_arr_.push_back(A);
//     }
// }

const double Cell::fitness()
{
    return fitness_;
}

// Initialize the cummulative gene length array
void Cell::FillGene_L() {
    int sum = 0;
    std::vector < Gene > ::iterator i;
    for (auto cell_it = Gene_arr_.begin(); cell_it != Gene_arr_.end(); ++cell_it) {
        sum += cell_it->length();
        Gene_L_.push_back(sum);
    }
}

// Return total mutation count
// spec:
//  - 0, Ns+Na
//  - 1, Ns
//  - 2. Na
int Cell::total_mutations(const int & spec) {
    assert((spec < 3) && (spec >= 0));

    int sa = 0;
    int s = 0;
    int a = 0;

    for (auto cell_it = Gene_arr_.begin(); cell_it != Gene_arr_.end(); ++cell_it) {
        int Ns = cell_it -> Ns();
        int Na = cell_it -> Na();
        s += Ns;
        a += Na;
        sa += (Ns + Na);
    }

    if (spec == 0) return sa;
    else if (spec == 1) return s;
    else return a;
}

//set to 0 the sum of fitness effect of the Cell and all its genes
void Cell::initializeCellCumulSumFitEffectMutCurrentGen() {
	this->setCellCumulSumFitEffectMutCurrentGen(0);
	for(auto gene_it = this->Gene_arr_.begin(); gene_it != this->Gene_arr_.end(); gene_it++) {
		gene_it->setCumulSumFitEffectMutCurrentGen(0);
	}
}
//Calculate the sum of the Cell genes mutations fitness effect
void Cell::updateCellCumulSumFitEffectMutCurrentGen() {
	for(auto gene_it = this->Gene_arr_.begin(); gene_it != this->Gene_arr_.end(); gene_it++) {
		this->setCellCumulSumFitEffectMutCurrentGen(this->getCellCumulSumFitEffectMutCurrentGen()+gene_it->getCumulSumFitEffectMutCurrentGen());
	}
}


