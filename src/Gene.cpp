// Gene.cpp
#include "Gene.h"

std::gamma_distribution<> Gene::gamma_ = std::gamma_distribution<>(1.0, 1.0);
std::normal_distribution<> Gene::normal_ = std::normal_distribution<>(1.0, 1.0);

//Input: gene file
Gene::Gene(std::fstream& gene_in,Cell *parent)
{
    myCell_ = parent;
    std::string line;
    // Read gene file line by line
    while (!gene_in.eof()){
        getline(gene_in,line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if (word == "Gene_NUM"){ 
            iss>>word; 
            g_num_ = atoi(word.c_str());
        }
        else if (word == "N_Seq"){ 
            iss>>nucseq_;
            ln_=nucseq_.length();
            if((ln_ % 3) != 0){
                  std::cerr << "Invalid length for nucleotide sequence: " << ln_ << std::endl;
                  exit(2);
            }
            else{
                  la_ = ln_/3;
                  std::string aaseq=GetProtFromNuc(nucseq_);

                  //check stop codons in midsequence
                  std::string::size_type loc = aaseq.find("X", 0 );
                  if(loc != std::string::npos){
                        std::cerr << "ERROR: DNA sequence has STOP codons in the middle"<<std::endl;
                        exit(2);
                    }           
            }
        }
        else if ( word=="E" ){
            iss >> word;
            e_ = atoi(word.c_str());
        }
        else if (word == "CONC"){
            iss >> word;
    	   conc_ = atof(word.c_str());
        }
        else if (word == "DG")
        { 
            iss>>word; 
      	    dg_ = atof(word.c_str());
            dg_ = exp(-dg_/kT);
        }
        else if (word == "F")
        { 
            iss>>word; 
            f_ = atof(word.c_str());
        }
        else if (word == "EFF")
        { 
            iss>>word; 
            eff_ = atof(word.c_str());
        }
        else if (word == "//"){;}//do nothing
    }
    Na_ = 0; //default
    Ns_ = 0;
    this->cumul_sum_fit_effect_mut_current_gen = 0;
}

// // copy constructor
Gene::Gene(const Gene& G)
{
    g_num_ = G.g_num_;
    ln_ = G.ln_;
    la_ = G.la_;
    nucseq_ = G.nucseq_;
    dg_ = G.dg_;
    f_ = G.f_;
    conc_ = G.conc_;
    e_ = G.e_;
    eff_ = G.eff_;
    Na_ = G.Na_;
    Ns_ = G.Ns_;
    this->cumul_sum_fit_effect_mut_current_gen = G.cumul_sum_fit_effect_mut_current_gen;
}

Gene::Gene()
{
	g_num_ = 0;
	ln_ = 0;
	la_ = 0;
	nucseq_ = "";
	dg_ = 0;
	f_ =5;
	conc_ = 0;
	e_ = 0;
	eff_ = 0;
	Na_ = 0;
	Ns_ = 0;
	this->cumul_sum_fit_effect_mut_current_gen = 0;
}

Gene::~Gene()
{
}

//Genes are equal if DNA sequence and concentration are equal.
bool Gene::operator== (Gene& G) 
{
    std::string temp = G.nseq();
    if ( (temp.compare(nucseq_) == 0) && (conc_ == G.conc_) ) return true;
    else return false;
}

// assignment overloading
Gene& Gene::operator=(const Gene& A)
{ 
    if (this != &A){
        this->g_num_ = A.g_num_;
        this->ln_ = A.ln_;
        this->la_ = A.la_;
        this->dg_ = A.dg_;
        this->f_ = A.f_;
        this->conc_ = A.conc_;
        this->e_ = A.e_;
        this->eff_ = A.eff_;
        this->Na_ = A.Na_;
        this->Ns_ = A.Ns_;
        (this->nucseq_).assign(A.nucseq_);
        this->cumul_sum_fit_effect_mut_current_gen = A.cumul_sum_fit_effect_mut_current_gen;
    }
    return *this;
}

/*
This version of the mutation function draws the DDG value from a gaussian distribution
distribution parameters are hardcoded for now
*/
double Gene::Mutate_Stabil_Gaussian(int i, int j)
{ 
    if(i>=ln_){
        std::cerr << "ERROR: Mutation site out of bounds. Gene "<<this->num()<<" length is : "<<this->length()<<" but selected site is "<<i<< std::endl;
        exit(2);
    }       

    //non-synonymous mutation
    if(randomNumber() <= fNs){

        double temp = Ran_Gaussian(1.0, 1.7);
        double x = exp(-temp/kT);

        dg_ *= x;
        Na_ += 1;

        return x;
    }
    else{
        Ns_ += 1;
        return 0;
    }
}

/*
This version of the mutation function gets the DDG value from the DDG matrix
input by the user.
INPUT: 
    i -> site to mutate
    j -> bp to mutate to
*/
std::string Gene::Mutate_Stabil(int i, int j)
{ 
    // extract codon to be mutated
    int cdn_ndx = (i%3);
    int cdn_start = i - cdn_ndx; 
    int resi = cdn_start/3;

    // fetch current codon
    std::string cdn_curr = nucseq_.substr(cdn_start, 3);
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string cdn_new = cdn_curr;

    std::string s = PrimordialAASeq.at(g_num_);
    // get amino acid from WT background
    int aa_primo = GetIndexFromAA(s.at(resi));

    // get mutated bp
    std::string bp = AdjacentBP(cdn_curr.substr(cdn_ndx, 1), j); //new BP
   
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);
    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    

    // get DDG value from matrix
    double x = matrix[g_num_][resi][aa_new-1];

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch primordial amino acid

    //Ignore mutations to and from CYSTEINE
    if( (aa_new==2) || (aa_curr==2)){
        return "CYSTEINE\tNA\tNA\tNA";
    }

    //Case unphysical DDG estimate
    if( x>DDG_min || x<DDG_max){
        return "UNPHYSICAL\tNA\tNA\tNA";
    }

    if( aa_curr == aa_new){//SILENT
          nucseq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else if(aa_primo == aa_new){//REVERT TO WT BACKGROUND

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 
          
          dg_ /= x_curr;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
    else{//TYPICAL NON-SYNONYMOUS

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 

          // assign new DG value
          // division accounts for mutation occuring on wildtype identity

          double diff = DDG_mean()-avg_DG;

          x *= exp(-diff/kT);

          //std::cout << -kT*log(dg_) << "\t" << -kT*log(x) << std::endl;

          dg_ /= x_curr;
          dg_ *= x;
          if(-kT*log(dg_) > 0){
            myCell_->ch_Fitness(0);
          }
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}

/*
This version of the mutation function draws the selection coefficient value from a gamma or normal distribution
*/
double Gene::Mutate_Select_Dist(int i, int j)
{ 
    if(i>=ln_){
    	std::cerr << "ERROR: Mutation site out of bounds. Gene "<<this->num()<<" length is : "<<this->length()<<" but selected site is "<<i<< std::endl;
        exit(2);
    }       
    // extract codon to be mutated
	int cdn_ndx = (i%3);
	int cdn_start = i - cdn_ndx;


	// fetch current codon
	std::string cdn_curr = nucseq_.substr(cdn_start, 3);
	std::string cdn_new = cdn_curr;
	// get mutated bp
	std::string bp = AdjacentBP( cdn_curr.substr(cdn_ndx, 1), j); //new BP
	// mutate codon
	cdn_new.replace(cdn_ndx, 1, bp);
	// check for stop codon
	cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);

	std::string the_gene_old_aa =this->getAAresidueFromCodonSequence(cdn_curr);
	std::string the_gene_new_aa =this->getAAresidueFromCodonSequence(cdn_new);

	if(the_gene_new_aa==the_gene_old_aa){ //Non-synonymous
		double s = RandomNormal();
		double wf = 1 + s;
		this->setCumulSumFitEffectMutCurrentGen(this->getCumulSumFitEffectMutCurrentGen()+s);
		f_ *= wf;
		nucseq_.replace(cdn_start, 3, cdn_new);
		Na_ += 1;
		return s;
	}else{ //Synonymous
		Ns_ += 1;
		return 0;
	}

}

/*
This version of the mutation function gets the selection coefficient value from the DMS matrix
input by the user.
INPUT: 
    i -> site to mutate
    j -> bp to mutate to
*/
std::string Gene::Mutate_Select(int i, int j)
{ 
    // extract codon to be mutated
    int cdn_ndx = (i%3);
    int cdn_start = i - cdn_ndx; 
    int resi = cdn_start/3;

    // fetch current codon
    std::string cdn_curr = nucseq_.substr(cdn_start, 3);
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string cdn_new = cdn_curr;

    std::string s = PrimordialAASeq.at(g_num_);     

    // get mutated bp
    std::string bp = AdjacentBP( cdn_curr.substr(cdn_ndx, 1), j); //new BP
   
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);
    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    
    // get selection coefficient from matrix
    double new_s = matrix[g_num_][resi][aa_new-1];

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch primordial amino acid

    if( aa_curr == aa_new){//SILENT
          nucseq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else{// NON-SYNONYMOUS 
          // assign new fitness value
          double new_f = f_ + new_s;
          this->setCumulSumFitEffectMutCurrentGen(this->getCumulSumFitEffectMutCurrentGen()+new_s);
          f_ = f_ * new_f;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}

void Gene::initGamma(double shape, double scale)
{
    Gene::gamma_.param(std::gamma_distribution<>::param_type(shape, scale));
}

void Gene::initNormal(double mean, double stddev)
{
    Gene::normal_.param(std::normal_distribution<>::param_type(mean, stddev));
}

double Gene::RandomGamma()
{
    return Gene::gamma_(g_rng);
}

double Gene::RandomNormal()
{
    return Gene::normal_(g_rng);
}

// Updates the current DNA sequence
void Gene::Update_Sequences(const std::string DNAsequence)
{ 
    int l = DNAsequence.length();

    if(l != ln_)
    {
        std::cerr << "ERROR: Replacing DNA sequence with a non-equal length DNA. "<< std::endl;
        std::cerr << "Make sure the gene list you provided matches the genes in the cell files."<< std::endl;
        exit(2);
    }       

    nucseq_ = DNAsequence;
}

double Gene::DDG_mean()
{
    return -0.3*-kT*log(dg_)-0.12;
}

// from Privalov 1979 (see also: Serohijos & Shakhnovich 2013)
// Boltzmann probability of the gene product to be in the native state
double Gene::Pnat()
{
    return dg_/(1+dg_);
}

// Number of functional copies in the cell
double Gene::functional()
{
    return conc_*Pnat();
}

// Number of misfolded copies in the cell
double Gene::misfolded()
{
    return conc_*(1-Pnat());
}

// Contribution to normalizing factor based on infinitely stable fold
double Gene::A_factor()
{
    return 1.0/conc_;
}

Cell *Gene::GetCell() const
{
    return myCell_;
}

const void Gene::setCell(Cell *C)
{
    myCell_ = C;
}

Gene::Gene(const Gene& G, Cell *p_new_Cell) {
	g_num_ = G.g_num_;
	ln_ = G.ln_;
	la_ = G.la_;
	nucseq_ = G.nucseq_;
	dg_ = G.dg_;
	f_ = G.f_;
	conc_ = G.conc_;
	e_ = G.e_;
	eff_ = G.eff_;
	Na_ = G.Na_;
	Ns_ = G.Ns_;
	this->cumul_sum_fit_effect_mut_current_gen =G.cumul_sum_fit_effect_mut_current_gen;
	myCell_= p_new_Cell;
}

double Gene::getCumulSumFitEffectMutCurrentGen() const {
	return cumul_sum_fit_effect_mut_current_gen;
}

void Gene::setCumulSumFitEffectMutCurrentGen(
		double cumulSumFitEffectMutCurrentGen) {
	cumul_sum_fit_effect_mut_current_gen = cumulSumFitEffectMutCurrentGen;
}

std::string Gene::getCodonSequence(const int i) {
	// extract codon to be mutated
	int cdn_ndx = (i%3);
	int cdn_start = i - cdn_ndx;

	// fetch current codon
	std::string cdn_curr = nucseq_.substr(cdn_start, 3);
	return cdn_curr;
}

std::string Gene::getAAresidueFromCodonSequence(const std::string& original_codon_seq) {
	std::string p_codon_seq = original_codon_seq;
	std::transform(p_codon_seq.begin(), p_codon_seq.end(), p_codon_seq.begin(), ::toupper);
	if (p_codon_seq == "GCT" or p_codon_seq == "GCC" or p_codon_seq == "GCA" or p_codon_seq == "GCG"){
		return "A";
	} else if (p_codon_seq == "TGT" or p_codon_seq == "TGC"){
		return "C";
	} else if (p_codon_seq == "GAT" or p_codon_seq == "GAC"){
		return "D";
	} else if (p_codon_seq == "GAA" or p_codon_seq == "GAG"){
		return "E";
	} else if (p_codon_seq == "TTT" or p_codon_seq == "TTC"){
		return "F";
	} else if (p_codon_seq == "GGT" or p_codon_seq == "GGC" or p_codon_seq == "GGA" or p_codon_seq == "GGG"){
		return "G";
	} else if (p_codon_seq == "CAT" or p_codon_seq == "CAC"){
		return "H";
	} else if (p_codon_seq == "ATT" or p_codon_seq == "ATC" or p_codon_seq == "ATA"){
		return "I";
	} else if (p_codon_seq == "AAA" or p_codon_seq == "AAG"){
		return "K";
	} else if (p_codon_seq == "TTA" or p_codon_seq == "TTG" or p_codon_seq == "CTT" or p_codon_seq == "CTC" or p_codon_seq == "CTA" or p_codon_seq == "CTG"){
		return "L";
	} else if (p_codon_seq == "ATG"){
		return "M";
	} else if (p_codon_seq == "AAT" or p_codon_seq == "AAC"){
		return "N";
	} else if (p_codon_seq == "CCT" or p_codon_seq == "CCC" or p_codon_seq == "CCA" or p_codon_seq == "CCG"){
		return "P";
	} else if (p_codon_seq == "CAA" or p_codon_seq == "CAG"){
		return "Q";
	} else if (p_codon_seq == "AGA" or p_codon_seq == "AGG" or p_codon_seq == "CGT" or p_codon_seq == "CGC" or p_codon_seq == "CGA" or p_codon_seq == "CGG"){
		return "R";
	} else if (p_codon_seq == "AGT" or p_codon_seq == "AGC" or p_codon_seq == "TCT" or p_codon_seq == "TCC" or p_codon_seq == "TCA" or p_codon_seq == "TCG"){
		return "S";
	} else if (p_codon_seq == "ACT" or p_codon_seq == "ACC" or p_codon_seq == "ACA" or p_codon_seq == "ACG"){
		return "T";
	} else if (p_codon_seq == "GTT" or p_codon_seq == "GTC" or p_codon_seq == "GTA" or p_codon_seq == "GTG"){
		return "V";
	} else if (p_codon_seq == "TGG"){
		return "W";
	} else if (p_codon_seq == "TAT" or p_codon_seq == "TAC"){
		return "Y";
	} else {
		return "STOP";
	}

}
