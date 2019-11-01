#include "Gene.h"

std::gamma_distribution<> Gene::gamma_ = std::gamma_distribution<>(1.0, 1.0);
std::normal_distribution<> Gene::normal_ = std::normal_distribution<>(1.0, 1.0);
std::exponential_distribution<> Gene::exponential_ = std::exponential_distribution<>(1.0);

//Input: gene file
Gene::Gene(std::ifstream& gene_in, Cell *parent)
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
            gene_idx_ = atoi(word.c_str());
        }
        else if (word == "N_Seq"){ 
            iss>>gene_seq_;
            gene_len_ = gene_seq_.length();
            if ((gene_len_ % 3) != 0){
                  std::cerr << "Invalid length for nucleotide sequence: " << gene_len_ << std::endl;
                  exit(2);
            }
            else{
                  prot_len_ = gene_len_/3;
                  std::string aaseq=GetProtFromNuc(gene_seq_);

                  //check stop codons in midsequence
                  std::string::size_type loc = aaseq.find('X', 0 );
                  if (loc != std::string::npos){
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
        else if (word == "DG"){
            iss>>word; 
      	    dg_ = atof(word.c_str());
            dg_ = exp(-dg_/kT);
        }
        else if (word == "F"){
            iss>>word; 
            f_ = atof(word.c_str());
        }
        else if (word == "EFF"){
            iss>>word; 
            eff_ = atof(word.c_str());
        }
        else if (word == "//")
            {;}//do nothing
    }
    Na_ = 0; //default
    Ns_ = 0;
    this->s_current_mutation = 0;
}

// copy constructor
Gene::Gene(const Gene& G)
{
    gene_idx_ = G.gene_idx_;
    gene_len_ = G.gene_len_;
    prot_len_ = G.prot_len_;
    gene_seq_ = G.gene_seq_;
    dg_ = G.dg_;
    f_ = G.f_;
    conc_ = G.conc_;
    e_ = G.e_;
    eff_ = G.eff_;
    Na_ = G.Na_;
    Ns_ = G.Ns_;
    this->s_current_mutation = G.s_current_mutation;
}

Gene::Gene()
{
    gene_idx_ = 0;
    gene_len_ = 0;
    prot_len_ = 0;
    gene_seq_ = "";
    dg_ = 0;
    f_ = 1;
    conc_ = 0;
    e_ = 0;
    eff_ = 0;
    Na_ = 0;
    Ns_ = 0;
    this->s_current_mutation = 0;
}

Gene::~Gene()
{
}

//Genes are equal if DNA sequence and concentration are equal.
bool Gene::operator== (Gene& G) 
{
    std::string temp = G.geneSeq();
    if ( (temp.compare(gene_seq_) == 0) && (conc_ == G.conc_) )
        return true;
    else
        return false;
}

// assignment overloading
Gene& Gene::operator=(const Gene& A)
{ 
    if (this != &A){
        this->gene_idx_ = A.gene_idx_;
        this->gene_len_ = A.gene_len_;
        this->prot_len_ = A.prot_len_;
        this->dg_ = A.dg_;
        this->f_ = A.f_;
        this->conc_ = A.conc_;
        this->e_ = A.e_;
        this->eff_ = A.eff_;
        this->Na_ = A.Na_;
        this->Ns_ = A.Ns_;
        (this->gene_seq_).assign(A.gene_seq_);
        this->s_current_mutation = A.s_current_mutation;
    }
    return *this;
}

/*
This version of the mutation function draws the DDG value from a gaussian distribution
distribution parameters are hardcoded for now
*/
double Gene::Mutate_Stabil_Gaussian(int i, int j)
{ 
    if (i>=gene_len_){
        std::cerr << "ERROR: Mutation site out of bounds."<< std::endl;
        exit(2);
    }       

    //non-synonymous mutation
    if (randomNumber() <= fNS){

        double temp = Ran_Gaussian(1.0, 1.7);
        double x = exp(-temp/kT);

        dg_ *= x;
        Na_ += 1;

        return x;
    }
    else{
        Ns_ += 1;
        return 1;
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

    // fetch current codon
    std::string cdn_curr = gene_seq_.substr(cdn_start, 3);

    // get mutated bp
    char bp = AdjacentBP(cdn_curr.at(cdn_ndx), j); //new BP
   
    // mutate codon
    std::string cdn_new = cdn_curr;
    cdn_new.replace(cdn_ndx, 1, 1, bp);

    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);

    // get amino acid from WT background
    int resi = cdn_start/3;

    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    
    // get DDG value from matrix
    double x = matrix[gene_idx_][resi][aa_new-1];

    std::string mutation = std::to_string(gene_idx_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string s = PrimordialAASeq.at(gene_idx_);
    int aa_primo = GetIndexFromAA(s.at(resi));

    //Ignore mutations to and from CYSTEINE
    if ( (aa_new==2) || (aa_curr==2)){
        return "CYSTEINE\tNA\tNA\tNA";
    }

    //Case unphysical DDG estimate
    if ( x>DDG_min || x<DDG_max){
        return "UNPHYSICAL\tNA\tNA\tNA";
    }

    if ( aa_curr == aa_new){//SILENT
          gene_seq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else if (aa_primo == aa_new){//REVERT TO WT BACKGROUND

          double x_curr = matrix[gene_idx_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 
          
          dg_ /= x_curr;
          gene_seq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
    else{//TYPICAL NON-SYNONYMOUS

          double x_curr = matrix[gene_idx_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 

          // assign new DG value
          // division accounts for mutation occuring on wildtype identity

          double diff = DDG_mean()-fold_DG;

          x *= exp(-diff/kT);

          //std::cout << -kT*log(dg_) << "\t" << -kT*log(x) << std::endl;

          dg_ /= x_curr;
          dg_ *= x;
          if (-kT*log(dg_) > 0){
            myCell_->ch_Fitness(0);
          }
          gene_seq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}

/*
This version of the mutation function draws the selection coefficient value from a gamma or normal distribution
*/
double Gene::Mutate_Select_Dist(int i, int j)
{ 
	if(i>=gene_len_){
		std::cerr << "ERROR: Mutation site out of bounds. Gene "<<this->num()<<" length is : "<<this->geneLength()<<" but selected site is "<<i<< std::endl;
		exit(2);
	}
	double s =0;//initialization of s
	// extract codon to be mutated
	int cdn_ndx = (i%3);
	int cdn_start = i - cdn_ndx;

	// fetch current codon
	std::string cdn_curr = gene_seq_.substr(cdn_start, 3);

	// get mutated bp
	char bp = AdjacentBP(cdn_curr.at(cdn_ndx), j); //new BP

	std::string cdn_new = cdn_curr;
	// mutate codon
	cdn_new.replace(cdn_ndx, 1, 1, bp);
	// check for stop codon
	cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);

	std::string the_gene_old_aa =this->getAAresidueFromCodonSequence(cdn_curr);
	std::string the_gene_new_aa =this->getAAresidueFromCodonSequence(cdn_new);

	if(the_gene_new_aa!=the_gene_old_aa){ //Non-synonymous
		if (randomNumber()<0.7){
			s = RandomNormal(); //little fitness effect (around neutral) represents 70% of the mutation fitness distribution
			Na_ += 1;
		}else{
			s = -1; //lethal mutations represents the rest of the mutation fitness distribution
			Na_ += 1;
		}

	}else{ //Synonymous
		s = 0;
		Ns_ += 1;
	}
	double wf = 1 + s;
	this->setS_current_mutation(s);
	f_ *= wf;
	gene_seq_.replace(cdn_start, 3, cdn_new);
	return s;
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

    std::string s = PrimordialAASeq.at(gene_idx_);

    // fetch current codon
    std::string cdn_curr = gene_seq_.substr(cdn_start, 3);

    // get mutated bp
    char bp = AdjacentBP(cdn_curr.at(cdn_ndx), j); //new BP
   
    std::string cdn_new = cdn_curr;
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);

    int resi = cdn_start/3;

    std::string mutation = std::to_string(gene_idx_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    // get selection coefficient from matrix
    double new_s = matrix[gene_idx_][resi][aa_new-1];
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);

    if ( aa_curr == aa_new){//SILENT
          new_s=0;
          gene_seq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          this->setS_current_mutation(new_s);
          return "SILENT\tNA\tNA\tNA";
    }
    else{// NON-SYNONYMOUS 
          // assign new fitness value
          double new_f = f_ + new_s;
          this->setS_current_mutation(new_s);
          f_ = f_ * new_f;
          gene_seq_.replace(cdn_start, 3, cdn_new);
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

void Gene::initExponential(double lambda)
{
    Gene::exponential_.param(std::exponential_distribution<>::param_type(lambda));
}

double Gene::RandomGamma()
{
    return Gene::gamma_(g_rng);
}

double Gene::RandomNormal()
{
    return Gene::normal_(g_rng);
}

double Gene::RandomExponential()
{
    return Gene::exponential_(g_rng);
}

// Updates the current DNA sequence
void Gene::Update_Sequences(const std::string DNAsequence)
{ 
    int l = DNAsequence.length();

    if (l != gene_len_){
        std::cerr << "ERROR: Replacing DNA sequence with a non-equal length DNA. "<< std::endl;
        std::cerr << "Make sure the gene list you provided matches the genes in the cell files."<< std::endl;
        exit(2);
    }       

    gene_seq_ = DNAsequence;
}

double Gene::DDG_mean() const
{
    return -0.3*-kT*log(dg_)-0.12;
}

// from Privalov 1979 (see also: Serohijos & Shakhnovich 2013)
// Boltzmann probability of the gene product to be in the native state
double Gene::Pnat() const
{
    return dg_/(1+dg_);
}

// Number of functional copies in the cell
double Gene::functional() const
{
    return conc_*Pnat();
}

// Number of misfolded copies in the cell
double Gene::misfolded() const
{
    return conc_*(1-Pnat());
}

// Contribution to normalizing factor based on infinitely stable fold
double Gene::A_factor() const
{
    return 1.0/conc_;
}

double Gene::getS_current_mutation() const {
    return s_current_mutation;
}

 void Gene::setS_current_mutation(
        double p_S_current_mutation) {
    s_current_mutation = p_S_current_mutation;
}

Gene::Gene(const Gene& G, Cell *p_new_Cell) {
    gene_idx_ = G.gene_idx_;
    gene_len_ = G.gene_len_;
    prot_len_ = G.prot_len_;
    gene_seq_ = G.gene_seq_;
    dg_ = G.dg_;
    f_ = G.f_;
    conc_ = G.conc_;
    e_ = G.e_;
    eff_ = G.eff_;
    Na_ = G.Na_;
    Ns_ = G.Ns_;
    this->s_current_mutation = G.s_current_mutation;
    myCell_= p_new_Cell;
}

Cell *Gene::GetCell() const
{
    return myCell_;
}

const void Gene::setCell(Cell *C)
{
    myCell_ = C;
}


std::string Gene::getCodonSequence(const int i) {
	// extract codon to be mutated
	int cdn_ndx = (i%3);
	int cdn_start = i - cdn_ndx;

	// fetch current codon
	std::string cdn_curr = gene_seq_.substr(cdn_start, 3);
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
