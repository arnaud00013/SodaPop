#include "PolyCell.h"

// By default the fitness function is set to neutral
int PolyCell::ff_ = 6;
bool PolyCell::useDist_ = false;
bool PolyCell::fromS_ = false;

PolyCell::PolyCell(){}
PolyCell::PolyCell(std::fstream& f) : Cell(f)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}    
PolyCell::PolyCell(std::fstream& f, const std::string& s) : Cell(f,s)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}

// Initialize the cumulative gene length array
void PolyCell::FillGene_L()
{
    int sum = 0;
    std::vector<Gene>::iterator i;
    for(i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        sum+= (i->length());
        Gene_L_.push_back(sum);
    }
}

void PolyCell::selectFitness()
{
    switch(PolyCell::ff_){
    	case 0: fit = &PolyCell::additive;
            break;
        case 1: fit = &PolyCell::fold;
            break;
        case 2: fit = &PolyCell::flux;
            break;
        case 3: fit = &PolyCell::toxicity;
            break;
        case 4: fit = &PolyCell::metabolicOutput;
            break;
        case 5: fit = &PolyCell::multiplicative;
            break;
        case 6: fit = &PolyCell::neutral;
            break;
        case 7: fit = &PolyCell::noMut;
            break;
        default:;
    }
}

// FOLDING-STABILITY BASED FLUX FITNESS FUNCTION
double PolyCell::fold()
{
    double sum_func = 0;
    //sum (concentration*Pnat) over all genes
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        sum_func += gene_it->functional();
    }
    return sum_func;
}

// METABOLIC FLUX FITNESS FUNCTION
double PolyCell::flux()
{
    double sum_func = 0;
    double a = 0;
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        // efficiency * Pnat * abundance
        sum_func += 1/(gene_it->eff()*gene_it->functional());
        a += gene_it->A_factor();
    }
    return static_cast<double>(a)/sum_func;
}

// MISFOLDING TOXICITY FITNESS FUNCTION
double PolyCell::toxicity()
{
    double f = 0;
    //sum (concentration*(1-Pnat)) over all genes
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        f += gene_it->misfolded();
    }
    return exp(-(COST*f));
}

// COMBINED METABOLIC OUTPUT FITNESS FUNCTION
double PolyCell::metabolicOutput()
{
    double flux = 0;
    double toxicity = 0;
    double a = 0;
    for (auto& it : Gene_arr_) {
        flux += 1 / (it.eff()*it.functional());
        toxicity += it.misfolded();
        a += it.A_factor();
    }

    flux = a/flux;
    toxicity = COST * toxicity;

    double fitness = flux - toxicity;
    // if toxicity > flux, return 0
    return (fitness < 0) ? 0 : fitness;
}

// MULTIPLICATIVE FITNESS FUNCTION
double PolyCell::multiplicative()
{
    double fitness = 0;
    for (auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it)
    {
        fitness += gene_it->f()*gene_it->e();
    }
    return (fitness/gene_count());
}

// NEUTRAL FITNESS FUNCTION
double PolyCell::neutral()
{
    return 1;
}

// NEUTRAL FITNESS FUNCTION
double PolyCell::noMut()
{
    return fitness();
}

void PolyCell::UpdateRates()
{
    fitness_ = (this->*fit)();
}

void PolyCell::ranmut_Gene(std::ofstream& log,int ctr)
{
    // get genome size
    int L = Gene_L_.back();

    // pick random site to mutate

    int site = (int) ( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    VectInt_citerator k = Gene_L_.begin();

    if(site >= (*k)){
    // random number generated is greater than
    // the first gene size
         for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
             if( site< (*k) ) break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    std::string mutation = "";

    int bp = (int) (3 * randomNumber());

    double wi = fitness();

    int the_gene_id = -1;//initialization with -1
	std::string the_gene_old_aa ="";
	std::string the_gene_new_aa ="";

    if(fromS_)
    {
        if(useDist_)
        {
            the_gene_id = j->num();
            the_gene_old_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        	(*j).Mutate_Select_Dist(site,bp);
        	the_gene_new_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        }
        else
        {
        	the_gene_id = j->num();
			the_gene_old_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        	mutation = j->Mutate_Select(site,bp);
        	the_gene_new_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        }
    }
    else
    {
        if(useDist_)
        {
        	the_gene_id = j->num();
        	the_gene_old_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
            (*j).Mutate_Stabil_Gaussian(site,bp);
            the_gene_new_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        }
        else
        {
        	the_gene_id = j->num();
			the_gene_old_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
            mutation = j->Mutate_Stabil(site,bp);
            the_gene_new_aa =j->getAAresidueFromCodonSequence(j->getCodonSequence(site));
        }
    }

    UpdateRates();
    double wf = fitness();
    double s = wf - wi;

    // save beneficial mutations to log
    // we could save all mutations with abs(s) >= some value x
    log << barcode().c_str() << "\t";
    log << std::fixed;
    log << mutation << "\t";
    log << s << "\t";
    log << ctr << std::endl;
    log<<this->barcode()<<"\t"<<the_gene_id<<"\t"<<the_gene_old_aa<<"\t"<<site<<"\t"<<the_gene_new_aa<<"\t"<<s<<"\t"<<ctr<<"\t"<<this->ID()<<std::endl;
}

void PolyCell::ranmut_Gene()
{
    // get genome size
    int L = Gene_L_.back();
    // pick random site to mutate

    int site = (int) ( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    VectInt_citerator k = Gene_L_.begin();


    if(site >= (*k)){
    // random number generated is greater than
    // the cumulative sum of genes
         for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
             if(site< (*k) ) break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    int bp = (int) (3 * randomNumber());
    // what is the input type?
    if(fromS_)
    {
        if(useDist_)
        {
            (*j).Mutate_Select_Dist(site,bp);
        }
        else
        {
            (*j).Mutate_Select(site,bp);
        }
    }
    else
    {
        if(useDist_)
        {
            (*j).Mutate_Stabil_Gaussian(site,bp);
        }
        else
        {
            (*j).Mutate_Stabil(site,bp);
        }
    }
         
    UpdateRates();
}

double PolyCell::normalizeFit(double fittest){
    if(fittest <= 0) {
        std::cerr << "Population collapse, average fitness is null.\n";
        exit(1);
    }
    double newfit = (Gene_arr_.begin()->f())/fittest;
    Gene_arr_.begin()->ch_f(newfit);
    UpdateRates();
    return newfit;
}

// Dump cell information to binary file
void PolyCell::dump(std::fstream& OUT, int cell_index)
{
    int x;
    double y;
    
    OUT.write((char*)(&cell_index),sizeof(int));
    
    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;		 	 
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());		 	 
    OUT.write((char*)(&x),sizeof(int));

   for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        int gene_nid = gene_it->num();
        double s = gene_it->e();
        double c = gene_it->conc();
        double eff = gene_it->eff();
        double dg = -kT*log(gene_it->dg());
        double f = gene_it->f();

        int Ns = gene_it->Ns();
        int Na = gene_it->Na();

        OUT.write((char*)(&gene_nid),sizeof(int));
        OUT.write((char*)(&s),sizeof(double));
        OUT.write((char*)(&c),sizeof(double));
        OUT.write((char*)(&eff),sizeof(double));
        OUT.write((char*)(&dg),sizeof(double));
        OUT.write((char*)(&f),sizeof(double));
        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene_it->nseq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell summary to binary file
void PolyCell::dumpShort(std::fstream& OUT)
{
    int x;
    double y;

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    x = Na();
    OUT.write((char*)(&x),sizeof(int));

    x = Ns();
    OUT.write((char*)(&x),sizeof(int));

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));
}

// Dump cell summary to binary file
void PolyCell::dumpSeq(std::fstream& OUT, int cell_index)
{
    int x;
    double y;

    OUT.write((char*)(&cell_index),sizeof(int));

    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));

    y = c_mrate_;            
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());             
    OUT.write((char*)(&x),sizeof(int));

    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){

        int Ns = gene_it->Ns();
        int Na = gene_it->Na();

        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene_it->nseq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell parent to binary file
void PolyCell::dumpParent(std::fstream& OUT)
{
    uint32_t a;
    a = parent();
    OUT.write((char*)(&a),sizeof(a));
}

double PolyCell::additive() {
	return this->fitness();

}

void PolyCell::dumpCellGeneContent(std::ofstream& OUT, int GEN_CTR) {
	for(auto gene_it = this->Gene_arr_.begin(); gene_it != this->Gene_arr_.end(); ++gene_it){
		OUT <<GEN_CTR<<"\t"<<this->ID()<<"\t"<<gene_it->num()<<std::endl;
	}
}

void PolyCell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        new_Na += gene_it->Na();
        new_Ns += gene_it->Ns();
    }
    Total_Na_ = new_Na;
    Total_Ns_ = new_Ns;
}

// Print cell information to stdout
void PolyCell::PrintCell(int cell_ndx)
{
      char buffer[140];
      sprintf(buffer,"C %6d %6d %12e %12e %d", cell_ndx, ID_, o_mrate_, mrate(), (int)Gene_arr_.size());  
      std::cout << buffer << std::endl;
      for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        int gene_nid = gene_it->num();
        double e = gene_it->e();
        double c = gene_it->conc();
        double dg = -kT*log(gene_it->dg());
        int Ns = gene_it->Ns();
        int Na = gene_it->Na();
           
        sprintf(buffer,"G %d% 2.2f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, Ns, Na);
        std::cout << buffer << std::endl;  
      }
      std::cout << std::endl;
}
