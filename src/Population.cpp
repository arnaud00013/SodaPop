#include "Population.h"

int Population::numberOfGenes = 0;
Input_Type Population::simType = Input_Type::selection_coefficient;
bool Population::noMut = false;
std::map<int,int> Population::reference_subpops_abundance_map={};

Population::Population():
    size_(0),
    sumFitness_(),
    mutationCounter_(0),
    generation_(0),
    subpopulations_abundance_()
{
    for (std::map<int,int>::iterator it_current_cell_type = reference_subpops_abundance_map.begin(); it_current_cell_type != reference_subpops_abundance_map.end();++it_current_cell_type){
        cells_[it_current_cell_type->first]={};
        cells_[it_current_cell_type->first].reserve(1000);
        resetSumFitness(it_current_cell_type->first);//initialize sumFitness for species
    }
}

Population::Population(int targetSize):
    size_(0),
    sumFitness_(),
    mutationCounter_(0),
    generation_(0),
    subpopulations_abundance_()
{
    size_ = 0;
    int sum_reserved = 0;
    for (std::map<int,int>::iterator it_current_cell_type = reference_subpops_abundance_map.begin(); it_current_cell_type != reference_subpops_abundance_map.end();++it_current_cell_type){//initialize cells map by reserving enough memory for species cells vector to avoid overflow and core_dumped_error
        cells_[it_current_cell_type->first]={};
        int current_species_pop_size = it_current_cell_type->second;
        int targetPopSize = current_species_pop_size < 10000 ? current_species_pop_size*5 : current_species_pop_size*2;
        cells_[it_current_cell_type->first].reserve(targetPopSize);
        sum_reserved += targetPopSize;
        resetSumFitness(it_current_cell_type->first);//initialize sumFitness for species
    }
    //make sure that the total memory reserved for each species cell vector is at least what is expected to avoid overflow
    if (sum_reserved < targetSize){
        std::cerr << "Logical Error! Risk of overflow"<<std::endl;
	exit(2);
    }

}


Population::Population(std::ifstream& startFile,const std::string & genesPath, int targetSize, Init_Pop popType):
    size_(0),
    sumFitness_(),
    mutationCounter_(0),
    generation_(0),
    subpopulations_abundance_()
{
    if(popType==Init_Pop::from_snapFile){
        initPolyclonal(startFile, genesPath, targetSize);
    }
    else{
        initMonoclonal(startFile, genesPath, targetSize);
    }
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        resetSumFitness(it_vec_cells_current_id->first);//initialize sumFitness for species
        for (auto& cell : it_vec_cells_current_id->second) {
            addSumFitness(it_vec_cells_current_id->first,cell.fitness());
        }
    }
    
    //update subpopulation_abundance_map based on current cell data
    this->update_subpop_abundance_map();
}

void Population::initLandscape(int fitArg, std::vector<std::string> matrixVec, std::string geneListFile, std::string genesPath){
    switch(Population::simType){
      case Input_Type::selection_coefficient:
          Cell::fromS_ = true;
          if (fitArg<5){
              // not allowed, should throw error
              std::cout << "Not a valid argument" << std::endl;
          }
          else Cell::ff_ = fitArg;

          InitMatrix();
          Population::numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);

          std::cout << "Gene count: " << numberOfGenes << std::endl;

          if(matrixVec.size()==1)
              ExtractDMSMatrix(matrixVec.front().c_str());
          else{
              Cell::useDist_ = true;
          }

        break;


      case Input_Type::stability:
          InitMatrix();
          Population::numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);
          Cell::ff_ = fitArg;
          // if DDG matrix is given
          if(matrixVec.size()){
              switch (matrixVec.size()){
                  case 2:
                      bind_DG = ExtractDDGMatrix(matrixVec.front().c_str(),Matrix_Type::is_binding);
                      std::cout << "-> Average ∆∆G_binding is " << bind_DG << " ..." << std::endl;
                  case 1:
                      fold_DG = ExtractDDGMatrix(matrixVec.front().c_str(),Matrix_Type::is_folding);
                      std::cout << "-> Average ∆∆G_folding is " << fold_DG << " ..." << std::endl;
                    break;
                  }
          }
          else{
              Cell::useDist_ = true;
          }

        break;
        
      default:
        break;
    }
}

void Population::initMonoclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
    std::cout << "Creating a population of " << targetSize   << " cells ..." << std::endl;
    //CMDLOG << "Creating a population of " << targetSize   << " cells ..." << std::endl;
    Cell A(startFile, genesPath);
    cells_[A.ID()] = std::vector <Cell>(targetSize , A);;
    size_ = targetSize;
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        for (auto& cell : it_vec_cells_current_id->second) {
            cell.ch_barcode(getBarcode());
            if (Cell::ff_ == 5 || Cell::ff_ == 9){
                cell.UpdateRates();
            }
            
        }
    }
}

void Population::initPolyclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
    // ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
    std::cout << "Creating population from file ..." << std::endl;
    int sum_reserved = 0;
    //reserve memory for each sell type vectors
    for (std::map<int,int>::iterator it_current_cell_type = reference_subpops_abundance_map.begin(); it_current_cell_type != reference_subpops_abundance_map.end();++it_current_cell_type){
        cells_.emplace(std::piecewise_construct, std::make_tuple(it_current_cell_type->first), std::make_tuple()); 
        int current_species_pop_size = it_current_cell_type->second;
        int targetPopSize = current_species_pop_size < 10000 ? current_species_pop_size*5 : current_species_pop_size*2;
        cells_[it_current_cell_type->first].reserve(targetPopSize);
        sum_reserved += targetPopSize;
    }
    //add cell in the appropriate cell vector of the map cells_
    int count = 0;
    while (count < Total_Cell_Count && !startFile.eof()){
        Cell current_cell(startFile, genesPath);
        cells_[current_cell.ID()].push_back(current_cell);
        ++count;  
    }
    
    if (Cell::ff_ == 5 || Cell::ff_ == 9){
        for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
            for (auto& cell : it_vec_cells_current_id->second) {
                cell.UpdateRates();
            }
        }    
            
    }
    
    size_ = count;
    //make sure that the total memory reserved for each species cell vector is at least what is expected to avoid overflow
    if (sum_reserved < targetSize){
        std::cerr << "Logical Error! Risk of overflow"<<std::endl;
	exit(2);
    }
}

void Population::divide(int targetBuffer, int targetSize, std::ofstream& LOG, bool use_mut_log,std::ofstream& pev_log, bool& track_pangenomes_evolution, double& lambda_plus, double& lambda_minus, double& r_prime, double& s_prime, int& GENERATION_CTR, int& DT){
    std::map<int,int> p_original_map_subpop_abundance(this->subpopulations_abundance_);
    // allocate space for temporary population
    Population newPopulation(targetBuffer);
    double relative_fitness(1);
    int n_progeny(0);
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        for (auto& cell : it_vec_cells_current_id->second) {
            // fitness of cell j with respect to sum of subpopulation fitness
            relative_fitness = cell.fitness()/getSumFitness(it_vec_cells_current_id->first);
            // probability parameter of binomial distribution
            std::binomial_distribution<> binCell(reference_subpops_abundance_map[it_vec_cells_current_id->first], relative_fitness);

            // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
            n_progeny = binCell(g_rng);
            
            // if nil, the cell will be wiped from the population
            if(n_progeny == 0) continue; 

            int old_size = newPopulation.cells_[cell.ID()].size();
            newPopulation.fill_n(n_progeny,cell);
            
            // iterator to current available position
            auto it = newPopulation.cells_[cell.ID()].begin() + old_size;//initialize to end of vector
            
            // iterator to end position of fill
            auto last = it + n_progeny;

            auto link = it;
            do{
                link->linkGenes();
                ++link;
            }while(link != last);
            if(!(Population::noMut)){
            // after filling with children, go through each one for mutation
                do{
                    std::binomial_distribution<> binMut(it->genome_size(), it->mrate());
                    int n_mutations = binMut(g_rng);
                    // attempt n mutations
                    for (int i=0;i<n_mutations;++i){
                        incrementMutationCount(1);
                        // change statement to switch
                        if (use_mut_log){
                            // mutate and write mutation to file
                            it->ranmut_Gene(LOG,getGeneration());
                        }
                        else{
                            it->ranmut_Gene();
                        }       
                    }
                    ++it;
                }while(it != last);
            }
        }
    }
    newPopulation.update_subpop_abundance_map();
    std::map<int, int>::iterator the_iter;
    for (the_iter=p_original_map_subpop_abundance.begin();the_iter!=p_original_map_subpop_abundance.end();++the_iter){
        while (newPopulation.get_abundance_by_cell_id(the_iter->first) < the_iter->second){
            //add random cell with cell_id
            auto cell_it = newPopulation.cells_[the_iter->first].begin(); //iterator of first element of species cells vector
            std::uniform_int_distribution<int> uniff_distr_cells_index(0, (newPopulation.cells_[the_iter->first]).size()-1);
            int random_cells_vec_index = uniff_distr_cells_index(g_rng);
            newPopulation.cells_[the_iter->first].push_back(*(cell_it + (random_cells_vec_index)));
            newPopulation.incrementSize(1);
            //update abundance map for cell_id
            newPopulation.increment_abundance_by_cell_id(the_iter->first);
        }
        if (newPopulation.get_abundance_by_cell_id(the_iter->first) > the_iter->second){
            //remove random cell with cell_id
            std::shuffle(newPopulation.cells_[the_iter->first].begin(), newPopulation.cells_[the_iter->first].end(), g_rng);
            newPopulation.cells_[the_iter->first].resize(the_iter->second);
            //update abundance map for cell_id and population size_ (included in the same function)
            newPopulation.update_subpop_abundance_map();
            newPopulation.setSize(targetSize);
        }
    }
    
    newPopulation.update_subpop_abundance_map();//update species abundance map and population size_ (included in the same function)

    Total_Cell_Count = newPopulation.getSize();
    assert (Total_Cell_Count == targetSize) ;
    
    // population species/cells map re-defined by newPopulation
    this->cells_=newPopulation.cells_;
    // update abundance map
    this->update_subpop_abundance_map();

    // reset and update sumFitness_
    // update Ns and Na for each cell
    std::map<int,double> fittest;
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        resetSumFitness(it_vec_cells_current_id->first); //reset sumFitness for species
        fittest[it_vec_cells_current_id->first] = 0;//initialize current species highest fitness to 0
        for (auto& cell : it_vec_cells_current_id->second) {
            double current = cell.abs_fitness();
            addSumFitness(it_vec_cells_current_id->first,current);
            if (current > fittest[it_vec_cells_current_id->first]){
                fittest[it_vec_cells_current_id->first] = current;
            }
            cell.UpdateNsNa();
        }
    }
    //normalize by fittest individual to prevent overflow
    if ((Population::simType == Input_Type::selection_coefficient && (Cell::ff_ == 5)) || (Population::simType == Input_Type::selection_coefficient && (Cell::ff_ == 9))){
        for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
            for (auto& cell : it_vec_cells_current_id->second) {
                cell.normalizeFit(fittest[cell.ID()]); //normalize fitness by subpopulation fittest
                //If the user activated the option to get pangenome evolution feedbacks, save Feedback on genome size ( x ), loss/gain rate ratio ( r_x ), loss rate Beta_x and gain rate Alpha_x in PANGENOME_LOG at each DT generations where DT is the time-step.
                if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
                    pev_log <<GENERATION_CTR<<"\t"<<cell.ID()<<"\t"<<(cell.gene_count())<<"\t"<<((1/s_prime)*r_prime*pow(cell.gene_count(),(lambda_minus-lambda_plus)))<<"\t"<<(r_prime*pow(cell.gene_count(),lambda_minus))<<"\t"<<s_prime*pow(cell.gene_count(),lambda_plus)<<"\t"<<cell.abs_fitness()<<std::endl;
                }
            }
        }
    }
}

void Population::fill_n(int n_progeny, const Cell& c)
{
    //avoid segmentation fault and add cells in the appropriate cell vector of the map cells_
    cells_[c.ID()].resize(cells_[c.ID()].size()+n_progeny, c);//resize vector
    incrementSize(n_progeny);
    update_subpop_abundance_map();
}

void Population::saveSnapshot(std::ofstream& toSnapshot, std::string dirName, int currentGen, Encoding_Type encoding){
    sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),dirName.c_str(), currentGen); 

    // Open snapshot file
    toSnapshot = std::ofstream(buffer, std::ios::out | std::ios::binary);
    if (toSnapshot.is_open()){
        //write
        writeSnapshotHeader(toSnapshot, encoding);
        writePop(toSnapshot, encoding);
        toSnapshot.close();
        // compress
        std::string command = "gzip -f ";
        command += buffer;
        const char *cmd = command.c_str();
        system(cmd);
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open file for output.");
    }
}

void Population::writeSnapshotHeader(std::ofstream& toSnapshot, Encoding_Type encoding)
{
    toSnapshot.write((char*)(&frame_time),sizeof(double));
    toSnapshot.write((char*)(&Total_Cell_Count),sizeof(int));
    toSnapshot.write((char*)(&encoding),sizeof(int));
}

void Population::writePop(std::ofstream& toSnapshot, Encoding_Type encoding){
    int idx = 1;
    switch (encoding){
        case Encoding_Type::by_default: //"normal" output format
        case Encoding_Type::full: 
            for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
                resetSumFitness(it_vec_cells_current_id->first); //initialize sumFitness for species
                for (auto& cell : it_vec_cells_current_id->second) {
                    addSumFitness(it_vec_cells_current_id->first,cell.fitness());
                    cell.dump(toSnapshot,idx++);
                } 
            }
            break;
        case Encoding_Type::no_sequence: //"short" output format
            for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
                resetSumFitness(it_vec_cells_current_id->first); //initialize sumFitness for species
                for (auto& cell : it_vec_cells_current_id->second) {
                    addSumFitness(it_vec_cells_current_id->first,cell.fitness());
                    cell.dumpShort(toSnapshot);
                }
            }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
            break;
    }
}

void Population::simul_pev_before_cell_division(std::ofstream& cells_gene_content_log, std::ofstream& gain_log,
        std::ofstream& loss_log, double& nb_gain_to_sim, double& nb_loss_to_sim,
        int& gain_event_ctr, int& loss_event_ctr,
        int& total_nb_event_pev_to_sim, double& ratio_gain_current_gen,
        double& lambda_plus, double& lambda_minus, double& r_prime, double& s_prime, double& b_for_s_x, bool& simul_pangenomes_evolution,
        bool& track_pangenomes_evolution, int& GENERATION_CTR, int& DT) {

    //Start
    std::vector<int> v_unique_species_ids = get_unique_species();
    double randnum = -1.0;//initialize a variable that will help us simulate gene gain and gene loss by respectio their ratio
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        for (auto& cell : it_vec_cells_current_id->second) {// for each cell
            //Gene gain through HGT between and within species + Gene loss event(s)
            if (simul_pangenomes_evolution){
                nb_gain_to_sim = 0;
                nb_loss_to_sim = 0;
                // number of gain event
                std::poisson_distribution<int> distribution_gain(((s_prime * pow(cell.gene_count(),lambda_plus))));
                nb_gain_to_sim = distribution_gain(g_rng);
                // number of loss event
                std::poisson_distribution<int> distribution_loss(((r_prime*pow(cell.gene_count(),lambda_minus))));
                nb_loss_to_sim = distribution_loss(g_rng);
                //Total number of pangenome evolution (gain and loss) events to simulate in the current generation
                total_nb_event_pev_to_sim = round(nb_gain_to_sim + nb_loss_to_sim);
                double nb_gain_for_ratio = nb_gain_to_sim;
                double nb_loss_for_ratio = nb_loss_to_sim;

                //Simulation of gene gain and gene loss events
                ratio_gain_current_gen = nb_gain_for_ratio / (nb_gain_for_ratio + nb_loss_for_ratio);
                for(int i=0;i<total_nb_event_pev_to_sim;i++){
                    //draw a number between 0 and 1
                    randnum = randomNumber();
                    if (randnum <= ratio_gain_current_gen && (nb_gain_to_sim>0)){
                        //choose random cell from which the gained gene will come (DON'T SHUFFLE cells_ here because you iterate through it)
                        std::uniform_int_distribution<int> uniff_distr_species_id(0, v_unique_species_ids.size()-1);
                        int random_index_species = uniff_distr_species_id(g_rng);
                        int the_random_species_id = *(v_unique_species_ids.begin()+random_index_species);
                        std::uniform_int_distribution<int> uniff_distr_cell_index(0, cells_[the_random_species_id].size()-1);
                        int random_index_cell = uniff_distr_cell_index(g_rng);
                        std::vector<Cell>::iterator cell_two = cells_[the_random_species_id].begin() + random_index_cell;
                        bool isAnyMobileGene_toGain = false;
                        isAnyMobileGene_toGain = cell_two->select_random_gene_gain();
                        if (!isAnyMobileGene_toGain){ //if there are no mobile genes in the genome, do not simulate any gain events
                            i = i - 1; // because the random cell selected did not have any mobile genes so the gain event was not simulated 
                            //std::cout << "Selected cell didn't have any mobile genes for HGT" << std::endl; //if this log is desired, uncomment this code line
                            continue;
                        }
                        int ID_gene_gained = cell.add_gene(b_for_s_x);//this command updates automatically Cell gene array and fitness after gain!
                        gain_event_ctr++;
                        //If the user activated the option to get pangenome evolution feedbacks, save feedback for the gain event at each DT generations where DT is the time-step
                        if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
                            gain_log <<gain_event_ctr<<"\t"<<GENERATION_CTR<<"\t"<<cell_two->ID()<<"\t"<<cell_two->barcode()<<"\t"<<cell.ID()<<"\t"<<cell.barcode()<<"\t"<<ID_gene_gained<<std::endl;
                        }
                    }else{
                        if (nb_loss_to_sim > 0){
                            if (cell.gene_count()>1){
                                int ID_gene_removed = cell.remove_rand_gene(b_for_s_x); ////this command updates automatically Cell gene array and fitness after loss!
                                loss_event_ctr++;
                                //If the user activated the option to get pangenome evolution feedbacks, save feedback for the loss event at each DT generations where DT is the time-step
                                if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
                                    loss_log <<loss_event_ctr<<"\t"<<GENERATION_CTR<<"\t"<<cell.ID()<<"\t"<<cell.barcode()<<"\t"<<ID_gene_removed<<std::endl;
                                }
                            }
                        }
                    }
                }

                //If the user activated the option to get pangenome evolution feedbacks, save cell gene content log with the appropriated fields (pangenome information are saved after fitness normalization in divide()).
                if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
                    cell.dumpCellGeneContent(cells_gene_content_log,GENERATION_CTR);
                }
            }
        }

    }
}

void Population::zip_PEV_Logs(std::string the_outputDir, std::string the_file, Encoding_Type encoding) {
    bool b_f_does_exit = false;
    sprintf(buffer,"%s/%s",the_outputDir.c_str(),the_file.c_str());


    if (FILE *file = fopen(buffer, "r")) {
        fclose(file);
        b_f_does_exit = true;
    }

    //verify if file exists
    if (b_f_does_exit){
        // compress
        std::string command = "gzip -f ";
        command += buffer;
        const char *cmd = command.c_str();
        system(cmd);

    }
    else{
        // error file does not exist -> throw exception
        std::string msg_err = "Error! The file you want to zip does not exist. The file you wanted to zip is: ";
        msg_err += buffer;
        throw std::runtime_error(msg_err);
    }

}

double Population::addSumFitness(const int& p_cell_id, double w){
    sumFitness_[p_cell_id] += w;
    return sumFitness_[p_cell_id];
}


void Population::update_subpop_abundance_map(){
    int tot_pop_size = 0;
    for (std::map<int,std::vector<Cell>>::iterator it_map_cells = cells_.begin();it_map_cells != cells_.end();++it_map_cells) {
        //update subpopulation_abundance_map based on current cell data
        subpopulations_abundance_[it_map_cells->first] = (it_map_cells->second).size();
        tot_pop_size += (it_map_cells->second).size();
    }
    this->setSize(tot_pop_size);
}   

bool Population::anyMobileGenesInPop(){
    for (std::map<int,std::vector<Cell>>::iterator it_vec_cells_current_id = cells_.begin();it_vec_cells_current_id!=cells_.end();++it_vec_cells_current_id){
        for (auto& cell : it_vec_cells_current_id->second) {
            if (cell.any_mobile_gene_present()){
                return true;
            }
        }
    }
    return false;
}

void Population::resetSumFitness(const int& p_cell_id) {
    sumFitness_[p_cell_id] = 0;
}

