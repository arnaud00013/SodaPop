#include "Population.h"

int Population::numberOfGenes = 0;
Input_Type Population::simType = Input_Type::selection_coefficient;
bool Population::noMut = false;

Population::Population():
    size_(0),
    sumFitness_(0),
	mutationCounter_(0),
	generation_(0)
{
	cells_.reserve(1000);
}

Population::Population(int targetSize):
    size_(0),
    sumFitness_(0),
	mutationCounter_(0),
	generation_(0)
{
    size_ = 0;
    cells_.reserve(targetSize);
}

Population::Population(std::ifstream& startFile,const std::string & genesPath, int targetSize, Init_Pop popType):
    size_(0),
    sumFitness_(0),
	mutationCounter_(0),
	generation_(0)
{
	if(popType==Init_Pop::from_snapFile){
		initPolyclonal(startFile, genesPath, targetSize);
	}
	else{
		initMonoclonal(startFile, genesPath, targetSize);
	}
    for (auto& cell : cells_) {
        addSumFitness(cell.fitness());
    }
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
    cells_ = std::vector <Cell>(targetSize , A);
    size_ = targetSize;
    for (auto& cell : cells_) {
        cell.ch_barcode(getBarcode());
    }
    if (Cell::ff_ == 5 || Cell::ff_ == 9){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
}

void Population::initPolyclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
	// ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
    std::cout << "Creating population from file ..." << std::endl;
    cells_.reserve(targetSize) ;
    int count = 0;
    while (count < Total_Cell_Count && !startFile.eof()){
        cells_.emplace_back(startFile, genesPath);
        ++count;  
    }
    if (Cell::ff_ == 5 || Cell::ff_ == 9){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
    size_ = count;
}

void Population::divide(int targetBuffer, int targetSize, std::ofstream& LOG, bool use_mut_log){
    // allocate space for temporary population
    Population newPopulation(targetBuffer);
    double relative_fitness(1);
    int n_progeny(0);
    for (auto& cell : cells_) {

        // fitness of cell j with respect to sum of population fitness
        relative_fitness = cell.fitness()/getSumFitness();

        // probability parameter of binomial distribution
        std::binomial_distribution<> binCell(targetSize, relative_fitness);

        // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
        n_progeny = binCell(g_rng);
            
        // if nil, the cell will be wiped from the population
        if(n_progeny == 0) continue; 

        // iterator to current available position
        auto it = std::end(newPopulation.cells_);

        // iterator to end position of fill
        auto last = it + n_progeny;

        //cell_it->setParent(cell_it - cells_.begin());
        
        newPopulation.fill_n(n_progeny,cell);

        auto link = it;

        do{
            link->linkGenes();
            ++link;
        }while(link < last);

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
                }while(it < last);
        }
    }

        // if the population is below N
        // randomly draw from progeny to pad
        while (newPopulation.getSize() < targetSize ){
            auto cell_it = newPopulation.cells_.begin();
            newPopulation.cells_.emplace_back(*(cell_it + randomNumber()*newPopulation.getSize()));
            newPopulation.incrementSize(1);
        }

        if (newPopulation.getSize() > targetSize ){
            std::shuffle(newPopulation.cells_.begin(), newPopulation.cells_.end(), g_rng);
            newPopulation.cells_.resize(targetSize);
            newPopulation.setSize(targetSize);
        }

        //alternative to shuffling
       /* while(v_size > N){
            int rand_idx = v_size*randomNumber();
            remove_at(newPopulation.cells_,rand_idx);
            v_size--;
        }*/

        Total_Cell_Count = newPopulation.getSize();
        assert (Total_Cell_Count == targetSize) ;
        
        // swap population with initial vector
        cells_.swap(newPopulation.cells_);

        // reset and update sumFitness_
        // update Ns and Na for each cell

        resetSumFitness();
        double fittest = 0;
        for (auto& cell : cells_) {
            double current = cell.fitness();
            addSumFitness(current);
            if (current > fittest) 
                fittest = current;
            cell.UpdateNsNa();
        }
        //normalize by fittest individual to prevent overflow
        if (Population::simType == Input_Type::selection_coefficient && (Cell::ff_ == 5)){
            for (auto& cell : cells_) {
                cell.normalizeFit(fittest);
            }
        }
}

void Population::fill_n(int n_progeny, const Cell& c)
{
    std::fill_n(std::back_inserter(cells_),n_progeny,c);
    incrementSize(n_progeny);
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
                for (const auto& cell : cells_) {
                    addSumFitness(cell.fitness());
                    cell.dump(toSnapshot,idx++);
                } 
            break;
        case Encoding_Type::no_sequence: //"short" output format
                for (const auto& cell : cells_) {
                    addSumFitness(cell.fitness());
                    cell.dumpShort(toSnapshot);
                }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
            break;
    }
}

void Population::simul_pev_before_cell_division(std::ofstream& pev_log,
		std::ofstream& cells_gene_content_log, std::ofstream& gain_log,
		std::ofstream& loss_log, int& nb_gain_to_sim, int& nb_loss_to_sim,
		int& gain_event_ctr, int& loss_event_ctr,
		int& total_nb_event_pev_to_sim, double& ratio_gain_current_gen,
		double& lambda_plus, double& lambda_minus, double& r_prime,
		double& a_for_s_x, double& b_for_s_x, bool& simul_pangenomes_evolution,
		bool& track_pangenomes_evolution, int& GENERATION_CTR, int& DT) {

	//Start
	double randnum = -1.0;//initialize a variable that will help us simulate gene gain and gene loss by respectio their ratio
	// for each cell
	for (auto& cell : cells_) {
		//Gene gain through HGT between and within species + Gene loss event(s)
		if (simul_pangenomes_evolution){
			nb_gain_to_sim = 0;
			nb_loss_to_sim = 0;
			// number of gain event
			nb_gain_to_sim = round((pow(cell.gene_count(),lambda_plus)));
			// number of loss event
			nb_loss_to_sim = round((r_prime*pow(cell.gene_count(),lambda_minus)));
			//Total number of pangenome evolution (gain and loss) events to simulate in the current generation
			total_nb_event_pev_to_sim = nb_gain_to_sim + nb_loss_to_sim;
			double nb_gain_for_ratio = nb_gain_to_sim;
			double nb_loss_for_ratio = nb_loss_to_sim;

			//Simulation of gene gain and gene loss events
			ratio_gain_current_gen = nb_gain_for_ratio / (nb_gain_for_ratio + nb_loss_for_ratio);
			for(int i=0;i<total_nb_event_pev_to_sim;i++){
				//draw a number between 0 and 1
				randnum = randomNumber();
				if (randnum <= ratio_gain_current_gen && (nb_gain_to_sim>0)){
					//choose random cell from which the gained gene will come (DON'T SHUFFLE cells_ here because you iterate through it)
					std::uniform_int_distribution<int> uniff_dist(0, cells_.size()-1);
					int random_index_cell = uniff_dist(g_rng);
					auto cell_two = cells_.begin() + random_index_cell;
					cell_two->select_random_gene();
					int ID_gene_gained = cell.add_gene(a_for_s_x,b_for_s_x);//this command updates automatically Cell gene array and fitness after gain!
					gain_event_ctr++;
					//If the user activated the option to get pangenome evolution feedbacks, save feedback for the gain event at each DT generations where DT is the time-step
					if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
						gain_log <<gain_event_ctr<<"\t"<<GENERATION_CTR<<"\t"<<cell_two->ID()<<"\t"<<cell_two->barcode()<<"\t"<<cell.ID()<<"\t"<<cell.barcode()<<"\t"<<ID_gene_gained<<std::endl;
					}
				}else{
					if (nb_loss_to_sim > 0){
						if (cell.gene_count()>1){
							int ID_gene_removed = cell.remove_rand_gene(a_for_s_x,b_for_s_x); ////this command updates automatically Cell gene array and fitness after loss!
							loss_event_ctr++;
							//If the user activated the option to get pangenome evolution feedbacks, save feedback for the loss event at each DT generations where DT is the time-step
							if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
								loss_log <<loss_event_ctr<<"\t"<<GENERATION_CTR<<"\t"<<cell.ID()<<"\t"<<cell.barcode()<<"\t"<<ID_gene_removed<<std::endl;
							}
						}
					}
				}
			}

			//If the user activated the option to get pangenome evolution feedbacks, save Feedback on genome size ( x ), loss/gain rate ratio ( r_x ), loss rate Beta_x and gain rate Alpha_x in PANGENOME_LOG at each DT generations where DT is the time-step. Do the same for cell gene content log with the appropriated fields.
			if (track_pangenomes_evolution && (((GENERATION_CTR % DT) == 0) || GENERATION_CTR==1)){
				pev_log <<GENERATION_CTR<<"\t"<<cell.ID()<<"\t"<<(cell.gene_count())<<"\t"<<(r_prime*pow(cell.gene_count(),(lambda_minus-lambda_plus)))<<"\t"<<(r_prime*pow(cell.gene_count(),lambda_minus))<<"\t"<<pow(cell.gene_count(),lambda_plus)<<"\t"<<cell.fitness()<<std::endl;
				cell.dumpCellGeneContent(cells_gene_content_log,GENERATION_CTR);
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
		// error opening file, throw exception
		throw std::runtime_error("Unable to open file for output.");
	}

}

double Population::addSumFitness(double w){
    sumFitness_ += w;
    return sumFitness_;
}
