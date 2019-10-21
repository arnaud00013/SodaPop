#include "Population.h"
#include <tclap/CmdLine.h>
#include <unistd.h>
#include "rng.h"

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

int main(int argc, char *argv[])
{
    // these variables will hold the parameters input (or not) by the user
    int currentGen = 1;
    int maxGen = currentGen + 1;
    Encoding_Type outputEncoding = Encoding_Type::by_default;
    unsigned int targetPopSize = 1;
    int timeStep = 1;

    bool enableAnalysis = false;
    bool trackMutations = false;
    Init_Pop createPop = Init_Pop::from_snapFile;
    bool noMut = false;


    /*## HGT ##*/
    double expected_nb_gain_events = 0;
    double expected_nb_loss_events = 0;
    int ctr_nb_gain_events = 0;
    int ctr_nb_loss_events = 0;
    int total_nb_event_pev_to_sim = 0;
    int nb_cpus_for_variant_analysis = 0;
    double ratio_gain_current_gen = 0;
    double lambda_plus = 0;
    double lambda_minus = 0;
    double r_prime = 0;
    double s_prime = 0;
    double a_for_s_x = 0;
    double b_for_s_x = 0;

    bool simul_pangenomes_evolution = false;
    bool track_pangenomes_evolution = false;
    bool execute_variant_analysis = false;
    /*## HGT ##*/
    
    bool simul_codon_usage_bias_effect = false;    

    std::string geneListFile;
    std::string genesPath;
    std::string outDir;
    std::string startSnapFile;
    std::vector<std::string> matrixVec;

    // Wrap everything in a try block
    // errors in input are caught and explained to user
    try { 
        // Define the command line object
        TCLAP::CmdLine cmd("SodaPop: a multi-scale model of molecular evolution", ' ', "v1.0");

        // Define value arguments
        TCLAP::ValueArg<int> maxArg("m","maxgen","Number of generations",false,10,"int");
        TCLAP::ValueArg<int> popArg("n","size","Initial population size",false,1,"int");
        TCLAP::ValueArg<int> dtArg("t","dt","Time interval for snapshots",false,1,"int");

        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");
        TCLAP::ValueArg<std::string> geneArg("g","gene-list","Gene list file",true,"null","filename");
        TCLAP::ValueArg<std::string> startArg("p","pop-desc","Population description file",true,"null","filename");
        TCLAP::ValueArg<std::string> libArg("l","gene-lib","Gene library directory",false,"files/genes/","filename");

        TCLAP::MultiArg<std::string> matrixArg("i","input","Input file(s) defining the fitness landscape",false,"filepath(s)");
        TCLAP::ValueArg<int> fitArg("f","fitness","Fitness function",false,5,"integer ID");
        TCLAP::ValueArg<std::string> inputArg("","sim-type","Define simulation type\n<s> (from selection coefficient, DMS or otherwise)\n<stability> (from DDG matrix or distribution)", false,"s","string");
        TCLAP::SwitchArg gammaArg("","gamma","Draw selection coefficients from gamma distribution", cmd, false);
        TCLAP::SwitchArg normalArg("","normal","Draw selection coefficients from normal distribution", cmd, false);
        TCLAP::ValueArg<double> alphaArg("","alpha","Alpha parameter of distribution\nGamma -> shape\nNormal -> mean",false,1,"double");
        TCLAP::ValueArg<double> betaArg("","beta","Beta parameter of distribution\nGamma -> scale\nNormal -> S.D.",false,1,"double");
        TCLAP::SwitchArg initArg("c","create-single","Create initial population on the fly", cmd, false);
        TCLAP::SwitchArg analysisArg("a","analysis","Enable analysis scripts", cmd, false);
        
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);
        TCLAP::ValueArg<int> seqArg("s","seq-output","Sequence output format",false,0,"integer");
        TCLAP::ValueArg<unsigned long> seedArg("", "seed", "Seed value for RNG.", false, 0, "unsigned int (64-bit)");

        /*## HGT ##*/
        // boolean switch to simulate neutral Pangenomes evolution (Gain and loss of genes)
        TCLAP::SwitchArg pangenomes_evo_Arg("V","pangenomes-evolution","simulate neutral Pangenomes evolution (random Gain and loss of genes)", cmd, false);
        
        // boolean switch to track pangenomes evolution events (Gain and loss of genes) and also track the evolution of genome size and loss rate / gain rate ratio
        TCLAP::SwitchArg track_pangenomes_evo_Arg("T","track-PanEv","track pangenomes evolution events (Gain and loss of genes) and also track the evolution of genome size and loss rate / gain rate ratio", cmd, false);

        //parameter a to calculate the selection coefficient s(x)
        TCLAP::ValueArg<double> a_Arg("","aForSx","parameter a to calculate s(x)",false,1,"double");

        //parameter b to calculate the selection coefficient s(x)
        TCLAP::ValueArg<double> b_Arg("","bForSx","parameter b to calculate s(x)",false,1,"double");

        //parameter r_prime to calculate r(x)
        TCLAP::ValueArg<double> r_prime_Arg("","rPrime","parameter rPrime to calculate beta(x)",false,1,"double");

        //parameter s_prime to calculate r(x)
        TCLAP::ValueArg<double> s_prime_Arg("","sPrime","parameter sPrime to calculate alpha(x))",false,1,"double");

        //parameter lambda_plus to calculate alpha(x)
        TCLAP::ValueArg<double> lambda_plus_Arg("","lambdaPlus","parameter lambdaPlus to calculate alpha(x)",false,1,"double");

        //parameter lambda_minus to calculate beta(x)***
        TCLAP::ValueArg<double> lambda_minus_Arg("","lambdaMinus","parameter lambdaMinus to calculate beta(x)",false,1,"double");
        
        //Execute variant analysis
        TCLAP::SwitchArg exec_variant_analysis_Arg("","execVA","Execute variant analysis", cmd, false);

        //number of cpus for variant analysis (if applicable)
        TCLAP::ValueArg<int> nbCpuVA_Arg("u","nbCpu-VA","Number of cpus used for variant analysis",false,1,"int");
        /*## HGT ##*/

        /*## SodaPop argument that records if the user want to simulate the fitness effect of codon usage bias ##*/
        TCLAP::SwitchArg simul_codon_usage_bias_Arg("","simulCUB","Simulate the fitness effect of codon usage bias", cmd, false);

        //Standard deviation of the codon usage bias Distribution of Fitness Effect (dfe)
        TCLAP::ValueArg<double> standard_dev_cub_dfe_Arg("","stdCubDfe","Standard deviation of the codon usage bias Distribution of Fitness Effect (dfe)",false,0,"double");

        // Add the arguments to the CmdLine object.
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(dtArg);
        cmd.add(prefixArg);
        cmd.add(geneArg);
        cmd.add(startArg);
        cmd.add(libArg);
        cmd.add(fitArg);
        cmd.add(seqArg);
        cmd.add(matrixArg);
        cmd.add(alphaArg);
        cmd.add(betaArg);
        cmd.add(inputArg);
        cmd.add(standard_dev_cub_dfe_Arg);

        /*## HGT ##*/
        cmd.add(a_Arg);
        cmd.add(b_Arg);
        cmd.add(r_prime_Arg);
        cmd.add(s_prime_Arg);
        cmd.add(lambda_plus_Arg);
        cmd.add(lambda_minus_Arg);
        cmd.add(nbCpuVA_Arg);
        /*## HGT ##*/

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get values from args. 
        maxGen = maxArg.getValue();
        targetPopSize = popArg.getValue();
        timeStep = dtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        Population::simType = stringToInput_Type(inputArg.getValue());

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        std::cout << "Begin ... " << std::endl;

        if (matrixArg.isSet()){
            matrixVec = matrixArg.getValue();
        }

        // This statement should be a switch
        // make enum for simType
        if (gammaArg.isSet()){
            Gene::initGamma(alphaArg.getValue(), betaArg.getValue());
        }
        else if (normalArg.isSet()){
            Gene::initNormal(alphaArg.getValue(), betaArg.getValue());
        }
        else{
            // default to gamma
            Gene::initGamma(alphaArg.getValue(), betaArg.getValue());
        }

        //call init here
        Population::initLandscape(fitArg.getValue(), matrixVec,geneListFile,genesPath);

        enableAnalysis = analysisArg.getValue();
        trackMutations = eventsArg.getValue();
	
	// if user want to simulate HGT and loss, set HGT parameters to commandline input
        if (pangenomes_evo_Arg.isSet()){
            /*## HGT ##*/
            simul_pangenomes_evolution = pangenomes_evo_Arg.getValue();
            track_pangenomes_evolution = track_pangenomes_evo_Arg.getValue();
            lambda_plus = lambda_plus_Arg.getValue();
            lambda_minus = lambda_minus_Arg.getValue();
            r_prime = r_prime_Arg.getValue();
            s_prime = s_prime_Arg.getValue();
            a_for_s_x = a_Arg.getValue();
            b_for_s_x = b_Arg.getValue();
            if (exec_variant_analysis_Arg.isSet()){
                execute_variant_analysis = exec_variant_analysis_Arg.getValue();
                nb_cpus_for_variant_analysis = nbCpuVA_Arg.getValue();
            }
            /*## HGT ##*/
        }
        if (simul_codon_usage_bias_Arg.isSet()){
            simul_codon_usage_bias_effect = simul_codon_usage_bias_Arg.getValue();
            if (simul_codon_usage_bias_effect){
                if (!standard_dev_cub_dfe_Arg.isSet()){
                    std::cerr << "Error! You chosed to simulate Codon Usage Bias fitness effect but you didn't input any value for the standard deviation of the fitness effect distribution." << std::endl;
                    exit(2);
                }else{ 
                    Gene::stddev_cub_dfe_ = standard_dev_cub_dfe_Arg.getValue();
                    
                }
                
            }
        }


        outputEncoding = intToEncoding_Type(seqArg.getValue());
        createPop = intToPop_Type(initArg.getValue());

    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;}

    /* general simulation initialization, can be put in a global method
    */
    if (Cell::ff_ == 7) {
        noMut = true;
        Population::noMut = true;
        std::cout << "Mutations are not active." << std::endl;
    }

    createOutputDir(outDir);

    /* OPEN SIMULATION LOG
    */
    std::ofstream CMDLOG;
    try{
        openCommandLog(CMDLOG, outDir, argv, argc);
    }catch (std::runtime_error &e) {}

    /* OPEN MUTATION LOG
    */
    std::ofstream MUTATIONLOG;
    if (trackMutations && !noMut){
        try{
            openMutationLog(MUTATIONLOG,outDir);
        }catch (std::runtime_error &e) {}
    }

    // OPEN PANGENOME EVOLUTION EVENTS LOGS if the -V and -T command-line options are enabled
    std::ofstream PANGENOMES_EVOLUTION_LOG;
    std::ofstream CELL_GENE_CONTENT_LOG;
	std::ofstream GENE_GAIN_EVENTS_LOG;
	std::ofstream GENE_LOSS_EVENTS_LOG;
	if(simul_pangenomes_evolution && track_pangenomes_evolution){
		try{
			openPevLog(PANGENOMES_EVOLUTION_LOG, outDir);
			openCellsgenecontentLog(CELL_GENE_CONTENT_LOG, outDir);
			openGainLog(GENE_GAIN_EVENTS_LOG, outDir);
			openLossLog(GENE_LOSS_EVENTS_LOG, outDir);
		}catch (std::runtime_error &er) {}
	}


    /*OPEN POPULATION SNAPSHOT
    //READ POPULATION SNAPSHOT:
    */

    std::ifstream startFile;
    try{
        openStartingPop(startSnapFile,startFile);
        readSnapshotHeader(startFile);
    }catch (std::runtime_error &e) {}

    if (simul_codon_usage_bias_effect){ //if codon usage bias fitness effects are simulated create class Cell map of codon usage maps
        //get the vector of unique species ID
        std::vector<int> vec_unique_species_id;
        char path[] = "files/start";
        vec_unique_species_id = get_unique_species(path);
        std::cout << "-> Loading species codon usage files ..." <<std::endl;
        //Create std::map<int, std::map<std::string, double>> Cell::map_codon_usage_maps_
        for (auto const &current_sp_id : vec_unique_species_id){
            Cell::map_codon_usage_maps_[current_sp_id]=get_codon_usage_map(current_sp_id, "files/start");
            normalize_CUF(current_sp_id, Cell::map_codon_usage_maps_);
        }
        std::cout << "-> All species codon usage files were opened successfully ..." <<std::endl;
    }


    //std::vector <Cell> Cell_arr;
    Population currentPop(startFile, genesPath, targetPopSize, createPop);

    /* should be handled by method initializing the population
    */
    startFile.close();

    Total_Cell_Count = currentPop.getSize();

    std::ofstream OUT;

    try{
        currentPop.saveSnapshot(OUT,outDir,currentGen,outputEncoding);
    }catch (std::runtime_error &e) {}

    /* MAIN SIMULATION BLOCK, can be put in a method outside main
    */

    int targetBuffer = targetPopSize < 10000 ? targetPopSize*5 : targetPopSize*2;

    std::cout << "Starting evolution ..." << std::endl;
    CMDLOG << "Starting evolution ..." << std::endl;
    
    Gene::simulCub_ = simul_codon_usage_bias_effect; //send to the class Gene the boolean telling if codon usage bias fitness effects are simulated

    // // PSEUDO WRIGHT-FISHER PROCESS
    while (currentGen < maxGen){
        printProgress(currentGen*1.0/maxGen);

        //Pangenome evolution events
        //If the user chose to simulate pangenome evolution select automatically the "multiplicative_without_genes_fit_mean" fitness function
        if (simul_pangenomes_evolution){
        	Cell::ff_ = 9;
		currentPop.simul_pev_before_cell_division(CELL_GENE_CONTENT_LOG, GENE_GAIN_EVENTS_LOG, GENE_LOSS_EVENTS_LOG, expected_nb_gain_events, expected_nb_loss_events, ctr_nb_gain_events, ctr_nb_loss_events, total_nb_event_pev_to_sim, ratio_gain_current_gen, lambda_plus, lambda_minus, r_prime, s_prime, a_for_s_x, b_for_s_x, simul_pangenomes_evolution, track_pangenomes_evolution,currentGen,timeStep);
        }
        currentPop.divide(targetBuffer, targetPopSize, MUTATIONLOG,(trackMutations && !noMut),PANGENOMES_EVOLUTION_LOG,track_pangenomes_evolution,lambda_plus, lambda_minus, r_prime, s_prime,currentGen,timeStep);

        // update generation counter
        ++currentGen;
        // save population snapshot every timeStep generations
        if( (currentGen % timeStep) == 0){
            try{
                currentPop.saveSnapshot(OUT,outDir,currentGen,outputEncoding);
            }catch (std::runtime_error &e) {}
         }
     }

    printProgress(currentGen/maxGen);
    std::cout << std::endl;
    MUTATIONLOG.close();
    if(simul_pangenomes_evolution && track_pangenomes_evolution){
		PANGENOMES_EVOLUTION_LOG.close();
		CELL_GENE_CONTENT_LOG.close();
		GENE_GAIN_EVENTS_LOG.close();
		GENE_LOSS_EVENTS_LOG.close();
		std::string sim_out_dir = "out/"+outDir;
		currentPop.zip_PEV_Logs(sim_out_dir,"PANGENOMES_EVOLUTION_LOG.txt",outputEncoding);
		std::ofstream OUT_CELL_GENE_CONTENT_LOG;
		currentPop.zip_PEV_Logs(sim_out_dir,"CELL_GENE_CONTENT_LOG.txt",outputEncoding);
		std::ofstream OUT_GENE_GAIN_EVENTS_LOG;
		currentPop.zip_PEV_Logs(sim_out_dir,"GENE_GAIN_EVENTS_LOG.txt",outputEncoding);
		std::ofstream OUT_GENE_LOSS_EVENTS_LOG;
		currentPop.zip_PEV_Logs(sim_out_dir,"GENE_LOSS_EVENTS_LOG.txt",outputEncoding);
	}
    std::cout << "Done." << std::endl;
    CMDLOG << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << currentPop.getMutationCount() << std::endl;
    CMDLOG.close();

    // if the user toggled analysis, call shell script
    if (enableAnalysis){
        std::string script = "tools/barcodes.sh";
        std::string command = "/bin/bash "+script+" "+outDir+" "+std::to_string(maxGen)+" "+std::to_string(targetPopSize ) +" "+std::to_string(timeStep)+" "+std::to_string(outputEncoding)+" "+std::to_string(Population::numberOfGenes);
        const char *cmd = command.c_str();
        system(cmd);
    }

    
    //Execute the variant analysis in case that the user mentioned it and if mutations were active during simulation. 
    //WARNING : Rscript directory path needs to be in user system environment to execute this script!!!
    if (execute_variant_analysis){
        if (Cell::ff_ == 7){
            std::cerr << "Cannot execute variant analysis when mutations are not active!";
            exit(1);
        }
        else{
            std::cout << "Genes variant analysis started!" << std::endl;
            std::string script = "tools/Post_simulation_genes_variants_analysis.r";
            std::string sim_wp = geneListFile; //absolute path of Sodapop workspace initialized with gene list file path
            std::string the_subst_to_remove = "files/genes/gene_list.dat"; //substring of geneListFile to remove to get workspace
            std::string::size_type the_pos_substr_to_remove = sim_wp.find(the_subst_to_remove);

            if (the_pos_substr_to_remove != std::string::npos){
                sim_wp.erase(the_pos_substr_to_remove, the_subst_to_remove.length());
            }
            //keep a_for_s_x and b_for_s_x digits precision
            char buffer_a[32];
            memset(buffer_a, 0, sizeof(buffer_a));
            snprintf(buffer_a, sizeof(buffer_a), "%g", a_for_s_x);
            std::string str_a_for_s_x(buffer_a);
            char buffer_b[32];
            memset(buffer_b, 0, sizeof(buffer_b));
            snprintf(buffer_b, sizeof(buffer_b), "%g", b_for_s_x);
            std::string str_b_for_s_x(buffer_b);

            std::string command = "Rscript "+script+" "+sim_wp+" "+std::to_string(targetPopSize)+" "+std::to_string(timeStep) +" "+outDir+" "+std::to_string(r_prime)+" "+std::to_string(s_prime)+" "+std::to_string(lambda_plus)+" "+std::to_string(lambda_minus)+" "+std::to_string(Cell::ff_)+" "+str_a_for_s_x+" "+str_b_for_s_x+" "+std::to_string(nb_cpus_for_variant_analysis);
            std::cout << "The command is :" << std::endl;
            std::cout << command << std::endl;
            const char *cmd = command.c_str();
            system(cmd);
            
        }
    }

    return 0;
}
