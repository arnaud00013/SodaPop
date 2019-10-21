############################################################## Script Header ##################################################
#@Author = Arnaud NG
#This script analyses a single simulation results and give useful plots to understand genes evolution and population dynamics during simulations
#Output will be in folder "Sim_Results/" of Sodapop workspace directory. It includes plots and a matrix summarizing genes variants analysis called mtx_variants_analysis_results_simulation_(SimulationName)!
library("ggplot2")
library("nlme")
library("session")
library("doBy")
library("reconPlots")
#library("Biostrings")
#library("msa")
library("gridExtra")
library("parallel")
library("foreach")
library("doParallel")
library("lmPerm")


#SCRIPT ARGUMENTS :
sodapop_workspace <- as.character(commandArgs(TRUE)[1]) #ABSOLUTE Path of the folder containing SodaPop files (should be the parent of directory of SodaPop "out/" or "src/")
sodapop_workspace <- ifelse(test=substr(x = sodapop_workspace,start = nchar(sodapop_workspace),stop = nchar(sodapop_workspace))=="/",yes=sodapop_workspace,no=paste0(sodapop_workspace,"/"))
sim_output_workpace <- paste0(sodapop_workspace,"out/") #ABSOLUTE path of the folder containing Sodapop simulations raw output repertory 
results_output_workspace <- paste0(sodapop_workspace,"Sim_Results/") #ABSOLUTE path of the folder containing the Post-simulation genes variants analysis results. Will be created if it does not exist yet.

nb_cells_simulated <- as.integer(commandArgs(TRUE)[2]) #Number of cells simulated
mu <- 5E-10 #(mutations per bp) in the order of average prokaryotic mutation rate from Drake et al. 1991 
dt <- as.integer(commandArgs(TRUE)[3]) #Timestep used to determine the simulation reports frequency 
simulation_name <- as.character(commandArgs(TRUE)[4]) #The simulations output directories name 
rPrime <- as.numeric(commandArgs(TRUE)[5]) #See Sela, Wolf & Koonin (2016) parameter r'
sPrime <- as.numeric(commandArgs(TRUE)[6]) #See Sela, Wolf & Koonin (2016) parameter s'
lambdaPlus <- as.numeric(commandArgs(TRUE)[7]) #See Sela, Wolf & Koonin (2016) parameter lambda+
lambdaMinus <- as.numeric(commandArgs(TRUE)[8]) #See Sela, Wolf & Koonin (2016) parameter lambda-
lambda <-lambdaMinus - lambdaPlus #See Sela, Wolf & Koonin (2016) parameter lambda
num_fitness_landscape <- as.integer(commandArgs(TRUE)[9]) #Fitness landscape function number in Sodapop documentation
a_for_sx <- as.numeric(commandArgs(TRUE)[10]) #See Sela, Wolf & Koonin (2016) Equation 10
b_for_sx <- as.numeric(commandArgs(TRUE)[11]) #See Sela, Wolf & Koonin (2016) Equation 10
nb_cores_variant_analysis <- as.integer(commandArgs(TRUE)[12]) #Number of CPUs used for genes variant analysis

############################################################## Script Body ###################################################
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("mkdir -p ",results_output_workspace),intern = FALSE,wait = TRUE)

system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("cat ",sodapop_workspace,"files/start/*.cell | grep -Fhc G > ", results_output_workspace,"start_total_nb_genes.csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("ls ",sodapop_workspace,"files/start/*.cell | wc -l >> ", results_output_workspace,"start_total_nb_genes.csv"),intern = FALSE,wait = TRUE)

v_tot_and_mean <- read.csv(file = paste0(results_output_workspace,"start_total_nb_genes.csv"),header = FALSE,sep = "\t",stringsAsFactors = FALSE)$V1
average_nb_genes_per_genomes_at_start <- floor(v_tot_and_mean[1]/v_tot_and_mean[2])

system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("rm -f ", results_output_workspace,"start_total_nb_genes.csv"),intern = FALSE,wait = TRUE)

#function for plotting linear model
ggplotRegression <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  require(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 2)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 2)
  p_val <- formatC(summary(fit)$coefficients[,3][2], format = "e", digits = 3)
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    xlab(xlabl)+
    ylab(ylabl)+
    labs(title = paste("Adj R2 = ",adj_r_sq,
                       " Slope =",slope,
                       " P =",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
   
    if (bool_gg_save){
        ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 20, units = "cm")
    }else{
        print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
    }
    #return result as the real float numbers
    adj_r_sq <- summary(fit)$adj.r.squared
    slope <-summary(fit)$coefficients[,1][2]
    p_val <- summary(fit)$coefficients[,3][2]
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}

wanted_x_0_at_equilibrium = ((rPrime/sPrime)^(-1/lambda))

#create a function that returns number of possible SINGLE-SITE synonymous mutations divided by 3 for a CODON
calculate_third_of_possible_ns_codon <- function(the_codon){
  the_codon <- toupper(the_codon)
  if (nchar(the_codon)!=3){
    stop("codon length should be 3!")
  }
  possible_single_site_mutated_codons <- rep("",9)
  num_mut_codon <-1
  for (pos_codon in 1:3){
    if (substr(the_codon,start = pos_codon,stop=pos_codon)=="A"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
      
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="T"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="C"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else{#G
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }
  }
  #count the number of synonymous mutations based on the genetic code
  nb_unique_syn_mut_codons <-0 #default initialization
  if (the_codon == "TTT") {
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTC"])
    
  } else if (the_codon == "TTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTT"])
    
  } else if (the_codon == "TTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTG","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCC","TCA","TCG","AGT","AGC")])
  } else if (the_codon == "TCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCC","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCC","AGT","AGC")])
    
  } else if (the_codon == "TAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAC")])
    
  } else if (the_codon == "TAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAT")])
    
  } else if (the_codon == "TGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGC")])
    
  } else if (the_codon == "TGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGT")])
    
  } else if (the_codon == "TGG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "CTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTC","CTA","CTG")])
    
  } else if (the_codon == "CTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTA","CTG")])
    
  } else if (the_codon == "CTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTG")])
    
  } else if (the_codon == "CTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTA")])
    
  } else if (the_codon == "CCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCC","CCA","CCG")])
    
  } else if (the_codon == "CCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCA","CCG")])
    
    
  } else if (the_codon == "CCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCG")])
    
  } else if (the_codon == "CCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCA")])
    
  } else if (the_codon == "CAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAC")])
    
  } else if (the_codon == "CAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAT")])
    
  } else if (the_codon == "CAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAG")])
    
  } else if (the_codon == "CAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAA")])
    
  } else if (the_codon == "CGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGC","CGA","CGG")])
    
  } else if (the_codon == "CGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGG")])
    
  } else if (the_codon == "CGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGC","CGG")])
    
  } else if (the_codon == "CGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGC")])
    
  } else if (the_codon == "ATT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATA")])
    
  } else if (the_codon == "ATC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATT","ATA")])
    
  } else if (the_codon == "ATA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATT")])
    
  } else if (the_codon == "ATG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "ACT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACC","ACA","ACG")])
    
    
  } else if (the_codon == "ACC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACA","ACG")])
    
  } else if (the_codon == "ACA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACG")])
    
    
  } else if (the_codon == "ACG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACA")])
    
  } else if (the_codon == "AAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAC")])
    
  } else if (the_codon == "AAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAT")])
    
  } else if (the_codon == "AAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAG")])
    
  } else if (the_codon == "AAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAA")])
    
  } else if (the_codon == "AGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGC","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGT","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGG")])
    
  } else if (the_codon == "AGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGA")])
    
  } else if (the_codon == "GTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTG")])
    
  } else if (the_codon == "GTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTT","GTA","GTG")])
    
  } else if (the_codon == "GTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTT","GTG")])
    
  } else if (the_codon == "GTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTT")])
    
  } else if (the_codon == "GCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCG")])
    
  } else if (the_codon == "GCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCT","GCA","GCG")])
    
  } else if (the_codon == "GCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCT","GCG")])
    
  } else if (the_codon == "GCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCT")])
    
  } else if (the_codon == "GAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAC")])
    
  } else if (the_codon == "GAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAT")])
    
  } else if (the_codon == "GAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAG")])
    
  } else if (the_codon == "GAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAA")])
    
  } else if (the_codon == "GGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGG")])
    
  } else if (the_codon == "GGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGT","GGA","GGG")])
    
  } else if (the_codon == "GGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGT","GGG")])
    
    
  } else if (the_codon == "GGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGT")])
  }
  return((nb_unique_syn_mut_codons/3))
}


#create a function that does all possible pariwise COMPARISON (not alignment) between all the sequences in a vector(so if n is the number of sequences there will be n*(n-1) comparisons between reads pair codons. The algorithm execution time is in O(n^2) because of the if condition and the 2 for loop on the same vectorof n sequences.)
#The output is the total number of pairwise differences. 
calculate_nb_pwd <- function(vector_read_codon_seqs){
  output <- 0 #initialize the variable that will contain the number of pairwise differences in single nucleotide sites
  
  pwdiff_with_others<- function(the_codon,lst_unique_codons){
    if (any(duplicated(lst_unique_codons))){
      stop("There are doublons in the list of unique codons.  Retry!")
    }else if (!all(nchar(c(the_codon,lst_unique_codons))==3)){
      stop("codon length should be 3!")
    }
    result<-rep(NA,length(lst_unique_codons))
    ind_res_list<-1
    for (current_diff_codon in lst_unique_codons){
      current_result<-0
      for (i in 1:3){
        if (substr(x = the_codon,i,i)!=substr(x = current_diff_codon,i,i)){
          current_result <- current_result +1
        }
      }
      result[ind_res_list] <- current_result
      ind_res_list <- ind_res_list + 1
    }
    return(result)
  }
  v_unique_freqs <- list(table(vector_read_codon_seqs))[[1]]
  
  for (current_unique_codon in names(v_unique_freqs)){
    output <- output + sum(v_unique_freqs[current_unique_codon]*pwdiff_with_others(the_codon = current_unique_codon,lst_unique_codons = names(v_unique_freqs)[which(names(v_unique_freqs)==current_unique_codon):length(names(v_unique_freqs))]))
  }
  return(output)
}  

get_consensus_sequence_codons <- function(vector_codon_seqs){
  if (!all(nchar(vector_codon_seqs)==3)){
    stop("codon length should be 3!")
  }
  res_pos_1 <- names(which.max(table(substr(x = vector_codon_seqs,start = 1,stop = 1))))
  res_pos_2 <- names(which.max(table(substr(x = vector_codon_seqs,start = 2,stop = 2))))
  res_pos_3 <- names(which.max(table(substr(x = vector_codon_seqs,start = 3,stop = 3))))
  return(paste0(res_pos_1,res_pos_2,res_pos_3))
}

get_nb_seg_sites_from_codon_seqs <- function(vector_codon_seqs){
  if (!all(nchar(vector_codon_seqs)==3)){
    stop("codon length should be 3!")
  }
  bool_pos_1_seg <- length(table(substr(x = vector_codon_seqs,start = 1,stop = 1)))>1
  bool_pos_2_seg <- length(table(substr(x = vector_codon_seqs,start = 2,stop = 2)))>1
  bool_pos_3_seg <- length(table(substr(x = vector_codon_seqs,start = 3,stop = 3)))>1
  return(sum(as.integer(x = c(bool_pos_1_seg,bool_pos_2_seg,bool_pos_3_seg))))
}

get_nb_bp_diffs_consensus_codon_vs_orig <- function(orig_codon,consensus_seq){
  if (!all(nchar(c(consensus_seq,orig_codon))==3)){
    stop("codon length should be 3!")
  }
  
  return(sum(as.integer(substr(x = orig_codon,start =1 ,stop = 1)!=substr(x = consensus_seq,start =1 ,stop = 1))+as.integer(substr(x = orig_codon,start =2 ,stop = 2)!=substr(x = consensus_seq,start =2 ,stop = 2))+as.integer(substr(x = orig_codon,start =3 ,stop = 3)!=substr(x = consensus_seq,start =3 ,stop = 3))))
}

get_synonymity_stats <- function(the_df_cell_content_log, the_gene_ID,generation_of_interest, with_indels=FALSE){
  #Initializations
  Nb_syn_sites <- 0 #number of synonymous sites
  Nb_nsyn_sites <- 0 #number of non-synonymous sites
  Nb_syn_mutations <- 0 #number of synonymous mutations COMPARED TO THE REFERENCE GENE SEQUENCE AT TIME 0 (GENERATION 0)
  Nb_nsyn_mutations <- 0 #number of non-synonymous mutations COMPARED TO THE REFERENCE GENE SEQUENCE AT TIME 0 (GENERATION 0)
  nb_segreg_sites <- 0 #number of segregating sites
  Nb_pwdiff_gene <- 0 #initialize the variable that will represent the number of pairwise differences in the gene
  
  df_the_gene_in_cells_at_generation_of_interest <- subset(x = the_df_cell_content_log,subset = (gene_ID==the_gene_ID)&(Generation_ctr==generation_of_interest))
  #handle the case when there are no copy of the gene in the current generation
  if (nrow(df_the_gene_in_cells_at_generation_of_interest)==0){
    return(list(nb_nsm=NA,nb_sm=NA,nb_nss=NA,nb_ss=NA,nb_segregative_sites=NA,nb_copy_g_in_G=0,nb_uniq_species_of_g_in_G=0,a1=NA,Ne_k_hat_gene_in_current_g=NA,Ne_S_Taj_gene_in_current_g=NA,D_Taj_gene_in_current_g=NA))
  }
  gene_seq <- read.csv(file = paste0(sodapop_workspace,"files/genes/",the_gene_ID,".gene"),header = FALSE,sep = "\t",stringsAsFactors = FALSE)$V3[4]
  gene_length <- nchar(gene_seq) 
  
  if (with_indels){
    v_gene_ref_and_sim_copies <- c(gene_seq,df_the_gene_in_cells_at_generation_of_interest$nucl_sequence)
    v_gene_ref_and_sim_copies <- v_gene_ref_and_sim_copies[!is.na(v_gene_ref_and_sim_copies)]
    currentAlignment <- msa(inputSeqs = v_gene_ref_and_sim_copies,type = "dna")
    consensus_seq_gene_copies_from_sim <- toupper(msaConsensusSequence(currentAlignment, type=c("upperlower"), thresh=c((100-(1E-16)),(1E-16)), ignoreGaps=TRUE))
    myMaskedAlignment <- currentAlignment
    colM <- IRanges(start=1, end=gene_length)
    colmask(myMaskedAlignment) <- colM
    unmasked_alignment <- unmasked(myMaskedAlignment)
    
    #go through the gene sequence and calculate Nb_syn_sites, Nb_syn_mutations, Nb_nsyn_mutations
    for (pos_in_gene in seq(from = 1,to = gene_length,by = 3)){
      current_codon_gene <- substr(x = gene_seq,start = pos_in_gene,stop=pos_in_gene+2)
      Nb_syn_sites <- Nb_syn_sites + calculate_third_of_possible_ns_codon(current_codon_gene)
      original_codon <- current_codon_gene
      all_codons <- substr(x = as.character(unmasked_alignment[,1]),start = pos_in_gene,stop=pos_in_gene+2)
      all_codons <- all_codons[!grepl(pattern = "-",x = all_codons)]
      consensus_codon_seq <- substr(x = consensus_seq_gene_copies_from_sim,start = pos_in_gene,stop=pos_in_gene+2)
      mutated_codons <- unique(all_codons[all_codons!=original_codon])
      if (any(nchar(c(original_codon, mutated_codons))!=3)){
        print(paste("problem with codon length for gene ",df_site$contig_name, " at 0-based position ",df_site$pos, ", original codon is ", original_codon, "mutated codons are ", mutated_codons," and gene sequence is ",gene_seq))
        stop()
      }
      nb_sm_consensus <- 0 #initialization 
      nb_nsm_consensus <- 0 #initialization 
      #count the number of synonymous and non-synonymous mutations based on the genetic code
      if (original_codon == "TTT") {
        
        
        nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTC"))
        nb_nsm_consensus <- 1 - nb_sm_consensus
        
      } else if (original_codon == "TTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTG","CTT","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TTG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","CTT","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCC","TCA","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCA","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCC","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCA","TCC","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAA"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TAG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TGT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TGA"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TGG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- as.integer(consensus_codon_seq!=original_codon)
        
      } else if (original_codon == "CTT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTC","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTG"){
        nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTC","CTA"))
        nb_nsm_consensus <- 1 - nb_sm_consensus
        
      } else if (original_codon == "CCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCC","CCA","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCA","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCC","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCC","CCA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGC","CGA","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGA","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGC","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGA","CGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATC","ATA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATT","ATA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATC","ATT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- as.integer(consensus_codon_seq!=original_codon)
        
      } else if (original_codon == "ACT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACC","ACA","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACA","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACC","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACC","ACA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGC","TCT","TCC","TCA","TCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGT","TCT","TCC","TCA","TCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTA","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTT","GTA","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTT","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTA","GTT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCA","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCT","GCA","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCT","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCA","GCT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGA","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGT","GGA","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGT","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGG"){
        
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGA","GGT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      }
      
      nb_pwdiff_codon <- calculate_nb_pwd(all_codons)
      Nb_pwdiff_gene <- Nb_pwdiff_gene + nb_pwdiff_codon
      nb_seg_site_current_codon <- get_nb_seg_sites_from_codon_seqs(all_codons)
      nb_segreg_sites <- nb_segreg_sites + nb_seg_site_current_codon
      nb_bp_diffs_current_codon <- get_nb_bp_diffs_consensus_codon_vs_orig(orig_codon = original_codon,consensus_seq = consensus_codon_seq)
      Nb_syn_mutations <- Nb_syn_mutations + (nb_sm_consensus*nb_bp_diffs_current_codon)
      Nb_nsyn_mutations <- Nb_nsyn_mutations + (nb_nsm_consensus*nb_bp_diffs_current_codon)
      
    }#end for loop codons
    
  }else{
    #go through the gene sequence and calculate Nb_syn_sites, Nb_syn_mutations, Nb_nsyn_mutations
    for (pos_in_gene in seq(from = 1,to = gene_length,by = 3)){
      current_codon_gene <- substr(x = gene_seq,start = pos_in_gene,stop=pos_in_gene+2)
      Nb_syn_sites <- Nb_syn_sites + calculate_third_of_possible_ns_codon(current_codon_gene)
      original_codon <- current_codon_gene
      all_codons <- substr(x = df_the_gene_in_cells_at_generation_of_interest$nucl_sequence,start = pos_in_gene,stop=pos_in_gene+2)
      consensus_codon_seq <- get_consensus_sequence_codons(all_codons)
      mutated_codons <- unique(all_codons[all_codons!=original_codon])
      if (any(nchar(c(original_codon, mutated_codons))!=3)){
        print(paste("problem with codon length for gene ",df_site$contig_name, " at 0-based position ",df_site$pos, ", original codon is ", original_codon, "mutated codons are ", mutated_codons," and gene sequence is ",gene_seq))
        stop()
      }
      nb_sm_consensus <- 0 #initialization 
      nb_nsm_consensus <- 0 #initialization 
      #count the number of synonymous and non-synonymous mutations based on the genetic code
      if (original_codon == "TTT") {
        
        
        nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTC"))
        nb_nsm_consensus <- 1 - nb_sm_consensus
        
      } else if (original_codon == "TTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTG","CTT","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TTG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","CTT","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCC","TCA","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCA","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCC","TCG","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TCT","TCA","TCC","AGT","AGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TAA"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TAG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TGT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "TGA"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- 0
        
      } else if (original_codon == "TGG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- as.integer(consensus_codon_seq!=original_codon)
        
      } else if (original_codon == "CTT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTC","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTA","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTC","CTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CTG"){
        nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("TTA","TTG","CTT","CTC","CTA"))
        nb_nsm_consensus <- 1 - nb_sm_consensus
        
      } else if (original_codon == "CCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCC","CCA","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCA","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCC","CCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CCT","CCC","CCA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGC","CGA","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGA","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGC","CGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "CGG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("CGT","CGA","CGC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATC","ATA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATT","ATA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ATC","ATT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ATG"){
        nb_sm_consensus <- 0
        nb_nsm_consensus <- as.integer(consensus_codon_seq!=original_codon)
        
      } else if (original_codon == "ACT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACC","ACA","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACA","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACC","ACG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "ACG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("ACT","ACC","ACA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGC","TCT","TCC","TCA","TCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGT","TCT","TCC","TCA","TCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "AGG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("AGA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTA","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTT","GTA","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTT","GTG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GTG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GTC","GTA","GTT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCA","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCT","GCA","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCT","GCG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GCG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GCC","GCA","GCT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAC"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
      } else if (original_codon == "GAG"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GAA"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGT"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGA","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGC"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGT","GGA","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGA"){
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGT","GGG"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      } else if (original_codon == "GGG"){
        
        if(consensus_codon_seq!=original_codon){#if there is a mutation
          nb_sm_consensus <- as.integer(consensus_codon_seq %in% c("GGC","GGA","GGT"))
          nb_nsm_consensus <- 1 - nb_sm_consensus
        }else{
          nb_sm_consensus <- 0
          nb_nsm_consensus <- 0
        }
        
        
      }
      
      nb_pwdiff_codon <- calculate_nb_pwd(all_codons)
      Nb_pwdiff_gene <- Nb_pwdiff_gene + nb_pwdiff_codon
      nb_seg_site_current_codon <- get_nb_seg_sites_from_codon_seqs(all_codons)
      nb_segreg_sites <- nb_segreg_sites + nb_seg_site_current_codon
      nb_bp_diffs_current_codon <- get_nb_bp_diffs_consensus_codon_vs_orig(orig_codon = original_codon,consensus_seq = consensus_codon_seq)
      Nb_syn_mutations <- Nb_syn_mutations + (nb_sm_consensus*nb_bp_diffs_current_codon)
      Nb_nsyn_mutations <- Nb_nsyn_mutations + (nb_nsm_consensus*nb_bp_diffs_current_codon)
      
    }#end for loop codons
  }
  
  #Calculate Nb_nsyn_sites from Nb_syn_sites 
  Nb_nsyn_sites <- (gene_length - Nb_syn_sites)
  #get a1 and the number of copy of a certain gene during a certain generation
  nb_copy_gene_in_current_gen=nrow(df_the_gene_in_cells_at_generation_of_interest)
  nb_uniq_species_of_gene_in_current_gen=len(uniq(df_the_gene_in_cells_at_generation_of_interest$cell_ID))
  if (nb_copy_gene_in_current_gen == 0){
    a1_g <- NA
  }else{
    a1_g <- sum((1:(nb_copy_gene_in_current_gen-1))^-1) 
  }
  #parameters to calculate expected sqrt_variance
  a1_g <- sum((1:(nb_copy_gene_in_current_gen-1))^-1)
  a2_current_gene_during_g <- sum((1:(nb_copy_gene_in_current_gen-1))^-2)
  b1_current_gene_during_g <- (nb_copy_gene_in_current_gen+1)/(3*(nb_copy_gene_in_current_gen-1))
  b2_current_gene_during_g <- (2*((nb_copy_gene_in_current_gen^2)+nb_copy_gene_in_current_gen+3))/((9*nb_copy_gene_in_current_gen)*(nb_copy_gene_in_current_gen-1))
  c1_current_gene_during_g <- b1_current_gene_during_g - (1/a1_g)
  c2_current_gene_during_g <- b2_current_gene_during_g - ((nb_copy_gene_in_current_gen+2)/(a1_g*nb_copy_gene_in_current_gen)) + (a2_current_gene_during_g/(a1_g^2))
  e1_current_gene_during_g <- c1_current_gene_during_g/a1_g
  e2_current_gene_during_g <- c2_current_gene_during_g/((a1_g^2)+a2_current_gene_during_g)
  
  #Find S , find a1_g, calculate expected sqrt_variance with formula form Tajima (1989) and calculate Tajima's D and Ne_S_Taj according to it
  sqrt_expected_variance_current_gene_during_g <- sqrt((e1_current_gene_during_g*nb_segreg_sites)+((e2_current_gene_during_g*nb_segreg_sites)*(nb_segreg_sites-1)))
  
  #Ne_k_hat, Ne_S_Taj and Tajima's D
  Ne_k_hat_g_during_g <- Nb_pwdiff_gene/(2*mu*choose(k = 2,n = nb_copy_gene_in_current_gen))
  Ne_S_Taj_g_during_g <- nb_segreg_sites/(2*mu*a1_g)
  D_Taj_g_during_g <- ((Ne_k_hat_g_during_g - Ne_S_Taj_g_during_g)*(2*mu))/sqrt_expected_variance_current_gene_during_g
  
  #return values of interest
  return(list(nb_nsm=Nb_nsyn_mutations,nb_sm=Nb_syn_mutations,nb_nss=Nb_nsyn_sites,nb_ss=Nb_syn_sites,nb_segregative_sites=nb_segreg_sites,nb_copy_g_in_G=nb_copy_gene_in_current_gen,nb_uniq_species_of_g_in_G=nb_uniq_species_of_gene_in_current_gen,a1=a1_g,Ne_k_hat_gene_in_current_g=Ne_k_hat_g_during_g,Ne_S_Taj_gene_in_current_g=Ne_S_Taj_g_during_g,D_Taj_gene_in_current_g=D_Taj_g_during_g))
}


save_plots_workspace_current_sim <- paste0(results_output_workspace,simulation_name) #####DON't FORGET that variable "save_plots_workspace_current_sim" does not have "/" at the end
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("rm -rf ",save_plots_workspace_current_sim,"/"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("mkdir ",save_plots_workspace_current_sim,"/"),intern = FALSE,wait = TRUE)

#Unzip PANGENOMES_EVOLUTION_LOG
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("gunzip -fk ",sim_output_workpace,simulation_name,"/PANGENOMES_EVOLUTION_LOG.txt.gz"),intern = FALSE,wait = TRUE)

df_pangenome_evol_sim <- read.csv(file = paste0(sim_output_workpace,simulation_name,"/PANGENOMES_EVOLUTION_LOG.txt"),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
df_pangenome_evol_sim$cell_ID <- as.character(df_pangenome_evol_sim$cell_ID)
df_pangenome_evol_sim <- subset(x = df_pangenome_evol_sim, subset= fitness!=0) #remove dead cells
simulation_time <- max(df_pangenome_evol_sim$Generation_ctr,na.rm=TRUE)
max_genome_size <- max(df_pangenome_evol_sim$x,na.rm=TRUE)
min_genome_size <- min(df_pangenome_evol_sim$x,na.rm=TRUE)
nb_species_simulated <- length(unique(x = na.omit(df_pangenome_evol_sim$cell_ID)))
ggplot(data=df_pangenome_evol_sim,aes(x = x,colour = cell_ID,fill=cell_ID)) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of species genome size for simulation ",simulation_name),subtitle = paste0(" (",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; xi = ",average_nb_genes_per_genomes_at_start,"; x0 = ",wanted_x_0_at_equilibrium,")")) + xlab("Genome size") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) + scale_x_continuous(limits =c(min_genome_size,max_genome_size),breaks = seq(from = min_genome_size,to = max_genome_size,by = 5)) 
ggsave(filename = "Genome_size_distribution.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#r(x) in function of x
ggplot(data = df_pangenome_evol_sim, mapping = aes(x = x,y = r_x,col=cell_ID)) + geom_line() + ggtitle(label = paste0("Gene loss rate/Gene gain rate ratio r(x) in function of genome size x for simulation ",simulation_name)) + xlab("Genome size x") + ylab("Gene loss rate/Gene gain rate ratio r(x)")+ geom_hline(yintercept = 1, color = "black",lty=2) + scale_x_continuous(limits =c(min_genome_size,max_genome_size),breaks = seq(from = min_genome_size,to = max_genome_size,by = 5)) 
ggsave(filename = "r_x_VS_x.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#Beta(x) and Alpha(x) in function of x
ggplot(data = df_pangenome_evol_sim) + geom_line(mapping = aes(x = x,y = Beta_x,col="Beta(x)")) + geom_line(mapping = aes(x = x,y = Alpha_x,col="Alpha(x)")) + ggtitle(label = paste0("Gene loss rate Beta(x) and gene gain rate Alpha(x) in function of genome size x for simulation ",simulation_name)) + xlab("Genome size x") + ylab("Rates")+ scale_x_continuous(limits =c(min_genome_size,max_genome_size),breaks = seq(from = min_genome_size,to = max_genome_size,by = 5)) 
ggsave(filename = "Beta_x_and_Alpha_x_VS_x.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#Fitness(x) in function of x
ggplot(data = df_pangenome_evol_sim, mapping = aes(x = x,y = fitness,col=cell_ID)) + geom_point() + ggtitle(label = paste0("Fitness in function of genome size x for simulation ",simulation_name)) + xlab("Genome size x") + ylab("Fitness") + scale_x_continuous(limits =c(min_genome_size,max_genome_size),breaks = seq(from = min_genome_size,to = max_genome_size,by = 5)) 
ggsave(filename = "Fitness_VS_x.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

if (as.integer(num_fitness_landscape) == 9){
  #Plot the expected curve of Alpha_x and Beta_x in function of x (i.e. number of genes) with Koonin's model equations and the parameters I used for the type1 and type3 simulations
  expected_intercept_current_sim <- curve_intersect(curve1 = function(x) sPrime*(x^lambdaPlus),curve2 = function(x) rPrime*(x^lambdaMinus) , empirical = FALSE, domain = c(min(c(round(wanted_x_0_at_equilibrium/1.25),round(average_nb_genes_per_genomes_at_start/1.25)),na.rm = TRUE),max(c(round(wanted_x_0_at_equilibrium*1.25),round(average_nb_genes_per_genomes_at_start*1.25)),na.rm = TRUE)))
  ggplot(data.frame(x=seq(from = min(c(round(wanted_x_0_at_equilibrium/1.25),round(average_nb_genes_per_genomes_at_start/1.25)),na.rm = TRUE), to = max(c(round(wanted_x_0_at_equilibrium*1.25),round(average_nb_genes_per_genomes_at_start*1.25)),na.rm = TRUE),by = 1)), aes(x)) + stat_function(fun=function(x) sPrime*(x^lambdaPlus),col="blue") + stat_function(fun=function(x) rPrime*(x^lambdaMinus),col="red") + scale_x_continuous(limits =c(min(c(round(wanted_x_0_at_equilibrium/1.25),round(average_nb_genes_per_genomes_at_start/1.25)),na.rm = TRUE),max(c(round(wanted_x_0_at_equilibrium*1.25),round(average_nb_genes_per_genomes_at_start*1.25)),na.rm = TRUE)),breaks = seq(from = min(c(round(wanted_x_0_at_equilibrium/1.25),round(average_nb_genes_per_genomes_at_start/1.25)),na.rm = TRUE),to = max(c(round(wanted_x_0_at_equilibrium*1.25),round(average_nb_genes_per_genomes_at_start*1.25)),na.rm = TRUE),by = 1)) + geom_vline(xintercept = expected_intercept_current_sim$x, color = "grey",lty=2) + ggtitle(paste0("Koonin model expected loss rate (red) and gain rate (blue) during simulations; Equilibrium at x = ",expected_intercept_current_sim$x))
  ggsave(filename = paste0("Expected_loss_and_gain_rate_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
  
}

#For Genome size (x), r_x, Beta_x, Alpha_x and fitness calculate POPULATIONS min, mean and max
df_pangenome_evol_sim$pop_x_min <- summaryBy(formula = x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) min(the_vector,na.rm=TRUE),fun.names = "min",order = FALSE)$x.min
df_pangenome_evol_sim$pop_x_mean <- summaryBy(formula = x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) mean(the_vector,na.rm=TRUE),fun.names = "mean",order = FALSE)$x.mean
df_pangenome_evol_sim$pop_x_max <- summaryBy(formula = x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) max(the_vector,na.rm=TRUE),fun.names = "max",order = FALSE)$x.max

df_pangenome_evol_sim$pop_r_x_min <- summaryBy(formula = r_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) min(the_vector,na.rm=TRUE),fun.names = "min",order = FALSE)$r_x.min
df_pangenome_evol_sim$pop_r_x_mean <- summaryBy(formula = r_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) mean(the_vector,na.rm=TRUE),fun.names = "mean",order = FALSE)$r_x.mean
df_pangenome_evol_sim$pop_r_x_max <- summaryBy(formula = r_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) max(the_vector,na.rm=TRUE),fun.names = "max",order = FALSE)$r_x.max

df_pangenome_evol_sim$pop_Beta_x_min <- summaryBy(formula = Beta_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) min(the_vector,na.rm=TRUE),fun.names = "min",order = FALSE)$Beta_x.min
df_pangenome_evol_sim$pop_Beta_x_mean <- summaryBy(formula = Beta_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) mean(the_vector,na.rm=TRUE),fun.names = "mean",order = FALSE)$Beta_x.mean
df_pangenome_evol_sim$pop_Beta_x_max <- summaryBy(formula = Beta_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) max(the_vector,na.rm=TRUE),fun.names = "max",order = FALSE)$Beta_x.max

df_pangenome_evol_sim$pop_Alpha_x_min <- summaryBy(formula = Alpha_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) min(the_vector,na.rm=TRUE),fun.names = "min",order = FALSE)$Alpha_x.min
df_pangenome_evol_sim$pop_Alpha_x_mean <- summaryBy(formula = Alpha_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) mean(the_vector,na.rm=TRUE),fun.names = "mean",order = FALSE)$Alpha_x.mean
df_pangenome_evol_sim$pop_Alpha_x_max <- summaryBy(formula = Alpha_x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) max(the_vector,na.rm=TRUE),fun.names = "max",order = FALSE)$Alpha_x.max

df_pangenome_evol_sim$pop_fitness_min <- summaryBy(formula = fitness~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) min(the_vector,na.rm=TRUE),fun.names = "min",order = FALSE)$fitness.min
df_pangenome_evol_sim$pop_fitness_mean <- summaryBy(formula = fitness~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) mean(the_vector,na.rm=TRUE),fun.names = "mean",order = FALSE)$fitness.mean
df_pangenome_evol_sim$pop_fitness_max <- summaryBy(formula = fitness~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) max(the_vector,na.rm=TRUE),fun.names = "max",order = FALSE)$fitness.max

#Calculate Cell_count
df_pangenome_evol_sim$pop_count <- summaryBy(formula = x~Generation_ctr+cell_ID,data = df_pangenome_evol_sim,full.dimension = TRUE,FUN = function(the_vector) length(the_vector),fun.names = "count",order = FALSE)$x.count

# #Fixed_Loss_rate(x) eq.4 Koonin(2016) and Fixed_gain_rate(x) eq.3 Koonin(2016) in function of x
# ggplot(data = df_pangenome_evol_sim) + geom_line(mapping = aes(x = x,y = Beta_x*((-mean((a_for_sx+(b_for_sx*x)),na.rm=TRUE)*exp(-pop_count*-mean((a_for_sx+(b_for_sx*x)),na.rm=TRUE)))/(1-exp(-pop_count*-mean((a_for_sx+(b_for_sx*x)),na.rm=TRUE)))),col="Fixed_Loss_rate(x)")) + geom_line(mapping = aes(x = x,y = Alpha_x*((mean((a_for_sx+(b_for_sx*x)),na.rm=TRUE))/(1-exp(-pop_count*mean((a_for_sx+(b_for_sx*x)),na.rm=TRUE)))),col="Fixed_Gain_rate(x)")) + ggtitle(label = paste0("Fixed Loss rate (eq.4 Koonin 2016) and Fixed gain rate (eq.3 Koonin2016) in function of genome size x for simulation ",simulation_name)) + xlab("Genome size x") + ylab("Rates")+ scale_x_continuous(limits =c(min_genome_size,max_genome_size),breaks = seq(from = min_genome_size,to = max_genome_size,by = 5)) 
# ggsave(filename = "Koonin_Fixed_Beta_x_and_Fixed_Alpha_x_VS_x.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
# 
#Calculate the increase / decrease ratio of genome size for the different values of population genome size mean recorded during the simulation
df_pangenome_evol_sim$int_pop_x_mean <- floor(x = df_pangenome_evol_sim$pop_x_mean)
df_species_int_x_mean <- unique(df_pangenome_evol_sim[,c("Generation_ctr","int_pop_x_mean")])
current_sim_int_x_mean <- sort(unique(x = df_species_int_x_mean$int_pop_x_mean),decreasing = FALSE)
v_increase_vs_decrease_ratio <- rep(NA,length(current_sim_int_x_mean))
v_nb_events <- rep(NA,length(current_sim_int_x_mean))
pos_in_v_ratio <- 1
for (a_certain_int_x_mean in current_sim_int_x_mean){
  v_pos_occ_current_int_x_mean <- which(df_species_int_x_mean$int_pop_x_mean==a_certain_int_x_mean)
  v_next_x_mean_from_current_int_x_mean <- df_species_int_x_mean$int_pop_x_mean[(v_pos_occ_current_int_x_mean[v_pos_occ_current_int_x_mean!=nrow(df_species_int_x_mean)])+1]
  nb_increase_from_current_int_x_mean <- length(which(v_next_x_mean_from_current_int_x_mean>a_certain_int_x_mean))
  nb_decrease_from_current_int_x_mean <- length(which(v_next_x_mean_from_current_int_x_mean<a_certain_int_x_mean))
  ratio_increase_to_decrease_int_x_mean <- nb_increase_from_current_int_x_mean/nb_decrease_from_current_int_x_mean
  if (!is.na(ratio_increase_to_decrease_int_x_mean)){
    v_increase_vs_decrease_ratio[pos_in_v_ratio] <- ratio_increase_to_decrease_int_x_mean
    v_nb_events[pos_in_v_ratio] <- nb_increase_from_current_int_x_mean + nb_decrease_from_current_int_x_mean
  }
  pos_in_v_ratio <- pos_in_v_ratio + 1
  
}

ggplot() + geom_line(mapping = aes(x=current_sim_int_x_mean,y=v_increase_vs_decrease_ratio)) + xlab("Integer part of species average genome size") + ylab(paste0("Increase to decrease ratio (timestep = ",dt,")"))  + geom_vline(xintercept = wanted_x_0_at_equilibrium, color = "black",lty=2)
ggsave(filename = "Increase_to_decrease_ratio_plot.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

ggplot() + geom_line(mapping = aes(x=current_sim_int_x_mean,y=v_nb_events)) + xlab("Integer part of species average genome size") + ylab(paste0("Number of events considered (timestep = ",dt,")")) 
ggsave(filename = "Nb_events_for_increase_to_decrease_ratio_plot.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
#Cell Count Time Series
ggplot(data = df_pangenome_evol_sim,mapping = aes(x = Generation_ctr,y = pop_count)) +
  geom_area(aes(col=cell_ID,fill=cell_ID),alpha=0.3,position = position_dodge(0.8)) +
  xlab("Generation") +
  ylab("Cell Count") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "Cell_count_time_series.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
#Cell Count (log10-SCALE) Time Series
ggplot(data = df_pangenome_evol_sim,mapping = aes(x = Generation_ctr,y = log10(pop_count+(1E-16)))) +
  geom_area(aes(col=cell_ID,fill=cell_ID),alpha=0.3,position = position_dodge(0.8)) +
  xlab("Generation") +
  ylab("log10(Cell Count)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "Cell_count_log_scale_time_series.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#define time series overlapping points DODGE
pd <- position_dodge(dt/2)

#TIME SERIES WITH RANGE BARS
#time series for x
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=pop_x_mean, colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=pop_x_min, ymax=pop_x_max), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("Genome Size (x)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
#time series for log10(x)
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=log10(pop_x_mean+(1E-16)), colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=log10(pop_x_min+(1E-16)), ymax=log10(pop_x_max+(1E-16))), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("log10(x)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_log10_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#time series for Beta_x
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=pop_Beta_x_mean, colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=pop_Beta_x_min, ymax=pop_Beta_x_max), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("Beta_x_rate (number of loss events per generation per cells)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_Beta_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
#time series for log10(Beta_x)
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=log10(pop_Beta_x_mean+(1E-16)), colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=log10(pop_Beta_x_min+(1E-16)), ymax=log10(pop_Beta_x_max+(1E-16))), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("log10(Beta_x_rate)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_log10_Beta_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#time series for Alpha_x
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=pop_Alpha_x_mean, colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=pop_Alpha_x_min, ymax=pop_Alpha_x_max), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("Alpha_x_rate (number of gain events per generation per cells)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_Alpha_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
#time series for log10(Alpha_x)
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=log10(pop_Alpha_x_mean+(1E-16)), colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=log10(pop_Alpha_x_min+(1E-16)), ymax=log10(pop_Alpha_x_max+(1E-16))), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("log10(Alpha_x_rate)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_log10_Alpha_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")


if (!all(is.nan(df_pangenome_evol_sim$pop_r_x_mean))){
  #r_x (Beta_x/Alpha_x)
  ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=pop_r_x_mean, colour=cell_ID, group=cell_ID)) + 
    geom_errorbar(aes(ymin=pop_r_x_min, ymax=pop_r_x_max), width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3)+
    xlab("Generation") +
    ylab("r(x)") +
    scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename = paste0("Time_series_r_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
  #log10(r_x)
  ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=log10(pop_r_x_mean+(1E-16)), colour=cell_ID, group=cell_ID)) + 
    geom_errorbar(aes(ymin=log10(pop_r_x_min+(1E-16)), ymax=log10(pop_r_x_max+(1E-16))), width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3)+
    xlab("Generation") +
    ylab("log10(r(x))") +
    scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename = paste0("Time_series_log10_r_x_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
  
}

#time series for fitness
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=pop_fitness_mean, colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=pop_fitness_min, ymax=pop_fitness_max), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("Cell fitness") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_fitness_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")

#time series for log10(fitness)
ggplot(subset(x = df_pangenome_evol_sim,subset = Generation_ctr%in%c(1,seq(dt,simulation_time,20))), aes(x=Generation_ctr, y=log10(pop_fitness_mean+(1E-16)), colour=cell_ID, group=cell_ID)) + 
  geom_errorbar(aes(ymin=log10(pop_fitness_min+(1E-16)), ymax=log10(pop_fitness_max+(1E-16))), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3)+
  xlab("Generation") +
  ylab("log10(Cell fitness)") +
  scale_x_continuous(breaks = seq(0,(simulation_time),round(0.1*simulation_time)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("Time_series_log10_fitness_sim_",simulation_name,".png"), path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")


if (as.integer(num_fitness_landscape) != 7){
  #df_mutations_log_current_sim <- read.csv(file = paste0(sim_output_workpace,simulation_name,"/MUTATION_LOG"),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  
  system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("gunzip -fk ",sim_output_workpace,simulation_name,"/CELL_GENE_CONTENT_LOG.txt.gz"),intern = FALSE,wait = TRUE)
  
  df_cell_gene_content_current_sim <- read.csv(file = paste0(sim_output_workpace,simulation_name,"/CELL_GENE_CONTENT_LOG.txt"),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  v_generations <- c(1,seq(from=dt,to = simulation_time,by=dt))
  v_unique_gene_IDs_current_sim <- sort(unique(df_cell_gene_content_current_sim$gene_ID))
  df_results_variants_analysis_current_sim <- data.frame(gene_ID=rep(v_unique_gene_IDs_current_sim,length(v_generations)),Generation_ctr=rep(x = v_generations,each=length(v_unique_gene_IDs_current_sim)),nb_nss=0,nb_ss=0,nb_nsm=0,nb_sm=0,dn=NA,ds=NA,dnds=NA,nb_seg_sites=0,nb_copy_gene=0,nb_uniq_species_of_gene=0,a1=NA,Ne_S_Taj=0,Ne_k_hat=0,D_Taj=NA,stringsAsFactors = FALSE)
  #Only calculate variant analysis results for few time points
  df_results_variants_analysis_current_sim <- subset(x = df_results_variants_analysis_current_sim,subset = Generation_ctr %in% c(1,seq(from=(dt),to = simulation_time,by=(dt))))
  #number of cpus for variant analysis. Let at least 1 cpus free
  nb_cores_variant_analysis <- ifelse(test=nb_cores_variant_analysis>=detectCores(),yes=detectCores()-1,no=nb_cores_variant_analysis)
  #Keep a number of rows in that is a multiple of nb_cores_variant_analysis 
  #while (nrow(df_results_variants_analysis_current_sim)%%nb_cores_variant_analysis!=0){
    #nb_cores_variant_analysis = nb_cores_variant_analysis - 1
  #}
  # Setup cluster
  cl <- makeCluster(nb_cores_variant_analysis,outfile=paste0("LOGFILE_VARIANT_ANALYSIS_sim_",simulation_name,".txt"))
  registerDoParallel(cl)
  #variants analysis main loop function
  main_loop_variants_analysis <- function(p_row_ind_begin,p_row_ind_end){
    system(paste0("touch ",save_plots_workspace_current_sim,"/LOGFILE_VARIANT_ANALYSIS_sim_",simulation_name,"_process",Sys.getpid(),".txt"))
    for (p_row_ind in p_row_ind_begin:p_row_ind_end){
      synonymity_infos <- get_synonymity_stats(the_df_cell_content_log = df_cell_gene_content_current_sim,the_gene_ID = df_results_variants_analysis_current_sim$gene_ID[p_row_ind],generation_of_interest = df_results_variants_analysis_current_sim$Generation_ctr[p_row_ind])
      nb_copy_and_a1_infos <- list(nb_copy_g_in_G=synonymity_infos$nb_copy_g_in_G,a1=synonymity_infos$a1)
      nb_uniq_species_of_current_g_in_current_g = synonymity_infos$nb_uniq_species_of_g_in_G
      #extract all the infos/variables from the two previous lists AND ASSIGN VARIABLES VALUES FOR THE ROW at p_row_ind  
      nb_nsm_current_g_during_g <- synonymity_infos$nb_nsm
      df_results_variants_analysis_current_sim$nb_nsm[p_row_ind] <- nb_nsm_current_g_during_g
      nb_sm_current_g_during_g <- synonymity_infos$nb_sm
      df_results_variants_analysis_current_sim$nb_sm[p_row_ind] <- nb_sm_current_g_during_g
      nb_nss_current_g_during_g <- synonymity_infos$nb_nss
      df_results_variants_analysis_current_sim$nb_nss[p_row_ind] <- nb_nss_current_g_during_g
      nb_ss_current_g_during_g <- synonymity_infos$nb_ss
      df_results_variants_analysis_current_sim$nb_ss[p_row_ind] <- nb_ss_current_g_during_g
      dn_current_g_during_g <- nb_nsm_current_g_during_g/nb_nss_current_g_during_g
      df_results_variants_analysis_current_sim$dn[p_row_ind] <- dn_current_g_during_g
      ds_current_g_during_g <- nb_sm_current_g_during_g/nb_ss_current_g_during_g
      df_results_variants_analysis_current_sim$ds[p_row_ind] <- ds_current_g_during_g
      
      if ((!is.na(nb_sm_current_g_during_g))&(nb_sm_current_g_during_g>0)){
        df_results_variants_analysis_current_sim$dnds[p_row_ind] <- df_results_variants_analysis_current_sim$dn[p_row_ind]/df_results_variants_analysis_current_sim$ds[p_row_ind]
      }
      
      nb_seg_sites_current_g_during_g <- synonymity_infos$nb_segregative_sites
      df_results_variants_analysis_current_sim$nb_seg_sites[p_row_ind] <- nb_seg_sites_current_g_during_g
      nb_copy_current_g_during_g <- nb_copy_and_a1_infos$nb_copy_g_in_G
      df_results_variants_analysis_current_sim$nb_copy_gene[p_row_ind] <- nb_copy_current_g_during_g
      a1_current_g_during_g <- nb_copy_and_a1_infos$a1
      df_results_variants_analysis_current_sim$a1[p_row_ind] <- a1_current_g_during_g
      df_results_variants_analysis_current_sim$nb_uniq_species_of_gene[p_row_ind] <- nb_uniq_species_of_current_g_in_current_g
      Ne_S_Taj_current_g_during_g <-synonymity_infos$Ne_S_Taj_gene_in_current_g
      df_results_variants_analysis_current_sim$Ne_S_Taj[p_row_ind] <- Ne_S_Taj_current_g_during_g
      df_results_variants_analysis_current_sim$Ne_k_hat[p_row_ind] <- synonymity_infos$Ne_k_hat_gene_in_current_g
      df_results_variants_analysis_current_sim$D_Taj[p_row_ind] <- synonymity_infos$D_Taj_gene_in_current_g
      msg_surrent_state<- (paste0(floor(x = ((p_row_ind-p_row_ind_begin)*100/(p_row_ind_end-p_row_ind_begin+1))),"% of the variants analysis has been done for simulation_",simulation_name," for process ",Sys.getpid()," !"))
      system(paste0("echo ",msg_surrent_state," >> ",save_plots_workspace_current_sim,"/LOGFILE_VARIANT_ANALYSIS_sim_",simulation_name,"_process",Sys.getpid(),".txt"))
      
    }
    return(df_results_variants_analysis_current_sim[p_row_ind_begin:p_row_ind_end,])
  }
  
  v_row_ind_begin <- seq(from=1,to=nrow(df_results_variants_analysis_current_sim),by=nrow(df_results_variants_analysis_current_sim)/nb_cores_variant_analysis)
  v_row_ind_end <- seq(from=nrow(df_results_variants_analysis_current_sim)/nb_cores_variant_analysis,to=nrow(df_results_variants_analysis_current_sim),by=nrow(df_results_variants_analysis_current_sim)/nb_cores_variant_analysis)
  df_results_variants_analysis_current_sim <- foreach(the_p_row_ind_begin = v_row_ind_begin,the_p_row_ind_end = v_row_ind_end, .combine = rbind,.export=c("df_cell_gene_content_current_sim","df_results_variants_analysis_current_sim"),.packages=c("ggplot2","nlme","session","doBy","reconPlots","gridExtra"))  %dopar%  main_loop_variants_analysis(p_row_ind_begin = the_p_row_ind_begin,p_row_ind_end = the_p_row_ind_end)
  
  stopCluster(cl)
  df_results_variants_analysis_current_sim
  #Plots :
  
  #Get gene variants information for the maximal generation at which we have gene sequence variants
  df_results_sim_at_max_generation_during_which_variants_exist <- subset(x = df_results_variants_analysis_current_sim,subset = Generation_ctr==max(subset(x=df_results_variants_analysis_current_sim,subset=((!is.na(nb_sm)) & (nb_sm!=0)))$Generation_ctr,na.rm=TRUE))
  
  nb_genes_variant <- nrow(subset(x = df_results_sim_at_max_generation_during_which_variants_exist,subset=((!is.na(nb_sm)) & (nb_sm!=0))))
  nb_genes_variant
  if (nb_genes_variant > 1){
    #dn/ds distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = dnds,colour = "red",fill="red")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes dN/dS at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("dN/dS") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "dnds_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #log10(dn/ds) distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = log10(dnds+(1E-16)),colour = "red",fill="red")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes log10(dN/dS) at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("log10(dN/dS)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "log10_dnds_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #Ne_S_Taj distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = Ne_S_Taj,colour = "blue",fill="blue")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes Ne_S_Taj at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("Ne_S_Taj") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "Ne_S_Taj_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #log10(Ne_S_Taj) distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = log10(Ne_S_Taj+(1E-16)),colour = "blue",fill="blue")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes log10(Ne_S_Taj) at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("log10(Ne_S_Taj)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "log10_Ne_S_Taj_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #Ne_k_hat distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = Ne_k_hat,colour = "blue",fill="blue")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes Ne_k_hat at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("Ne_k_hat") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "Ne_k_hat_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #log10(Ne_k_hat) distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = log10(Ne_k_hat+(1E-16)),colour = "blue",fill="blue")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes log10(Ne_k_hat) at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("log10(Ne_k_hat)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "log10_Ne_k_hat_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
    #D_Taj distribution AT GENERATION == maximum of Generation during which gene variants exist
    ggplot(data=df_results_sim_at_max_generation_during_which_variants_exist,aes(x = D_Taj,colour = "blue",fill="blue")) + geom_density(alpha=0.1) + ggtitle(paste0("Distribution of genes D_Taj at the end of simulation_",simulation_name),subtitle = paste0(" (Simulation duration = ",simulation_time," generations; ",nb_species_simulated," species; N = ",nb_cells_simulated," cells; ",average_nb_genes_per_genomes_at_start," genes per cell in average at the start of the simulation)")) + xlab("D_Taj") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=12)) #+ scale_x_continuous(limits =c(,),breaks = seq(from = ,to = ,by = )) 
    ggsave(filename = "D_Taj_distribution_AT_THE_END_OF_SIMULATION.png", path=save_plots_workspace_current_sim, width = 15, height = 20, units = "cm")
    
  }
  ##########
  #initialization
  Slope_lm_Ne_k_hat_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  Slope_lm_Ne_S_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  Slope_lm_D_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  pvalue_lm_Ne_k_hat_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  pvalue_lm_Ne_S_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  pvalue_lm_D_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  adj_r_square_lm_Ne_k_hat_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  adj_r_square_lm_Ne_S_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  adj_r_square_lm_D_Taj_vs_Gene_Mobility <- rep(NA,length(v_unique_gene_IDs_current_sim))
  p_val_threshold <- 0.05/length(v_unique_gene_IDs_current_sim)
  #####DON't FORGET that variable "save_plots_workspace_current_sim" does not have "/" at the end
  for (v_g_ind in 1:length(v_unique_gene_IDs_current_sim)){
    list_results_current_lm_Ne_k_hat_vs_Gene_Mobility <- list(adj_r_sq_current_lm = NA,slope_current_lm = NA,p_val_current_lm=NA)
    list_results_current_lm_Ne_S_Taj_vs_Gene_Mobility <- list(adj_r_sq_current_lm = NA,slope_current_lm = NA,p_val_current_lm=NA)
    list_results_current_lm_D_Taj_vs_Gene_Mobility <- list(adj_r_sq_current_lm = NA,slope_current_lm = NA,p_val_current_lm=NA)
    
    save_plots_workspace_current_gene <- paste0(save_plots_workspace_current_sim,"/gene_",v_unique_gene_IDs_current_sim[v_g_ind])
    system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("rm -rf ",save_plots_workspace_current_gene,"/"),intern = FALSE,wait = TRUE)
    system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("mkdir ",save_plots_workspace_current_gene,"/"),intern = FALSE,wait = TRUE)
    
    df_current_gene_in_time <- subset(x = df_results_variants_analysis_current_sim,subset = gene_ID==v_unique_gene_IDs_current_sim[v_g_ind])
    
    if (length(na.omit(df_current_gene_in_time$dnds))>=1){
      #dN/dS in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_dnds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(df_current_gene_in_time$dnds~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene dN/dS",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$dnds,na.rm = TRUE),max(df_current_gene_in_time$dnds,na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
      #log10(dN/dS) in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_log10_dnds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(log10(df_current_gene_in_time$dnds+(1E-16))~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene log10(dN/dS)",col="red",xlim=c(1,simulation_time+1),ylim=c(min(log10(df_current_gene_in_time$dnds+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$dnds+(1E-16)),na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
    }
    
    if (length(na.omit(df_current_gene_in_time$dn))>=1){
      #dn in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_dn_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(df_current_gene_in_time$dn~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene dn",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$dn,na.rm = TRUE),max(df_current_gene_in_time$dn,na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
      #log10(dn) in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_log10_dn_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(log10(df_current_gene_in_time$dn+(1E-16))~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene log10(dn)",col="red",xlim=c(1,simulation_time+1),ylim=c(min(log10(df_current_gene_in_time$dn+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$dn+(1E-16)),na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
    }
    
    if (length(na.omit(df_current_gene_in_time$ds))>=1){
      #ds in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_ds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(df_current_gene_in_time$ds~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene ds",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$ds,na.rm = TRUE),max(df_current_gene_in_time$ds,na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
      #log10(ds) in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_log10_ds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(log10(df_current_gene_in_time$ds+(1E-16))~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene log10(ds)",col="red",xlim=c(1,simulation_time+1),ylim=c(min(log10(df_current_gene_in_time$ds+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$ds+(1E-16)),na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
      if (length(intersect(which(is.finite(df_current_gene_in_time$dn)),which(is.finite(df_current_gene_in_time$ds))))>=2){
        #dn vs ds 
        v_pos_for_current_reg <-intersect(which(is.finite(df_current_gene_in_time$dn)),which(is.finite(df_current_gene_in_time$ds)))
        tryCatch(expr = ggplotRegression(fit=lmp(dn ~ ds, data = df_current_gene_in_time[v_pos_for_current_reg,],na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("dn_vs_ds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "ds",ylabl = "dn"),error=function(e) print(paste0("Regression skipped because of following error : ",e)))
      
      }
    }

    if (length(na.omit(df_current_gene_in_time$Ne_S_Taj))>=1){
   
        #Ne_S_Taj in function of Generation plot
        png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_Ne_S_Taj_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
        plot(df_current_gene_in_time$Ne_S_Taj~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene Ne_S_Taj",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$Ne_S_Taj,na.rm = TRUE),max(df_current_gene_in_time$Ne_S_Taj,na.rm=TRUE)))
        axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
        dev.off()
        #log10(Ne_S_Taj) in function of Generation plot
        png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_log10_Ne_S_Taj_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
        plot(replace(x = log10(df_current_gene_in_time$Ne_S_Taj+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_S_Taj+(1E-16))==-16),values = NA)~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene log10(Ne_S_Taj)",col="red",xlim=c(1,simulation_time+1),ylim=c(min(replace(x = log10(df_current_gene_in_time$Ne_S_Taj+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_S_Taj+(1E-16))==-16),values = NA),na.rm = TRUE),max(replace(x = log10(df_current_gene_in_time$Ne_S_Taj+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_S_Taj+(1E-16))==-16),values = NA),na.rm=TRUE)))
        axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
        dev.off()
      
    }
    
    #Ne_S_Taj vs Gene_Mobility (STANDARDIZED)
    current_y <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$Ne_S_Taj)
    current_x <-scale(center=TRUE,scale=TRUE,df_current_gene_in_time$nb_uniq_species_of_gene)
    v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))

    if (length(intersect(which(is.finite(df_current_gene_in_time$Ne_S_Taj)),which(is.finite(df_current_gene_in_time$nb_uniq_species_of_gene))))>=2){
      tryCatch(expr ={list_results_current_lm_Ne_S_Taj_vs_Gene_Mobility=ggplotRegression(fit=lmp(Ne_S_Taj ~ nb_uniq_species_of_gene, data = as.data.frame(scale(center=TRUE,scale=TRUE,df_current_gene_in_time)[v_pos_for_current_reg,]),na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("standardized_Ne_S_Taj_vs_Gene_Mobility_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "standardized Gene_Mobility",ylabl = "standardized Ne_S_Taj")},error=function(e) print(paste0("Regression skipped because of following error : ",e)))
    }
    if (length(na.omit(df_current_gene_in_time$Ne_k_hat))>=1){
      
        #Ne_k_hat in function of Generation plot
        png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_Ne_k_hat_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
        plot(df_current_gene_in_time$Ne_k_hat~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene Ne_k_hat",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$Ne_k_hat,na.rm = TRUE),max(df_current_gene_in_time$Ne_k_hat,na.rm=TRUE)))
        axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
        dev.off()
        #log10(Ne_k_hat) in function of Generation plot
        png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_log10_Ne_k_hat_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
        plot(replace(x = log10(df_current_gene_in_time$Ne_k_hat+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_k_hat+(1E-16))==-16),values = NA)~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene log10(Ne_k_hat)",col="red",xlim=c(1,simulation_time+1),ylim=c(min(replace(x = log10(df_current_gene_in_time$Ne_k_hat+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_k_hat+(1E-16))==-16),values = NA),na.rm = TRUE),max(replace(x = log10(df_current_gene_in_time$Ne_k_hat+(1E-16)),list = which(log10(df_current_gene_in_time$Ne_k_hat+(1E-16))==-16),values = NA),na.rm=TRUE)))
        axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
        dev.off()
      
      
    }
      
      #Ne_k_hat vs Gene_Mobility (STANDARDIZED)
      current_y <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$Ne_k_hat)
      current_x <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$nb_uniq_species_of_gene)
      v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))

    if (length(intersect(which(is.finite(df_current_gene_in_time$Ne_k_hat)),which(is.finite(df_current_gene_in_time$nb_uniq_species_of_gene))))>=2){
      tryCatch(expr ={list_results_current_lm_Ne_k_hat_vs_Gene_Mobility=ggplotRegression(fit=lmp(Ne_k_hat ~ nb_uniq_species_of_gene, data = as.data.frame(scale(center=TRUE,scale=TRUE,df_current_gene_in_time)[v_pos_for_current_reg,]),na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("standardized_Ne_k_hat_vs_Gene_Mobility_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "standardized Gene_Mobility",ylabl = "standardized Ne_k_hat")},error=function(e) print(paste0("Regression skipped because of following error : ",e)))
    }
    
    
    if (length(na.omit(df_current_gene_in_time$D_Taj))>=1){
      #D_Taj in function of Generation plot
      png(filename = paste0(save_plots_workspace_current_gene,"/Time_series_D_Taj_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
      plot(df_current_gene_in_time$D_Taj~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene D_Taj",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$D_Taj,na.rm = TRUE),max(df_current_gene_in_time$D_Taj,na.rm=TRUE)))
      axis(side=1, at=seq(0, simulation_time+dt, by=(2*dt)),las=2)
      dev.off()
      if (length(na.omit(df_current_gene_in_time$dnds))>=1){
        #Time series of dn/ds, Tajima's D and Gene_Mobility (Number of unique species containing the gene)
        png(filename = paste0(save_plots_workspace_current_gene,"/dnds_Tajima_D_and_Gene_Mobility_combined_time_series_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
        par(mfrow=c(3,1))
        plot(df_current_gene_in_time$dnds~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene dN/dS",col="blue",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$dnds,na.rm = TRUE),max(df_current_gene_in_time$dnds,na.rm=TRUE)))
        axis(side=1, c(1,seq(from=(dt*25),to = simulation_time,by=(dt*25))),las=2)
        plot(df_current_gene_in_time$D_Taj~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene Tajima's D",col="red",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$D_Taj,na.rm = TRUE),max(df_current_gene_in_time$D_Taj,na.rm=TRUE)))
        axis(side=1, at=c(1,seq(from=(dt*25),to = simulation_time,by=(dt*25))),las=2)
        plot(df_current_gene_in_time$nb_uniq_species_of_gene~df_current_gene_in_time$Generation_ctr, type="p", pch=19,xaxt="n",xlab="Generation",ylab="Gene Mobility",col="green",xlim=c(1,simulation_time+1),ylim=c(min(df_current_gene_in_time$nb_uniq_species_of_gene,na.rm = TRUE),max(df_current_gene_in_time$nb_uniq_species_of_gene,na.rm=TRUE)))
        axis(side=1, at=c(1,seq(from=(dt*25),to = simulation_time,by=(dt*25))),las=2)
        dev.off()
        
        par(mfrow=c(1,1))
        if (length(intersect(which(is.finite(df_current_gene_in_time$D_Taj)),which(is.finite(df_current_gene_in_time$nb_uniq_species_of_gene))))>=2){
          ##D_Taj in function of Gene_Mobility (STANDARDIZED) plot 
          # png(filename = paste0(save_plots_workspace_current_gene,"/Tajima_D_VS_log10_Gene_Mobility_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
          # plot(df_current_gene_in_time$D_Taj~log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16)), type="p", pch=19,xlab="log10(Gene_Mobility)",ylab="Tajima_D",col="red",xlim=c(min(log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16)),na.rm=TRUE)),ylim=c(min(df_current_gene_in_time$D_Taj,na.rm = TRUE),max(df_current_gene_in_time$D_Taj,na.rm=TRUE)))
          # dev.off()
          current_y <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$D_Taj)
          current_x <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$nb_uniq_species_of_gene)
          v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))
          tryCatch(expr ={list_results_current_lm_D_Taj_vs_Gene_Mobility=ggplotRegression(fit=lmp(D_Taj~nb_uniq_species_of_gene , data = as.data.frame(scale(center=TRUE,scale=TRUE,df_current_gene_in_time)[v_pos_for_current_reg,]),na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("standardized_Tajima_D_VS_Gene_Mobility_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "standardized Gene_Mobility",ylabl = "standardized Tajima_D")},error=function(e) print(paste0("Regression skipped because of following error : ",e)))
        }
      
        if (length(intersect(which(is.finite(df_current_gene_in_time$D_Taj)),which(is.finite(df_current_gene_in_time$dnds))))>=2){
          ##D_Taj in function of dN/dS plot (STANDARDIZED)
          # png(filename = paste0(save_plots_workspace_current_gene,"/Tajima_D_VS_log10_dnds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
          # plot(df_current_gene_in_time$D_Taj~log10(df_current_gene_in_time$dnds+(1E-16)), type="p", pch=19,xlab="log10(dN/dS)",ylab="Tajima_D",col="red",xlim=c(min(log10(df_current_gene_in_time$dnds+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$dnds+(1E-16)),na.rm=TRUE)),ylim=c(min(df_current_gene_in_time$D_Taj,na.rm = TRUE),max(df_current_gene_in_time$D_Taj,na.rm=TRUE)))
          # dev.off()
          df_current_gene_in_time <- subset(x = df_results_variants_analysis_current_sim,subset = (gene_ID==v_unique_gene_IDs_current_sim[v_g_ind])&(!is.na(dnds)))
          current_y <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$D_Taj)
          current_x <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$dnds)
          v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))
          tryCatch(expr = ggplotRegression(fit=lmp(D_Taj~dnds, data = as.data.frame(scale(center=TRUE,scale=TRUE,df_current_gene_in_time)[v_pos_for_current_reg,]),na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("standardized_Tajima_D_VS_dnds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "standardized dN/dS",ylabl = "standardized Tajima_D"),error=function(e) print(paste0("Regression skipped because of following error : ",e)))
        }
        if (length(intersect(which(is.finite(df_current_gene_in_time$nb_uniq_species_of_gene)),which(is.finite(df_current_gene_in_time$dnds))))>=2){
          ##Gene_Mobility in function of dN/dS plot (STANDARDIZED)
          # png(filename = paste0(save_plots_workspace_current_gene,"/Gene_Mobility_VS_dnds_Original_values_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
          # plot(df_current_gene_in_time$nb_uniq_species_of_gene~df_current_gene_in_time$dnds, type="p", pch=19,xlab="dN/dS",ylab="Gene_Mobility",col="red",ylim=c(min(df_current_gene_in_time$nb_uniq_species_of_gene,na.rm = TRUE),max(df_current_gene_in_time$nb_uniq_species_of_gene,na.rm=TRUE)))
          # dev.off()
          current_y <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$nb_uniq_species_of_gene)
          current_x <- scale(center=TRUE,scale=TRUE,df_current_gene_in_time$dnds)
          v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))
          tryCatch(expr =ggplotRegression(fit=lmp(nb_uniq_species_of_gene~dnds, data = as.data.frame(scale(center=TRUE,scale=TRUE,df_current_gene_in_time)[v_pos_for_current_reg,]),na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("standardized_Gene_Mobility_VS_dnds_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "dN/dS",ylabl = "Gene_Mobility"),error=function(e) print(paste0("Regression skipped because of following error : ",e)))
        }
        #if (length(intersect(which(is.finite(df_current_gene_in_time$nb_uniq_species_of_gene)),which(is.finite(df_current_gene_in_time$dnds))))>=2){
          ##log10(Gene_Mobility) in function of log10(dN/dS) plot
          # png(filename = paste0(save_plots_workspace_current_gene,"/Gene_Mobility_VS_dnds_log10_scale_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),width = 1000, height = 700, units = "px")
          # plot(log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16))~log10(df_current_gene_in_time$dnds+(1E-16)), type="p", pch=19,xlab="log10(dN/dS)",ylab="log10(Gene_Mobility)",col="red",ylim=c(min(log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16)),na.rm = TRUE),max(log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16)),na.rm=TRUE)))
          # dev.off()
          # current_y <- log10(df_current_gene_in_time$nb_uniq_species_of_gene+(1E-16))
          # current_x <- log10(df_current_gene_in_time$dnds+(1E-16))
          # v_pos_for_current_reg <-intersect(which(is.finite(current_y)),which(is.finite(current_x)))
          # tryCatch(expr = ggplotRegression(fit=lmp(nb_uniq_species_of_gene ~ dnds, data = log10(df_current_gene_in_time+(1E-16))[v_pos_for_current_reg,],na.action=na.omit),ggsave_path=save_plots_workspace_current_gene,the_filename=paste0("Gene_Mobility_VS_dnds_in_log10_scale_sim_",simulation_name,"_gene_",v_unique_gene_IDs_current_sim[v_g_ind],".png"),xlabl = "log10(dN/dS)",ylabl = "log10(Gene_Mobility)"),error=function(e) print(paste0("Regression skipped because of following error : ",e)))
        #}
      }
    }
    #Linear regressions explained by Gene Mobility
    adj_r_square_lm_Ne_k_hat_vs_Gene_Mobility[v_g_ind] <-list_results_current_lm_Ne_k_hat_vs_Gene_Mobility$adj_r_sq_current_lm
    adj_r_square_lm_Ne_S_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_Ne_S_Taj_vs_Gene_Mobility$adj_r_sq_current_lm
    adj_r_square_lm_D_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_D_Taj_vs_Gene_Mobility$adj_r_sq_current_lm
    Slope_lm_Ne_k_hat_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_Ne_k_hat_vs_Gene_Mobility$slope_current_lm
    Slope_lm_Ne_S_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_Ne_S_Taj_vs_Gene_Mobility$slope_current_lm
    Slope_lm_D_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_D_Taj_vs_Gene_Mobility$slope_current_lm
    pvalue_lm_Ne_k_hat_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_Ne_k_hat_vs_Gene_Mobility$p_val_current_lm
    pvalue_lm_Ne_S_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_Ne_S_Taj_vs_Gene_Mobility$p_val_current_lm
    pvalue_lm_D_Taj_vs_Gene_Mobility[v_g_ind] <- list_results_current_lm_D_Taj_vs_Gene_Mobility$p_val_current_lm
  }
  
  #Summary Figures of significant Gene_Mobility effects
  ggplot() + geom_density(mapping = aes(x = Slope_lm_Ne_k_hat_vs_Gene_Mobility[pvalue_lm_Ne_k_hat_vs_Gene_Mobility<=(p_val_threshold)],colour="Ne_k_hat")) + geom_density(mapping = aes(x = Slope_lm_D_Taj_vs_Gene_Mobility[pvalue_lm_D_Taj_vs_Gene_Mobility<=(p_val_threshold)],colour="D_Taj")) + geom_density(mapping = aes(x = Slope_lm_Ne_S_Taj_vs_Gene_Mobility[pvalue_lm_Ne_S_Taj_vs_Gene_Mobility<=(p_val_threshold)],colour="Ne_S_Taj")) + ggtitle(label = "Distribution of Slopes for regression models in which Gene Mobility is the explanatory variable") + xlab("Slope") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=9),plot.subtitle = element_text(hjust=0,size=8)) + geom_vline(xintercept = 0,mapping = aes(x = Slope_lm_Ne_k_hat_vs_Gene_Mobility[pvalue_lm_Ne_k_hat_vs_Gene_Mobility<=(p_val_threshold)])) + scale_color_manual(name="X ~ Gene Mobility",values = c(Ne_k_hat="red",Ne_S_Taj="blue",D_Taj="black")) + xlim(-3,3) + ylim(0,10)
  ggsave(filename = "Gene_Mobility_significant_slopes.png",path = save_plots_workspace_current_sim,width = 18, height = 10, units = "cm")
  
  ggplot() + geom_boxplot(aes(y=c(adj_r_square_lm_Ne_k_hat_vs_Gene_Mobility,adj_r_square_lm_Ne_S_Taj_vs_Gene_Mobility,adj_r_square_lm_D_Taj_vs_Gene_Mobility),x=rep(x = c("Ne_k_hat","Ne_S_Taj","D_Taj"),each=length(adj_r_square_lm_Ne_k_hat_vs_Gene_Mobility)))) + xlab("Response variable") + ylab("Adjusted R_squared")
  ggsave(filename = "Gene_Mobility_significant_adj_R_squared.png",path = save_plots_workspace_current_sim,width = 18, height = 10, units = "cm")  
  
  #save the result matrix of the CURRENT SIMULATION
  mtx_results_current_sim <- as.matrix(df_results_variants_analysis_current_sim)
  write.table(mtx_results_current_sim, file = paste0(save_plots_workspace_current_sim,"/mtx_variants_analysis_results_simulation_",simulation_name,".csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
  
}
#Interesting genes to look as variant analysis results
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("rm -f ",save_plots_workspace_current_sim,"/intersting_genes_variant_analysis.txt"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("touch ",save_plots_workspace_current_sim,"/intersting_genes_variant_analysis.txt"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("echo \'",paste(v_unique_gene_IDs_current_sim[which(pvalue_lm_Ne_S_Taj_vs_Gene_Mobility<=(p_val_threshold))],collapse=';'),"\' >> ",save_plots_workspace_current_sim,"/intersting_genes_variant_analysis.txt"),intern = FALSE,wait = TRUE)


#Only keep .gz files
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("rm -f ",sim_output_workpace,simulation_name,"/*.txt"),intern = FALSE,wait = TRUE)

########
#Save RSession
save.session(paste0(results_output_workspace,"simulation_",simulation_name,"_RSession.Rda"))


# #Restore Session
# library("session")
# results_output_workspace <- "/home/arnaudng/SodaPop_big_sims_with_HGT_and_loss/Sim_Results/"
# restore.session(paste0(results_output_workspace,"simulation_",simulation_name,"_RSession.Rda"))
