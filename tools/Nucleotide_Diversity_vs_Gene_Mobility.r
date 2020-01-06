############################################################## Script Header ##################################################
#@Author = Arnaud NG
#This script produces plots to assess the effect of the "HGT space" on the correlation between nucleotide diversity vs gene mobility

#LOAD PACKAGES
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
library("gplots")

#INITIALIZATION
sodapop_workspace <- as.character(commandArgs(TRUE)[1]) #ABSOLUTE Path of the folder containing SodaPop files (should be the parent of directory of SodaPop "out/" or "src/")
sodapop_workspace <- ifelse(test=substr(x = sodapop_workspace,start = nchar(sodapop_workspace),stop = nchar(sodapop_workspace))=="/",yes=sodapop_workspace,no=paste0(sodapop_workspace,"/"))
sim_output_workpace <- paste0(sodapop_workspace,"out/") #ABSOLUTE path of the folder containing Sodapop simulations raw output repertory 
results_output_workspace <- paste0(sodapop_workspace,"Sim_Results/") #ABSOLUTE path of the folder containing the Post-simulation genes variants analysis results. Will be created if it does not exist yet.
hgt_rates <- c("1mu","10mu","100mu","1000mu")
lambda_exp_distr_s_hgt <- c("neutral","1E3","1E6","1E9")
mtx_slopes_Ne_k_hat_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_slopes_Ne_k_hat_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_slopes_Ne_k_hat_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
mtx_slopes_Ne_S_Taj_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_slopes_Ne_S_Taj_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_slopes_Ne_S_Taj_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
mtx_rsq_Ne_k_hat_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_rsq_Ne_k_hat_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_rsq_Ne_k_hat_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
mtx_rsq_Ne_S_Taj_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_rsq_Ne_S_Taj_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_rsq_Ne_S_Taj_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
mtx_pval_Ne_k_hat_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_pval_Ne_k_hat_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_pval_Ne_k_hat_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
mtx_pval_Ne_S_Taj_vs_Gene_Mobility <- matrix(data = NA,nrow = 4,ncol = 4)
rownames(mtx_pval_Ne_S_Taj_vs_Gene_Mobility) <- hgt_rates
colnames(mtx_pval_Ne_S_Taj_vs_Gene_Mobility) <- lambda_exp_distr_s_hgt
lst_current_sim_Ne_k_hat_vs_Gene_mobility <- list(adj_r_sq_current_lm=NA,slope_current_lm=NA,p_val_current_lm=NA)
lst_current_sim_Ne_S_Taj_vs_Gene_mobility <- list(adj_r_sq_current_lm=NA,slope_current_lm=NA,p_val_current_lm=NA)

for (i in 1:length(hgt_rates)){
  for (j in 1:length(lambda_exp_distr_s_hgt)){
    simulation_name <- paste0("simTenSpeciesHGT_",hgt_rates[i],"_",lambda_exp_distr_s_hgt[j])
    df_current_sim <- read.csv(file = paste0("D:/Bureau/SodaPop_test_new/Sim_Results/mtx_variants_analysis_results_simulation_",simulation_name,".csv"),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    df_subset_near_end_sim_current_reg <- subset(df_current_sim,Generation_ctr%in%(sort(unique(Generation_ctr),decreasing = TRUE)[sample(1:3,1)]))
    df_subset_near_end_sim_current_reg <- subset(df_subset_near_end_sim_current_reg,nb_uniq_species_of_gene!=0)
    df_subset_near_end_sim_current_reg$Ne_k_hat <- log10(df_subset_near_end_sim_current_reg$Ne_k_hat+(1e-16)) 
    df_subset_near_end_sim_current_reg$nb_uniq_species_of_gene <- log10(df_subset_near_end_sim_current_reg$nb_uniq_species_of_gene+(1e-16))
    df_subset_near_end_sim_current_reg$Ne_S_Taj <- log10(df_subset_near_end_sim_current_reg$Ne_S_Taj+(1e-16))
    tryCatch(expr ={lst_current_sim_Ne_k_hat_vs_Gene_mobility <- ggplotRegression(fit=lmp(Ne_k_hat ~ nb_uniq_species_of_gene, data = df_subset_near_end_sim_current_reg,na.action=na.omit,center=FALSE),ggsave_path=results_output_workspace,the_filename=paste0("Multiple_gene_near_endsim_Ne_k_hat_vs_Gene_Mobility_sim_",simulation_name,".png"),xlabl = "log10(Gene_Mobility+(1e-16))",ylabl = "log10(Ne_k_hat+(1e-16))")},error=function(e) print(paste0("Regression skipped because of following error : ",e)))
    tryCatch(expr ={lst_current_sim_Ne_S_Taj_vs_Gene_mobility <-ggplotRegression(fit=lmp(Ne_S_Taj ~ nb_uniq_species_of_gene, data = df_subset_near_end_sim_current_reg,na.action=na.omit,center=FALSE),ggsave_path=results_output_workspace,the_filename=paste0("Multiple_gene_near_endsim_Ne_S_Taj_vs_Gene_Mobility_sim_",simulation_name,".png"),xlabl = "log10(Gene_Mobility+(1e-16))",ylabl = "log10(Ne_S_Taj+(1e-16))")},error=function(e) print(paste0("Regression skipped because of following error : ",e)))
    mtx_slopes_Ne_k_hat_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_k_hat_vs_Gene_mobility$slope_current_lm
    mtx_rsq_Ne_k_hat_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_k_hat_vs_Gene_mobility$adj_r_sq_current_lm
    mtx_pval_Ne_k_hat_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_k_hat_vs_Gene_mobility$p_val_current_lm
    mtx_slopes_Ne_S_Taj_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_S_Taj_vs_Gene_mobility$slope_current_lm
    mtx_rsq_Ne_S_Taj_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_S_Taj_vs_Gene_mobility$adj_r_sq_current_lm
    mtx_pval_Ne_S_Taj_vs_Gene_Mobility[i,j] <- lst_current_sim_Ne_S_Taj_vs_Gene_mobility$p_val_current_lm
    
  }
  
}


#plot heatmap Ne_k_hat vs Gene Mobility (WITHOUT FDR FILTER)
png(filename = paste0(results_output_workspace,"Heatmap_Ne_k_hat_vs_Gene_Mobility_WITHOUT_FDR_FILTER.png"),width = 1000, height = 700, units = "px")
heatmap.2(x = mtx_slopes_Ne_k_hat_vs_Gene_Mobility, dendrogram = "none", Rowv = FALSE, Colv = FALSE,key = TRUE,trace="none",scale="none", col = colorRampPalette(c("red","blue"), space = "rgb")(2))
dev.off()

#plot heatmap Ne_k_hat vs Gene Mobility (FDR FILTERED)
mtx_slopes_Ne_k_hat_vs_Gene_Mobility_fdr_filtered <- mtx_slopes_Ne_k_hat_vs_Gene_Mobility
mtx_slopes_Ne_k_hat_vs_Gene_Mobility_fdr_filtered[matrix(p.adjust(mtx_pval_Ne_k_hat_vs_Gene_Mobility,method="fdr"),nrow=length(hgt_rates))>0.05] <- NA
png(filename = paste0(results_output_workspace,"Heatmap_Ne_k_hat_vs_Gene_Mobility_FDR_FILTERED.png"),width = 1000, height = 700, units = "px")
heatmap.2(x = mtx_slopes_Ne_k_hat_vs_Gene_Mobility_fdr_filtered, dendrogram = "none", Rowv = FALSE, Colv = FALSE,na.color = "black",key = TRUE,trace="none",scale="none", col = colorRampPalette(c("blue","red"), space = "rgb")(2))
dev.off()

#plot heatmap Ne_S_Taj vs Gene Mobility (WITHOUT FDR FILTER)
png(filename = paste0(results_output_workspace,"Heatmap_Ne_S_Taj_vs_Gene_Mobility_WITHOUT_FDR_FILTER.png"),width = 1000, height = 700, units = "px")
heatmap.2(x = mtx_slopes_Ne_S_Taj_vs_Gene_Mobility, dendrogram = "none", Rowv = FALSE, Colv = FALSE,key = TRUE,trace="none",scale="none", col = colorRampPalette(c("red","blue"), space = "rgb")(2))
dev.off()

#plot heatmap Ne_S_Taj vs Gene Mobility (FDR FILTERED)
mtx_slopes_Ne_S_Taj_vs_Gene_Mobility_fdr_filtered <- mtx_slopes_Ne_S_Taj_vs_Gene_Mobility
mtx_slopes_Ne_S_Taj_vs_Gene_Mobility_fdr_filtered[matrix(p.adjust(mtx_pval_Ne_S_Taj_vs_Gene_Mobility,method="fdr"),nrow=length(hgt_rates))>0.05] <- NA
png(filename = paste0(results_output_workspace,"Heatmap_Ne_S_Taj_vs_Gene_Mobility_FDR_FILTERED.png"),width = 1000, height = 700, units = "px")
heatmap.2(x = mtx_slopes_Ne_S_Taj_vs_Gene_Mobility_fdr_filtered, dendrogram = "none", Rowv = FALSE, Colv = FALSE,na.color = "black",key = TRUE,trace="none",scale="none", col = colorRampPalette(c("blue","red"), space = "rgb")(2))
dev.off()

########
#Save RSession
save.session(paste0(results_output_workspace,"Nucleotide_Diversity_vs_Gene_Mobility_RSession.Rda"))


# #Restore Session
# library("session")
# sodapop_workspace <- "D:/Bureau/SodaPop_test_new/" #ABSOLUTE Path of the folder containing SodaPop files (should be the parent of directory of SodaPop "out/" or "src/")
# sodapop_workspace <- ifelse(test=substr(x = sodapop_workspace,start = nchar(sodapop_workspace),stop = nchar(sodapop_workspace))=="/",yes=sodapop_workspace,no=paste0(sodapop_workspace,"/"))
# sim_output_workpace <- paste0(sodapop_workspace,"out/") #ABSOLUTE path of the folder containing Sodapop simulations raw output repertory 
# results_output_workspace <- paste0(sodapop_workspace,"Sim_Results/") #ABSOLUTE path of the folder containing the Post-simulation genes variants analysis results. Will be created if it does not exist yet.
# restore.session(paste0(results_output_workspace,"Nucleotide_Diversity_vs_Gene_Mobility_RSession.Rda"))

