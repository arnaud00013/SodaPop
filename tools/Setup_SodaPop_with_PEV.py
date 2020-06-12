#!/usr/bin/python3.6

import sys
import os
from random import seed
from random import sample
# seed random number generator
seed(12345)

"""
@Author = Arnaud NG

-This script creates the input files for Sodapop simulations using the pangenome evolution module

-Arguments : The absolute path of the SodaPop workspace, the name of the .gene header file, the name of the .cell header file, the name of the population_data header file,the number of species simulated, the size of the whole community pool of accessory genes, the number of core genes at the start of the simulation and the number of accessory genes that are mobile.

-All header files should be in a folder called headers/ stored in the files/ repertory of SodaPop!!!
"""

print("######################## SIMULATION SETUP SCRIPT ################")
print("@Author = Arnaud NG \n -This script creates the input files for Sodapop simulations using the pangenome evolution module \n -Arguments : The absolute path of the SodaPop workspace, the name of the .gene header file, the name of the .cell header file, the name of the population_data header file, the number of species simulated, the size of the whole community pool of accessory genes, the number of housekeeping/non-consitutive genes and the number of accessory genes that are mobile.")
init_inp=input("-All header files should be in a folder called headers/ stored in the files/ repertory of SodaPop!!! Press Enter to continue")

#Save script arguments values
sodapop_workspace = str(sys.argv[1])
if (sodapop_workspace[-1]!="/"): #make sure path finish by "/"
    sodapop_workspace = sodapop_workspace + "/"
gene_header_file = str(sys.argv[2])
cell_header_file = str(sys.argv[3])
pop_data_header_file = str(sys.argv[4])
nb_species_simulated = int(sys.argv[5])
nb_unique_community_acc_genes = int(sys.argv[6]) 
nb_housekeeping_genes = int(sys.argv[7])
nb_mobile_acc_genes = int(sys.argv[8])

#initialization
lst_housekeeping_genes = list(range(0,nb_housekeeping_genes)) #list of IDs for genes that are present in all species
lst_unique_community_acc_genes = list(range(nb_housekeeping_genes,nb_unique_community_acc_genes+nb_housekeeping_genes)) #list of IDs for genes that are not present in all species
lst_mobile_genes = sample(lst_unique_community_acc_genes, nb_mobile_acc_genes) #determine randomly which community accessory genes are mobile by respecting the number of mobile genes given by the user in input

#Make sure that the number of mobile genes in input is smaller than the initial number of genes
if (nb_mobile_acc_genes > nb_unique_community_acc_genes+nb_housekeeping_genes):
    raise Exception("Input Error : The number of community accessory genes that are mobile cannot exceed the number of accessory genes!")

#Initialize the gene list file
os.system("echo 'Gene_Count\t{1}' > {0}files/genes/gene_list.dat".format(sodapop_workspace,nb_unique_community_acc_genes+nb_housekeeping_genes))

for current_gene_id in range (0,nb_unique_community_acc_genes+nb_housekeeping_genes):
    #initialize file from header
    os.system("/bin/cp -f {0}files/header/{2} {0}files/genes/{1}.gene".format(sodapop_workspace,current_gene_id,gene_header_file))
    #replace Gene_NUM by appropriate gene id
    os.system("sed -i.bu '1s/.*/Gene_NUM    {1}/' {0}files/genes/{1}.gene".format(sodapop_workspace,current_gene_id))
    #Delete macOs sed backup file
    os.system("rm -f {0}files/genes/{1}.gene.bu".format(sodapop_workspace,current_gene_id))
    #Add file in gene list file
    os.system("echo 'G\t{1}.gene' >> {0}files/genes/gene_list.dat".format(sodapop_workspace, current_gene_id))
    if (current_gene_id in lst_mobile_genes): #identify mobile genes in .gene files
        os.system("echo Mobile >> {0}files/genes/{1}.gene".format(sodapop_workspace,current_gene_id))

#Initialize the population data file
os.system("/bin/cp -f {0}files/header/{1} {0}files/start/population.dat".format(sodapop_workspace,pop_data_header_file))
for current_sp_id in range (1,nb_species_simulated+1):
    #Ask user for the number of cells for this species
    nb_cells_for_current_species = input("How many cells for species {0} : ".format(current_sp_id))
    while ((not nb_cells_for_current_species.isdigit()) or ("." in nb_cells_for_current_species) or ("," in nb_cells_for_current_species) or ("-" in nb_cells_for_current_species)): #make sure input is an integer
        nb_cells_for_current_species = input("How many cells for species {0} (Please enter a positive integer!) : ".format(current_sp_id))
    nb_cells_for_current_species = int(nb_cells_for_current_species)

    #Ask the user for the number of genes for this species
    nb_genes_for_current_species = input("How many genes for species {0} : ".format(current_sp_id))
    while ((not nb_genes_for_current_species.isdigit()) or ("." in nb_genes_for_current_species) or ("," in nb_genes_for_current_species) or ("-" in nb_genes_for_current_species) or (int(nb_genes_for_current_species)<nb_housekeeping_genes)): #make sure input is a positive integer higher or equal to the number of housekeeping genes
        nb_genes_for_current_species = input("How many genes for species {0}? Please enter a positive integer higher or equal to the number of housekeeping genes ({1}) : ".format(current_sp_id,nb_housekeeping_genes))
    nb_genes_for_current_species = int(nb_genes_for_current_species)
    #initialize file from header
    os.system("/bin/cp -f {0}files/header/{2} {0}files/start/{1}.cell".format(sodapop_workspace, current_sp_id,cell_header_file))
    #Replace org_id line with appropriate value
    os.system("sed -i.bu '2s/.*/org_id  {1}/' {0}files/start/{1}.cell".format(sodapop_workspace, current_sp_id))
    #Delete macOs sed backup file
    os.system("rm -f {0}files/start/{1}.cell.bu".format(sodapop_workspace, current_sp_id))
    
    #Add genes list in current cell file
    for current_gene_id in (lst_housekeeping_genes + sample(lst_unique_community_acc_genes,nb_genes_for_current_species-nb_housekeeping_genes)):
        # Add file in gene list file
        os.system("echo 'G\t{1}' >> {0}files/start/{2}.cell".format(sodapop_workspace, current_gene_id,current_sp_id))

    #Add file in the population data file
    os.system("echo 'C\t{2}\tfiles/start/{1}.cell' >> {0}files/start/population.dat".format(sodapop_workspace, current_sp_id,nb_cells_for_current_species))


print("######################## END OF .cell FILES AUTOMATIC SETUP ################")
print("Don't forget to add the species codon usage data in the .cell files before executing Sodasumm if you want SodaPop to consider codon usage bias fitness effect. You can use the script add_Codon_Usage_data_into_Sodapop_cell_file_from_Kasuza_db_like_file.sh to do that (one species at the time for each species you simulate!). The script takes the path of a Kazusa database-like file first and the path of the corresponding .cell file. If you don't need this feature, you can execute Sodasumm now!!!")
print("############################################################################")
