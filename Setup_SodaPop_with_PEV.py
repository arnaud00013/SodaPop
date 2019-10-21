#!/usr/bin/python3.6

import sys
import os

"""
@Author = Arnaud NG

-This script creates the input files for Sodapop simulations using the pangenome evolution module

-Arguments : The name of the SodaPop workspace, the name of the .gene header file, the name of the .cell header file,
the name of the population_data header file, the number of cells simulated in each species, the number of species
simulated and the INITIAL number of genes in each species cells.

-All header files should be in a folder called headers/ stored in the files/ repertory of SodaPop!!!
"""

#Save script arguments values
sodapop_workspace = str(sys.argv[1])
gene_header_file = str(sys.argv[2])
cell_header_file = str(sys.argv[3])
pop_data_header_file = str(sys.argv[4])
nb_cells_per_species = int(sys.argv[5])
nb_species_simulated = int(sys.argv[6])
init_nb_genes_per_cells = int(sys.argv[7])

#Initialize the gene list file
os.system("echo 'Gene_Count\t{1}' > {0}files/genes/gene_list.dat".format(sodapop_workspace,init_nb_genes_per_cells))

for current_gene_id in range (0,init_nb_genes_per_cells):
    #initialize file from header
    os.system("/bin/cp -f {0}files/header/{2} {0}files/genes/{1}.gene".format(sodapop_workspace,current_gene_id,gene_header_file))
    #replace Gene_NUM by appropriate gene id
    os.system("sed -i.bu '1s/.*/Gene_NUM    {1}/' {0}files/genes/{1}.gene".format(sodapop_workspace,current_gene_id))
    #Delete macOs sed backup file
    os.system("rm -f {0}files/genes/{1}.gene.bu".format(sodapop_workspace,current_gene_id))
    #Add file in gene list file
    os.system("echo 'G\t{1}.gene' >> {0}files/genes/gene_list.dat".format(sodapop_workspace, current_gene_id))

#Initialize the population data file
os.system("/bin/cp -f {0}files/header/{1} {0}files/start/population.dat".format(sodapop_workspace,pop_data_header_file))
for current_sp_id in range (1,nb_species_simulated+1):
    #initialize file from header
    os.system("/bin/cp -f {0}files/header/{2} {0}files/start/{1}.cell".format(sodapop_workspace, current_sp_id,cell_header_file))
    #Replace org_id line with appropriate value
    os.system("sed -i.bu '2s/.*/org_id  {1}/' {0}files/start/{1}.cell".format(sodapop_workspace, current_sp_id))
    #Delete macOs sed backup file
    os.system("rm -f {0}files/start/{1}.cell.bu".format(sodapop_workspace, current_sp_id))
    #Add initial genes list in current cell file
    for current_gene_id in range(0, init_nb_genes_per_cells):
        # Add file in gene list file
        os.system("echo 'G\t{1}' >> {0}files/start/{2}.cell".format(sodapop_workspace, current_gene_id,current_sp_id))

    #Add file in the population data file
    os.system("echo 'C\t{2}\tfiles/start/{1}.cell' >> {0}files/start/population.dat".format(sodapop_workspace, current_sp_id,nb_cells_per_species))

print("######################## END OF .cell FILES AUTOMATIC SETUP ################")
print("Don't forget to add the species codon usage data in the .cell files before executing Sodasumm if you want SodaPop to consider codon usage bias fitness effect. You can use the script add_Codon_Usage_data_into_Sodapop_cell_file_from_Kasuza_db_like_file.sh to do that (one species at the time for each species you simulate!). The script takes the path of a Kazusa database-like file first and the path of the corresponding .cell file. If you don't need this feature, you can execute Sodasumm now!!!")
print("############################################################################")
