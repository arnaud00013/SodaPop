#!/bin/bash
#this script takes in input the path of a file in Kazusa database codon usage-like format (1st input argument), convert it and add the data in the input .cell file (2nd input argument). It also remove all stop codons data, which are useless for CAI calculations.

sed -e 's/([^()]*)/)/g' $1 | sed 's/)  /\n/g' | sed 's/)//g' | sed '/^$/d' | sed '/UAG/d' | sed '/UGA/d' | sed '/UAA/d' >> $2

echo "##################################################################"
echo "DON'T FORGET TO EXECUTE Sodasumm after Codon usage configurations are done for every species you want to simulate!!!"
echo "##################################################################"
