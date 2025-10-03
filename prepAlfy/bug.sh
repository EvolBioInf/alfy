#Create the output files for comparison
prepAlfy -q query -s subject -n > prepAlfy.out
alfy64 -i A+DQ083238.fasta -j hiv42.fasta -o alfy.out > /dev/null

#Extract the interval of interest
cutSeq -r 114-154 A+DQ083238.fasta > query.fasta

#Blast interval of interest to subject fasta file
blastn -query query.fasta -subject hiv42.fasta -outfmt 6 | head


