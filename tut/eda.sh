cp ../data/A+DQ083238.fasta .
cp ../data/hiv42.fasta .
alfy -i A+DQ083238.fasta -j hiv42.fasta
alfy -i A+DQ083238.fasta -j hiv42.fasta  | 
    awk -f quantifyGenotypes.awk
alfy -i A+DQ083238.fasta -j hiv42.fasta  | 
    sed 's/+.*//' | 
    awk -f quantifyGenotypes.awk
alfy -w 400 -i A+DQ083238.fasta -j hiv42.fasta  | 
    sed 's/+.*//' | 
    awk -f quantifyGenotypes.awk
alfy -i A+DQ083238.fasta -j hiv42.fasta -M
alfy -i A+DQ083238.fasta -j hiv42.fasta -M |
    sed 's/+.*//' | 
    awk -f quantifyGenotypes.awk
cutSeq -r 6556-6705 A+DQ083238.fasta > nh.fasta
blastn -query nh.fasta -subject hiv42.fasta -outfmt 6  
alfy -f 150 -i A+DQ083238.fasta -j hiv42.fasta -M |
    sed 's/+.*//' | 
    awk -f quantifyGenotypes.awk
alfy -i A+DQ083238.fasta -j hiv42.fasta -M -P 0.05 | 
    sed 's/+.*//' | 
    awk -f quantifyGenotypes.awk
