bash neiSim.sh -n 3 \
     -t 100 \
     -r 1 \
     -l 10000
alfy -i query.fasta -j subjects.fasta > obs.txt
ms2nn -q 1 haplotypes.ms > exp.txt
awk -f accuracy.awk -v e=exp.txt -v o=obs.txt
