bash cross.sh -h
bash cross.sh
alfy -i q.fasta -j s.fasta
for a in $(seq 20); do
    bash cross.sh
    alfy -i q.fasta -j s.fasta |
          tail -n 2
done
bash neiSim.sh -h
bash neiSim.sh -s 1
alfy -i query.fasta -j subjects.fasta
ms2nn -q 1 haplotypes.ms
blastn -query query.fasta -subject subjects.fasta \
         -outfmt 6 |
    column -t
bash neiSim.sh -s 2
alfy -i query.fasta -j subjects.fasta |
    tee alfy.out
awk -f quantifyGenotypes.awk alfy.out
ms2nn -q 1 haplotypes.ms
printf "Alfy\tBlast\n"
tail -n +2 alfy.out |
    while read start end score anno; do
          printf "$anno\t"
          cutSeq -r $start-$end query.fasta > f.fasta
          blastn -query f.fasta -subject \
                 subjects.fasta -outfmt 6 |
              head -n 1 |
              cut -f 2
    done
bash neiSim.sh -r 1 -s 3
ms2nn -q 1 haplotypes.ms |
    tee exp.txt
alfy -i query.fasta -j subjects.fasta |
    tee obs.txt
awk -f accuracy.awk -v e=exp.txt -v o=obs.txt
for a in $(seq 100); do
    bash testAlfy.sh
done
