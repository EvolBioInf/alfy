for a in 2 5 10 20; do
    echo -n $a " "
    ms $a 1 -t 10000 -r 0 5000000 |
	ms2dna > t.fasta
    getSeq "S1$" t.fasta > q.fasta
    getSeq -c "S1$" t.fasta > s.fasta
    /usr/bin/time -f "%U %E %M" ./alfy -i q.fasta -j s.fasta > /dev/null
done
