l=10000
s=$RANDOM
while getopts "hl:s:" opt; do
    case $opt in
          l) l=$OPTARG;;
          s) s=$OPTARG;;
          h)
              printf "Usage: bash cross.sh [-h] [option]...\n"
              printf "Simulate single midpoint crossover between one\n"
              printf "   query and two subjects in q.fasta and s.fasta.\n"
              printf "Example: bash cross.sh -s 3\n"
              printf "   -l int\n\tsequence length (default $l)\n"
              printf "   -s int\n\tseed (default random integer)\n"
              exit 0;;
          ?)
              printf "Usage: bash cross.sh [-h] [option]...\n"
              printf "Simulate single midpoint crossover between one\n"
              printf "   query and two subjects in q.fasta and s.fasta.\n"
              printf "Example: bash cross.sh -s 3\n"
              printf "   -l int\n\tsequence length (default $l)\n"
              printf "   -s int\n\tseed (default random integer)\n"
              exit 1;;
    esac
done
ranseq -l $l -s $s |
    sed 's/>.*/>s_1/' > s1.fasta
s=$(($s + 1))
mutator -m 0.1 -s $s s1.fasta |
    sed 's/>.*/>s_2/' > s2.fasta
s=$(($s + 1))
cutSeq -r 1-5000 s1.fasta |
    mutator -m 0.01 -s $s |
    sed 's/>.*/>q/' > q.fasta
s=$((s + 1))
cutSeq -r 5001-10000 s2.fasta |
    mutator -m 0.01 -s $s |
    tail -n +2 >> q.fasta
cat s[12].fasta > s.fasta
