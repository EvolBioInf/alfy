n=3
t=100
r=0
l=10000
s=$RANDOM
while getopts "hn:t:r:l:s:" opt; do
    case $opt in
          n) n=$OPTARG;;
          t) t=$OPTARG;;
          r) r=$OPTARG;;
          l) l=$OPTARG;;
          s) s=$OPTARG;;
          h)
            printf "Usage: bash neiSim.sh [-h] [option]...\n"
            printf "Neighbors simulation as one query in\n"
            printf "\tquery.fasta and one or more subjects\n"
            printf "\tin subjects.fasta; the underlying\n"
            printf "\thaplotypes are in haplotypes.ms.\n"
            printf "Example: bash neiSim.sh -n 5 -t 50 -r 1\n"
            printf "   -n int\n\tsample size (default $n)\n"
            printf "   -t float\n\ttheta (default $t)\n"
            printf "   -r float\n\trho (default $r)\n"
            printf "   -l int\n\tsequence length (default $l)\n"
            printf "   -s int\n\tseed (default random integer)\n"
            exit 0;;
          ?)
            printf "Usage: bash neiSim.sh [-h] [option]...\n"
            printf "Neighbors simulation as one query in\n"
            printf "\tquery.fasta and one or more subjects\n"
            printf "\tin subjects.fasta; the underlying\n"
            printf "\thaplotypes are in haplotypes.ms.\n"
            printf "Example: bash neiSim.sh -n 5 -t 50 -r 1\n"
            printf "   -n int\n\tsample size (default $n)\n"
            printf "   -t float\n\ttheta (default $t)\n"
            printf "   -r float\n\trho (default $r)\n"
            printf "   -l int\n\tsequence length (default $l)\n"
            printf "   -s int\n\tseed (default random integer)\n"
            exit 1;;
    esac
done
ms $n 1 -t $t -r $r $l -seeds 0 0 $s -T > haplotypes.ms
ms2dna -s $s haplotypes.ms |
    tr -d S > tmp.fasta
getSeq 1 tmp.fasta > query.fasta
getSeq -c 1 tmp.fasta > subjects.fasta
rm tmp.fasta
