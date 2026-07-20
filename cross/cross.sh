l=10000
m=0.1
M=0.01
s=$RANDOM
while getopts "hl:m:M:s:" opt; do
    case $opt in
        l) l=$OPTARG;;
        m) m=$OPTARG;;
        M) M=$OPTARG;;
        s) s=$OPTARG;;
        h)
            echo "Usage: bash cross.sh [-h] [option]..."
            echo "Simulate single midpoint crossover between"
            msg="two subjects in s.fasta and one query in q.fasta."
            echo "   $msg"
            echo "Example: bash cross.sh -s 3"
            msg="-l int\n\tsequence length (default $l)\n"
            printf "   $msg"
            msg="-m float\n\tmutation rate before crossover"
            msg="$msg (default $m)\n"
            printf "   $msg"
            msg="-M float\n\tmutation rate after crossover"
            msg="$msg (default $M)\n"
            printf "   $msg"
            msg="-s int\n\tseed (default random integer)\n"
            printf "   $msg"
            exit 0;;
        ?)
        echo "Usage: bash cross.sh [-h] [option]..."
        echo "Simulate single midpoint crossover between"
        msg="two subjects in s.fasta and one query in q.fasta."
        echo "   $msg"
        echo "Example: bash cross.sh -s 3"
        msg="-l int\n\tsequence length (default $l)\n"
        printf "   $msg"
        msg="-m float\n\tmutation rate before crossover"
        msg="$msg (default $m)\n"
        printf "   $msg"
        msg="-M float\n\tmutation rate after crossover"
        msg="$msg (default $M)\n"
        printf "   $msg"
        msg="-s int\n\tseed (default random integer)\n"
        printf "   $msg"
        exit 1;;
    esac
done
ranseq -l $l -s $s |
    sed 's/>.*/>s_1/' > s1.fasta
s=$(($s + 1))
mutator -m $m -s $s s1.fasta |
    sed 's/>.*/>s_2/' > s2.fasta
mid=$(($l / 2))
s=$(($s + 1))
cutSeq -r 1-$mid s1.fasta |
    mutator -m $M -s $s |
    sed 's/>.*/>q/' > q.fasta
s=$((s + 1))
cutSeq -r $(($mid + 1))-$l s2.fasta |
    mutator -m $M -s $s |
    tail -n +2 >> q.fasta
cat s[12].fasta > s.fasta
rm s[12].fasta
