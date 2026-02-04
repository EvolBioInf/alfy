bash ./neiSim.sh -s 3
cat query.fasta subjects.fasta haplotypes.ms > get.txt
d=$(diff get.txt want.txt)
if [ "$d" != "" ]; then
    echo "Fail: $d"
else
    echo "Pass"
fi
rm get.txt
