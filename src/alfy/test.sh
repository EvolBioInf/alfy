./alfy -i ../../data/A+DQ083238.fasta -j ../../data/hiv42.fasta > test.txt
d=$(diff test.txt r.txt)
if [[ $d != "" ]]; then
    echo "alfy: failed: $d"
else
    echo "alfy: passed"
fi
rm test.txt
