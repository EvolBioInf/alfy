./alfy -i ../../data/A+DQ083238.fasta -j ../../data/hiv42.fasta > test.txt
d=$(diff test.txt r.txt)
if [[ d -ne "" ]]; then
    echo "alfy: failed"
else
    echo "alfy: passed"
fi
rm test.txt
