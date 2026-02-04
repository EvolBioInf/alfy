bash ./cross.sh -s 3
cat q.fasta s.fasta > get.txt
d=$(diff get.txt want.txt)
if [ "$d" != "" ]; then
    echo "Fail: $d"
else
    echo "Pass"
fi
rm get.txt
