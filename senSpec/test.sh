awk -f senSpec.awk -v e=exp.txt -v o=obs.txt > get.txt
d=$(diff get.txt want.txt)
if [ "$d" != "" ]; then
    echo "Fail: $d"
else
    echo "Pass"
fi
rm get.txt
