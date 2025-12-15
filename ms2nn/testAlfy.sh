n=5
theta=100
rho=1
len=10000
ms $n 1 -t $theta -r $rho $len -T > ms.out
ms2nn -q 1 ms.out > ms.nn
ms2dna ms.out | tr -d S > ms.fasta
getSeq 1 ms.fasta > q.fasta
getSeq -c 1 ms.fasta > s.fasta
alfy -M -i q.fasta -j s.fasta > ms.alfy
paste ms.alfy ms.nn
