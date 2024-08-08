trimal -in $1.reordered.fas -out $1.trimal.fas -gt 0.4
iqtree2 -s $1.trimal.fas -n 4 -m GTR+FO*H4
