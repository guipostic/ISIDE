for i in b c d e f g h i j l k m n o p; do grep ^A _block.pdb | grep -v H > test_.pdb; editconf -f test_.pdb -o _block.pdb ; done
for i in a b c d e f g h i j l k m n o p; do sed -i /'OXT'/d p_block.pdb ; done
