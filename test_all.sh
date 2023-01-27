#!/bin/sh
color="\e[0;32m"

echo "\n========================================================================="
echo "${color}[Compiling]\e[0m"
cd /pbc/scripts && make
echo "\n${color}[Testing]\e[0m"
echo "[Testing]" > ../output.txt
/pbc/test /pbc/param/a.param > temp.txt
cat temp.txt && cat temp.txt >> ../output.txt
echo "\n${color}[CPABE_OPTION]\e[0m"
echo "\n[CPABE_OPTION]" >> ../output.txt
/pbc/cpabe_option /pbc/param/a.param > temp.txt
cat temp.txt && cat temp.txt >> ../output.txt
echo "\n${color}[DABE]\e[0m"
echo "\n[DABE]" >> ../output.txt
/pbc/dabe /pbc/param/a.param > temp.txt
cat temp.txt && cat temp.txt >> ../output.txt
echo "\n${color}[DABE_OPTION]\e[0m"
echo "\n[DABE_OPTION]" >> ../output.txt
/pbc/dabe_option /pbc/param/a.param > temp.txt
cat temp.txt && cat temp.txt >> ../output.txt
echo "\n${color}[Ending]\e[0m"
make clean
rm -f temp.txt
echo "=========================================================================\n"
