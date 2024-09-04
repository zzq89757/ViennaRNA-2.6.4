# cd ../ && make -s > make.log&& cd bin && ./RNAcofold -d 2 -p --noPS --noGU --noconv --noClosingGU -P /data/ntc/Repository/ViennaRNA-2.6.4_make/misc/dna_mathews2004.par < fold_test.seq > cofold.out
#  ./RNAcofold -d 2 -p --noPS --noGU --noconv --noClosingGU  < fold_test.seq >> cofold.out
 ./RNAcofold -a -d 0 -p --noPS --noGU --noconv --noClosingGU  --noLP < fold_test.seq > cofold.out
 ./RNAcofold -a -d 2 -p --noPS --noGU --noconv --noClosingGU  -P /data/ntc/Repository/ViennaRNA-2.6.4_make/misc/dna_mathews2004.par< fold_test.seq >> cofold.out
# /data/ntc/Project/VB/ViennaRNA-2.6.4/src/bin/RNAcofold -a -d 0  --noPS --noGU  --noClosingGU  --noLP -P /data/ntc/Repository/ViennaRNA-2.6.4_make/misc/dna_mathews2004.par  --noconv< fold_test.seq >> cofold.out
/data/ntc/Repository/ViennaRNA-2.6.4/src/bin/RNAcofold -a -d 0  --noPS --noGU  --noClosingGU  --noLP -P /data/ntc/Repository/ViennaRNA-2.6.4_make/misc/dna_mathews2004.par  --noconv< fold_test.seq >> cofold.out
# python cofold.py
rm -rf *ps