cd ../ && make -s > make.log&& cd bin && ./RNAcofold -d 0 -p --noPS --noGU --noconv --noClosingGU -P /data/ntc/Repository/ViennaRNA-2.6.4_make/misc/dna_mathews2004.par < fold_test.seq > cofold.out
 ./RNAcofold -d 0 -p --noPS --noGU --noconv --noClosingGU  < fold_test.seq >> cofold.out
# python cofold.py
rm -rf *ps