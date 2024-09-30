gcc -g Dsuite.c -o Dsuite -lm
./Dsuite > a.out
# echo "\n"
RNAcofold -p -d0 --noLP --noClosingGU --noPS  --noconv --noGU < /data/ntc/Repository/ViennaRNA-2.6.4_make/src/bin/fold_test.seq >> a.out
cd ../../../ && make -s > make.log && cd -
../../RNAcofold -p -d0 --noLP --noClosingGU --noconv --noGU --noPS < /data/ntc/Repository/ViennaRNA-2.6.4_make/src/bin/fold_test.seq >> a.out
/data/ntc/Repository/ViennaRNA-2.6.4/src/bin/RNAcofold -p -d0 --noLP --noClosingGU --noconv --noGU --noPS < /data/ntc/Repository/ViennaRNA-2.6.4_make/src/bin/fold_test.seq >> a.out
# RNAcofold -p -d0 --noLP --noClosingGU --noGU --noconv  --noPS < /data/ntc/Repository/ViennaRNA-2.6.4_make/src/bin/fold_test.seq
# RNAcofold -p -d0 --noLP --noClosingGU --noGU --noconv < /data/ntc/Repository/ViennaRNA-2.6.4_make/src/bin/fold_test.seq