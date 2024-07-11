cd ../ && make -s > make.log&& cd bin && ./RNAcofold -d 0 -p --noPS --noGU --noconv --noClosingGU < fold_test.seq > cofold.out
python cofold.py
rm -rf *ps