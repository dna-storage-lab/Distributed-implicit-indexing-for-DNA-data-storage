cd ./src/
gcc -o ../bin/SubsampleFastq SubsampleFastq.c
g++ -O3 -pthread -lm PrimerIdentification.c ../lib/edlib/edlib.cpp -o ../bin/PrimerIdentification
nvcc -O3 -arch=sm_75 -Xcompiler -fopenmp SlidingCorrelation.cu -o ../bin/SlidingCorrelation -lcudart
gcc -o ../bin/MajorityVoting MajorityVoting.c
gcc -o ../bin/PostDecHammingDistance PostDecHammingDistance.c