cd ./src/

# ----------------------------------------------------------
# Compile Programs - StageI
# ----------------------------------------------------------
gcc -o ../bin/SubsampleFastq SubsampleFastq.c
g++ -O3 -pthread -lm PrimerIdentification.c ../lib/edlib/edlib.cpp -o ../bin/PrimerIdentification
nvcc -O3 -arch=sm_75 -Xcompiler -fopenmp SlidingCorrelation.cu -o ../bin/SlidingCorrelation -lcudart
gcc -o ../bin/FilterThreshold FilterThreshold.c
gcc -O3 -fopenmp MultiReadMerging_StageI.c -o ../bin/MultiReadMerging_StageI
gcc -o ../bin/PostDecHammingDistance PostDecHammingDistance.c
gcc -o ../bin/ImageParsing ImageParsing.c

# ----------------------------------------------------------
# Compile Programs - StageII
# ----------------------------------------------------------
gcc -o ../bin/MajorityVoting MajorityVoting.c
gcc -o ../bin/ClusterSplitter ClusterSplitter.c
gcc -O3 -fopenmp MultiReadMerging_StageII.c -o ../bin/MultiReadMerging_StageII