MATROOT=/usr/local/MATLAB/R2019a

cd bin
mkl_num_threads=10 KMP_AFFINITY=compact,granularity=fine LD_PRELOAD=${MATROOT}/sys/os/glnxa64/libstdc++.so.6 ./SOBI
