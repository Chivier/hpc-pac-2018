PCG_HOME=`pwd`
BIN_DIR=./bin
SRC_DIR=./
MATROOT= /usr/local/MATLAB/R2019a
MATINC_PATH=${MATROOT}/extern/include
MKLINC_PATH=${MKLROOT}/include
LIBRARY_PATH=${MATROOT}/bin/glnxa64

TARGET = ${BIN_DIR}/SOBI

CC = icc

OBJ = SOBI.o
CFLAGS = -std=c99 -march=native -O2 -qopenmp 
LINKLINE = -Wl,--start-group ${MATROOT}/sys/os/glnxa64/libstdc++.so.6 ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -lmx -lmat -Wl,-rpath=${MATROOT}/bin/glnxa64 -Wl,-rpath-link=${MATROOT}/bin/glnxa64

${TARGET}:${OBJ}
	$(CC) ${OBJ} -L${LIBRARY_PATH} ${LINKLINE} -o ${TARGET}

%.o:${SRC_DIR}/%.c
	$(CC) ${CFLAGS} -I${MATINC_PATH} -I${MKLINC_PATH} -c $<


clean:
	@rm -f $(BIN_DIR)/SOBI *.o

