CC   = gcc-9
BIN  = FHT.exe
OBJ  = FHT.o lalgebra.o testingstuff.o 
FFTW_INCDIR = /usr/local/include
FFTW_LIBDIR = /usr/local/lib

LIBS =  -L${FFTW_LIBDIR} -lfftw3
INCS =  -I. -I${FFTW_INCDIR}
CFLAGS = -O3 ${INCS} -Wno-implicit-function-declaration

all: ${BIN}

$(BIN): $(OBJ)
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f ${BIN} ${OBJ}
