MPI_CC=CC
ifeq ($(PE_ENV),PGI)
	MPI_FLAGS=-O3 -fast -mp=nonuma
endif

ifeq ($(PE_ENV),INTEL)
        MPI_FLAGS=-O3 -openmp -Minfo=all
endif

ifeq ($PE_ENV),CRAY)
        MPI_FLAGS=-O3 -h omp
endif

ifeq ($PE_ENV),GNU)
        MPI_FLAGS=-O3 -fopenmp
endif


SOURCES= wmBrick3D.cpp OpenChannel3D.cpp WallMountedBrick.cpp
OBJECTS=wmBrick3D.o OpenChannel3D.o WallMountedBrick.o workArounds.o

TARGET=WMBrick3D

%.o: %.cpp
	$(MPI_CC) $(MPI_FLAGS) -c $^

$(TARGET): $(OBJECTS)
	$(MPI_CC) $(MPI_FLAGS) -o $@ $^

clean:
	rm -f *.o $(TARGET) *~
