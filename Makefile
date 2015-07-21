MPI_CC=CC
ifeq ($(PE_ENV),PGI)
	MPI_FLAGS=-O3 -fast -acc -Minfo=acc -Mnoopenmp
	ifeq ($(CRAYPAT_COMPILER_OPTIONS),1)
		MPI_FLAGS+= -DCRAYPAT
	endif
else
	MPI_FLAGS=-O3 -hnoomp -hacc -hlist=m
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
