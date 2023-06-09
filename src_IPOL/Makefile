# C source code

# C++ source code
CXXSRC	= main.cpp \
    ./LibImages/LibImages.cpp \
    ./LibImages/LibImagesThreads.cpp \
    ./LibImages/mt19937ar.c \
    ./TWI/TWI.cpp \
    ./TWI/ConnectedComponent.cpp \
    ./LEM/LEM.cpp \
    ./LEM/LandscapeEvolution.cpp \
    ./LEM/WaterEvolution.cpp \
    ./LEM/CreepEvolution.cpp \
    ./LEM/WaterAndSedimentationEvolution.cpp \
    ./Utilities/Memory.cpp \
    ./Utilities/Parameters.cpp \
    ./Utilities/ProgressBar.cpp \
    ./Utilities/Utilities.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= lem_bin

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops -g -mavx
#
# Miguel: sin optims
#COPT	= -O0

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	= -lpng -ltiff -lm -lrt -lfftw3f -lfftw3f_threads -lpthread -ggdb

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)
