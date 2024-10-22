# Makefile for the CollCalc program

# Compiler
CC = g++

INCLUDE = -I include 

RM = rm

GSL_LIBS = -lgsl -lgslcblas
#EXTRA_LIBS = -lc++ -framework Foundation
#GSL_LIBS = -L/usr/lib/x86-64_linux-gnu -lgsl -lgslcblas -lm

# Optimization flags
OPTIM = -Ofast 
#OPTIM = -O3

# Parallelization flags
MP = -fopenmp # requires OpenMP

CFLAGS = $(OPTIM) $(MP) # remove the MP flag if multithreading is not required

# Extra flags for the GSL package
#EXTRA_CFLAGS = -stdlib=libstdc++ -stdlib=libc++

# Source files
MAIN = main/CollCalc.cpp
SRC = source/Annihilation.cpp
# Directory with the process files
PROCESS_DIR = processes
# Directory to save the executables
TARGET_DIR = bin


# Ensure the target directory exists
$(TARGET_DIR):
	mkdir -p $(TARGET_DIR)

#all:

.PHONY: all clean

# CONDFLAG is the conditional compiler flag (see Collcalc.cpp)

# standard 2->2 annihilation (2 X states)
annihilation: CONDFLAG = -DINT3 # 3-dimensional integral
annihilation: $(TARGET_DIR)/$(process)_ann

# one X state
co-annihilation: CONDFLAG = -DINT4 # 4-dimensional integral
co-annihilation: $(TARGET_DIR)/$(process)_coann

$(TARGET_DIR)/%: $(MAIN) $(SRC) $(PROCESS_DIR)/$(process).cpp | $(TARGET_DIR)
	$(CC) $(INCLUDE) $^ $(CFLAGS) -o $@ $(GSL_LIBS) $(CONDFLAG) # -pg

# make annihilation process=file_name (without the extension)	

clean:
	@ $(RM) $(TARGET_DIR)/*
