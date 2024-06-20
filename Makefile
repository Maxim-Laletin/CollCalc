# Makefile for the CollCalc program

# Compiler
CC = g++

# INCDIR = 
# LIBDIR = 

# Extra flags for the GSL package
EXTRA_CFLAGS=-stdlib=libstdc++ -stdlib=libc++
#INCLUDE = -I/usr/include/gsl

INCLUDE = -I include 

RM = rm

# Source files
SRC = main/CollCalc.cpp
OBJ = $(SRC:.cpp=.o)

# Directory to save the executables
TARGET_DIR = bin


#EXTRA_LIBS = -lc++ -framework Foundation
#GSL_LIBS = -L/usr/lib/x86-64_linux-gnu -lgsl -lgslcblas -lm
GSL_LIBS = -lgsl -lgslcblas


# Ensure the target directory exists
$(TARGET_DIR):
	mkdir -p $(TARGET_DIR)

#all:

# standard 2->2 annihilation
Annihilation: main/CollCalc.cpp source/Annihilation.cpp processes/Axion_coannihilation.cpp | $(TARGET_DIR)
	$(CC) $(INCLUDE) $^ -Ofast -fopenmp -o $@ ${GSL_LIBS} -DINT3 #-pg -fopenmp

# one DM state (like for axions)
Co-annihilation: main/CollCalc.cpp source/Annihilation.cpp processes/Axion_coannihilation.cpp | $(TARGET_DIR)
	$(CC) $(INCLUDE) $^ -Ofast -fopenmp -o $@ ${GSL_LIBS} -DINT4 #-pg -fopenmp

%.o : %.cpp
	$(CC) $(INCLUDE) -c $< -O3 -fopenmp -o $@ ${GSL_LIBS} #${INCLUDE}
	#-L${LIBDIR} ${WSTP_LIB} ${EXTRA_LIBS} -o $@ ${GSL_LIBS}

clean : 
	@ ${RM} -rf *.o #*tm.cpp C0 
