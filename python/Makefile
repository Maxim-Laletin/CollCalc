# Makefile for CollCalc implementation to Python

# Compiler
CC = g++

INCLUDE = -I ../include -I.
# CollCalc_to_python.h is located in the current dir

RM = rm

# GSL libraries
GSL_LIBS = -lgsl -lgslcblas

# Optimization flags
OPTIM = -Ofast 

CFLAGS = -shared -fPIC $(OPTIM) $(INCLUDE) $(PYTHON_INCLUDES)

# Source files
MAIN = CollCalc_to_python.cpp
SRC = ../source/Annihilation.cpp
PROCESS_DIR = ../processes
TARGET_DIR = py

MAIN_SRC = $(MAIN) $(SRC) $(PROCESS_DIR)/$(process).cpp

# SWIG settings
SWIG = swig
SWIG_FLAGS = -python -c++
PYTHON_INCLUDES = $(shell python3-config --includes)
PYTHON_LIBS = $(shell python3-config --ldflags)

# Ensure the target directory exists
$(TARGET_DIR):
	mkdir -p $(TARGET_DIR)

.PHONY: all clean

# The name is provided by the user
process ?= default

SWIG_INTERFACE = CollCalc_$(process).i
WRAP_CPP = CollCalc_wrap_$(process).cxx
PYTHON_MODULE = _collcalc_$(process).so

# Build the shared object for Python
to_python: $(TARGET_DIR)/$(PYTHON_MODULE) $(TARGET_DIR)/collcalc_$(process).py

# Generate the SWIG wrapper
$(WRAP_CPP): $(SWIG_INTERFACE)
	$(SWIG) $(SWIG_FLAGS) -o $@ $<

CollCalc_$(process).i: 
	@echo "Generating CollCalc_$(process).i..."
	@bash generate_interface.sh $(process)
	@echo "CollCalc_$(process).i created succesfully"

# Build the shared library
$(TARGET_DIR)/$(PYTHON_MODULE): $(WRAP_CPP) $(MAIN_SRC)
	$(CC) $(CFLAGS) $(WRAP_CPP) $(MAIN_SRC) $(PYTHON_LIBS) -o $@ $(GSL_LIBS)

# Move the Python file to the target directory
$(TARGET_DIR)/collcalc_$(process).py: $(WRAP_CPP)
	mv collcalc_$(process).py $(TARGET_DIR)

# Clean
clean:
	@ $(RM) -f $(TARGET_DIR)/* CollCalc_wrap_*.cxx CollCalc_*.py CollCalc_*.i
