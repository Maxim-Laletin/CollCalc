#!/bin/bash
# Generates the interface file for SWIG

process=$1

cat << EOF > CollCalc_${process}.i
%module collcalc_${process} 

%{
#include "CollCalc_to_python.h"
%}

/* Expose the integrate function */
double CoAnnihilation(double x, double q, double rel_acc);

double Annihilation(double x, double q, double p, double rel_acc);
EOF