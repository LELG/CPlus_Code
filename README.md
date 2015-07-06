# README #

The C++ Version requires 

include <assert.h>        
include <iostream>         
include <memory>           
include <string>           
include <sstream>          
include <stdexcept>        
include <stdlib.h>      
include <unordered_map>
include <unistd.h> 		
include <random> 			
include <fstream>			
include <omp.h>
include <tuple>
include <time.h>
include <ctime>
include <random>
include <cmath>
include <iomanip>
include <gsl/gsl_randist.h>
include <gsl/gsl_ieee_utils.h>
include <gsl/gsl_randist.h>
include <gsl/gsl_rng.h>
 
Be aware of the compilation process as:

/usr/local/bin/g++ -Wall -I/usr/local/include -c V1_Main_Tumour_Evoulution.cpp --std=c++11
/usr/local/bin/g++ -L/usr/local/lib V1_Main_Tumour_Evoulution.o -fopenmp -lgsl -lgslcblas -lm
rm V1_Main_Tumour_Evoulution.o
mv a.out V1_Main_Tumour_Evoulution
./V1_Main_Tumour_Evoulution

######### Notes

The garbage collector is coded as unique pointers.
These version do not require yet the HPC layer.
The models are described in the pdf

Do not hesitate  to contact me for questions !

