# README #

The C++ Version requires 

#include <assert.h>         // assert
#include <iostream>         // std::cout,std::endl
#include <memory>           // std::unique_ptr
#include <string>           // std::string
#include <sstream>          // std::istringstream
#include <stdexcept>        // std::exception, std::runtime_error
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <unordered_map>	// std::unorder_map
#include <unistd.h> 		// std::abs
#include <random> 			// std::generator
#include <fstream>			// std::ofstream
#include <omp.h>			// Hyperthreading commands
#include <tuple>
#include <time.h>
#include <ctime>
#include <random>
#include <cmath>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

Be aware of the compilation process as:

/usr/local/bin/g++ -Wall -I/usr/local/include -c V1_Main_Tumour_Evoulution.cpp --std=c++11
/usr/local/bin/g++ -L/usr/local/lib V1_Main_Tumour_Evoulution.o -fopenmp -lgsl -lgslcblas -lm
rm V1_Main_Tumour_Evoulution.o
mv a.out V1_Main_Tumour_Evoulution
./V1_Main_Tumour_Evoulution



This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact