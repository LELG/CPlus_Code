/***************************************************************************************************

                       
    				The University of Melbourne & Peter MacCallum Cancer Centre

							Version 1
       
      
 Author:       Luis E Lara, ID: 657070, llara@student.unimelb.edu.au
Supervisor:    David Goode, PhD
 
g++ -c -I/$HOME/include ranlib_prb.cpp


	Compile as: 

	/usr/local/bin/g++ -Wall -I/usr/local/include -c V1_Main_Tumour_Evoulution.cpp --std=c++11
	/usr/local/bin/g++ -L/usr/local/lib V1_Main_Tumour_Evoulution.o -fopenmp -lgsl -lgslcblas -lm
	rm V1_Main_Tumour_Evoulution.o
	mv a.out V1_Main_Tumour_Evoulution
	./V1_Main_Tumour_Evoulution



****************************************************************************************************/

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

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/format.hpp>

// includes for cmd line parsing
#include <boost/program_options.hpp>
#include "parser.hpp"
namespace po = boost::program_options;

#include "Random.hpp"


/******************************* Constant Varaible Definitions *******************************/

#define NUM_THREADS 4
#define CPP_NO_COPYING_OF( Clazz )      \
    Clazz( Clazz const& );              \
    Clazz& operator=( Clazz const& )

/************ MODEL PARAMETERS ************/
const double P_Ben = 0.001;
const double P_Del = 0.001;
const double P_Inc = 0.5;
const double P_Dec = 0.5;
const double sel_pres = 0.01;
const double mut_rate = 0.0002002;//0.0002002
const double death_rate = 0.02;
const double quiescent_rate = 0.005;
const double prolif_rate = 0.03;
const unsigned long long Pop_Size =   7500000000; //1000000000;
const unsigned long long Pop_Size_p = 1000000000; //1000000000;
const unsigned int stop_simulation_when_plateaus = 1;
const double driver_quantile = 4.264891;//4.264891; // 99%
const double killer_quantile = -1.281552; // 10%
const double beneficial_quantile = 1.644854;// 95%
const double deleterious_quantile = -0.39;//35%
const unsigned int time_of_recovery_default = 24;
const bool Save_Clones = false;


const size_t K = 4;

/************ MODEL PARAMETERS Don't modify these ************/
const double  PS = (double) Pop_Size;
const double diff = prolif_rate - death_rate ;


const unsigned int dt = 3600;


std::random_device rd;
std::default_random_engine uniforms(rd());
std::default_random_engine new_uniforms(rd());
std::mt19937 binomials(rd());
std::mt19937 poissons(rd());
std::mt19937 Normals(rd());

gsl_rng *r_global ;



//std::default_random_engine generator;






/*********************** WARNING DEFINITIONS ********************************/ 
/*
	This section declares the throwing exeptions for the execution of the 
	code.
*/
/****************************************************************************/
namespace cpp {
/****************************************************************************/

    using namespace std;

    /*** Main Fnctions ***/
    bool hopefully( bool const c ) { return c; }
    bool throwX( string const& s ) { throw runtime_error( s ); }

    /**

		This Funtion is designed to throw an exeption when reading from
		command line.

		Input: Void
		Output: string containing the read from commnad line.

    *********************/
    string lineFromInput()
    {
        string result;
        getline( cin, result )
            || throwX( "lineFromInput: std::getline failed (EOF?)" );
        return result;
    }

    /**
			| --> Polymorphic to lineFromInput(void)
		This Funtion is designed to throw an exeption when reading from
		command line.

		Input: a constant reference
		Output: a string containing the read val form input

    *******************************************/
    string lineFromInput( string const& prompt )
    {
        cout << prompt;
        return lineFromInput();
    }

    /**
			| --> Polymorphic to lineFromInput(void)
		This Funtion is designed to throw an exeption when reading from
		command line.

		Input: a constant reference
		Output: an integer number from the read line input

    *******************************************/
    int intFromInput( string const& prompt )
    {
        istringstream   stream( lineFromInput( prompt ) );
        int             result;

        stream >> result
            || throwX( "intFromInput: input line was not a valid number spec" );
        return result;
    }
}  // namespace cpp


/***********************  CLONAL EXPANSION CORE ********************************/ 
/*
	This section declares all the Data Structures and Functions of the basic 
	clonal expansion model.

	Contains:
			(DS)	Clone + Cons + Decons
			  \---- get_Clone_DS(...)

			(DS)	Clonal_Expansion + Cons + Decons
			  \---- get_Clonnal_Expansion_DS(...)
			  
			(Fun)
			  |---- print_Clonal_Expansion_DS(...)
			  |---- carcinogenesis(...)
			  |---- print_Clone_DS(...)
			  |---- Division_Model(...)
			  |---- Waiting_until_next_division_round(...)
			  |---- Local_Clone_Expansion(...)
			  |---- Basic_Clonal_Expansion(...)
			  |---- Compute_Total_Population(...)
			  |---- map_Feedback(...)
			  \---- compute_Tumour_Evolution(...)
*/
/****************************************************************************/
namespace core {
/****************************************************************************/

    using namespace std;
    using namespace cpp;

    /* DS CLONE */
    struct Clone
    /***********/
    {
        CPP_NO_COPYING_OF( Clone );

        /**
        	NOTE: Defining as a const makes it read only
        	consider this to model only non - Stem Cells.

        ************************************************************/
    	double const P_Beneficial;				// This is for the probability of acquiring a beneficial mutation
		double const P_Deleterious;				// This  is for the probability of acquirind a deleterious mutation
		double const P_Increment_Mutation_Rate; // This controls the Muation Rate for some clones
		double const P_Decrement_Mutation_Rate;	// This contorls the Decrement of the mutation rate for some cloens
		unsigned int const Number_Dec_Time;


		
		double  Death_Rate;
		double  Division_Rate;
		double  Proliferation_Rate;
		unsigned int Time_to_Division; 	// This value is in hours
		unsigned int  Number_of_Memebers_to_Start_Heterogeneity;
		
		unsigned int Number_of_Deleterious_Muatations;
		unsigned int Number_of_Beneficial_Mutations;
		unsigned int Number_of_Passenger_Mutations; 
		unsigned int Number_of_Driver_Mutations_in_Cell_Division;
		unsigned int Generation_ID_Counter;

		

		bool clone_extinct;
		

		double  Mutation_Rate;
       
       	/**
			The Following variables can be modified at runtime.	

       	********************************************************/

       	double  Mutation_Effect;
       	unsigned long long int Number_of_Mutations;

    	unsigned int Time_of_Recovery;
		unsigned int Remaining_Time_in_Recovery_Phase; // Hours to wait until next division
		unsigned int Clone_Size; 					  // Number of cells in teh clone 
		unsigned int Clone_Fitness;

		bool Initiall_Expasion_Period;
		bool Clone_Ready_to_Compete;

		bool Clone_in_Recovery_Mode;
		bool Progenitor_Ready_to_Reproduce;
		bool Clone_Induction_of_Heterogeneity;

		unsigned int  Counter_of_Del_Time;
		
		//Mitosis values
		bool In_G0_Phase;
		bool In_G1_Phase;
		bool In_S_Phase;
		bool In_G2_Phase;
		bool In_M_Phase;

		double P_Expansion[3]; 


		unsigned int Remaining_Time_in_G1_Phase;
		unsigned int Remaining_Time_in_S_Phase;
		unsigned int Remaining_Time_in_G2_Phase;
		unsigned int Remaining_Time_in_M_Phase;

		string Generation_ID;

		/** 
			CONSTRUCTORS: 
						Clone(void) -- > default
						clone(P_Ben) --> if desired

		*****************************************************/
        Clone(po::variables_map params):
				P_Beneficial(params["prob_mut_pos"].as<double>()),
				P_Deleterious(params["prob_mut_neg"].as<double>()),
				P_Increment_Mutation_Rate(params["prob_inc_mut"].as<double>()),
				P_Decrement_Mutation_Rate(params["prob_dec_mut"].as<double>()),
				Number_Dec_Time(10),
				Death_Rate(params["die"].as<double>()),
				Division_Rate(params["pro"].as<double>()),
				Proliferation_Rate(params["pro"].as<double>()),
				Time_to_Division(24),
				Number_of_Memebers_to_Start_Heterogeneity(27),
				Number_of_Deleterious_Muatations(0),
				Number_of_Beneficial_Mutations(0),
				Number_of_Passenger_Mutations(0),
				Number_of_Driver_Mutations_in_Cell_Division(0),
				Generation_ID_Counter(0),
				clone_extinct(0),
				Mutation_Rate(params["mut"].as<double>()),
				Mutation_Effect(0.0),
				Number_of_Mutations(0),
				Time_of_Recovery(48),
				Remaining_Time_in_Recovery_Phase(0),
				Clone_Size(0),
				Clone_Fitness(0),
				Initiall_Expasion_Period(0),
				Clone_Ready_to_Compete(0),
				Clone_in_Recovery_Mode(0),
				Progenitor_Ready_to_Reproduce(0),
				Clone_Induction_of_Heterogeneity(0),
				Counter_of_Del_Time(0),
				In_G0_Phase(0),
				In_G1_Phase(0),
				In_S_Phase(0),
				In_G2_Phase(0),
				In_M_Phase(0),
				P_Expansion({Death_Rate, Proliferation_Rate, 1.0-(Death_Rate + Proliferation_Rate)}),
				Remaining_Time_in_G1_Phase(0),
				Remaining_Time_in_S_Phase(0),
				Remaining_Time_in_G2_Phase(0),
				Remaining_Time_in_M_Phase(0),
				Generation_ID("")
				{  							}

		~Clone()
		{}
    }; // end DS

    /**	FUNCTION get_Clone_DS()
		
		This function returns a 
		unique pointer of type Clone. Remmeber
		that unique pointers are memory efficient
		and guaratee a unique scope of the DS.

    ********************************************/
    unique_ptr<Clone> get_Clone_DS(po::variables_map params)
    {
        return unique_ptr<Clone>( new Clone(params) );
    } // end function

    /* DS Clonal_Expansion */
	struct Clonal_Expansion
	/**********************/
	{
		double Selective_Pressure;
		unsigned long long int Population_Size;
		//unsigned long long int Quiescent_Size;
		bool FD;
		
        po::variables_map params; // map to store all parameters of simulation
		double feedback; // This is the output of a proportional control system
		
		/**
			This is the main Data Structure to store all clones, this will
			reflect the heterogeneity.

			NOTE: This is a map of of unorder maps, this will allow us to model
			micro enviroments. The layer can be computed in MPI
		*/
		//unordered_map<string, unique_ptr<Clone> > Tumour;
		//unordered_map<int, struct Clone*> Tumour;
		/* We can create as many constructures as we want */

		 vector<unique_ptr<Clone> > *Tumour = new vector<unique_ptr<Clone> >;

		/* CONSTRUCTOR */
		Clonal_Expansion(po::variables_map params) :
				 				Selective_Pressure			( 0.0 ),
				 				Population_Size 			( 0   ),
			 					FD(0),
                                params                          ( params ),
								feedback 					( 0.0 )
			 					{ 									}
		~Clonal_Expansion() 
		{}

        void reduce_population(double sample_frac, Random r);
			 				
	};// end DS

	/**	FUNCTION get_Clonnal_Expansion_DS()
		
		This function returns a 
		unique pointer of type Clonal_Expansion. 

    ******************************************************/
	unique_ptr<Clonal_Expansion> get_Clonnal_Expansion_DS(po::variables_map params)
	{
		return unique_ptr<Clonal_Expansion>( new Clonal_Expansion(params) );
	}// end function


    /*
     * Reduce tumour population (in-place) by a fixed percentage.
     * Randomly select that percentage of cells from the population,
     * and discard the remainder.
     */
    void Clonal_Expansion::reduce_population(double sample_proportion, Random r)
    {
        if (sample_proportion <= 0 || 1 < sample_proportion)
        {
            throwX( "argument to reduce_population must be a proportion, a float in (0,1]" );
        }
        else if (sample_proportion == 1)
        {
            return;
        }

        int total_cells_to_sample = sample_proportion * Tumour->size();
        int cells_sampled_so_far = 0;

		vector<unique_ptr<Clone>> *tumour_sample = new vector<unique_ptr<Clone>>;

        for(std::vector<unique_ptr<Clone>>::iterator clone_it = Tumour->begin(); clone_it != Tumour->end(); ++clone_it) {
            int clone_sample_size = r.binomial_sample((*clone_it)->Clone_Size, sample_proportion);

            if (clone_sample_size > 0)
            {
                (*clone_it)->Clone_Size = clone_sample_size;
                tumour_sample->push_back(std::move(*clone_it));
                cells_sampled_so_far += clone_sample_size;
            }

            if (cells_sampled_so_far >= total_cells_to_sample)
            {
                break;
            }
        }

        Tumour = tumour_sample;
    }


	Random r;



	/************************** CORE FUNCTIONS ************************************/
	
	/**	FUNCTION print_Clonal_Expansion_DS(..) ------- # 1
		
		A basic Print of the DS.

		Input: unique constat pointer reference of type clonal Expansion.
		Ouput: Void

    ************************************************************************/
	void print_Clonal_Expansion_DS(unique_ptr<Clonal_Expansion> const  & CE) 
	{
		if(CE) // if NULL this will crash
		{
			cout << "\n\n##################################\n\n CLONAL EXPANSION DS \n{ " << endl;
			cout << "\tSelective_Pressure: " << CE -> Selective_Pressure << endl;
			cout << "\tPopulation_Size: " << CE -> Population_Size << endl;
			cout << "\tTumour size: " << CE -> Tumour -> size() << endl;
			cout << "}\n#################################" << endl;
		}
		else
			throwX( "WARNING: You send to print_Clonal_Expansion_DS() a NULL pointer, send a Clonal_Expansion pointer " );
	}// end function
	

	/**	FUNCTION carcinogenesis(..) ------- # 2
		
		This function introduces a new "progenitor cell into the DS". 
		This will seed a new tumour subpopulation with given 
		hallmarks acquired by the Clone DS.

		Input: (1) unique constat pointer reference of type clonal expansion
			   (2) int the layer, this is the key of the DS

		Ouput: Void

		NOTE: In further models the values may be change, and this function
		may be a good candidate for modications.  Specifically, in the 
		framework of Stem Cells. 

    ************************************************************************/
	void carcinogenesis(unique_ptr<Clonal_Expansion> const  & CE)
	{
		
		if(CE)// is the pointer NUL?
		{	
		
			CE -> Tumour -> push_back( get_Clone_DS(CE->params) );// We update the size that way
			CE -> Tumour -> back() -> Progenitor_Ready_to_Reproduce = false;
			CE -> Tumour -> back() -> Initiall_Expasion_Period = true;
			CE -> Tumour -> back() -> Clone_Size = 1;
			CE -> Tumour -> back() -> Mutation_Effect = 0.0;
			CE -> Tumour -> back() -> Mutation_Effect = 0.0;
			CE -> Tumour -> back() -> Time_of_Recovery = r.Recovery_After_Replication() + 48; // this should be a function
			
			CE -> Tumour -> back() -> In_G0_Phase = false;
			CE -> Tumour -> back() -> In_G1_Phase = true;
			CE -> Tumour -> back() -> In_S_Phase = false;
			CE -> Tumour -> back() -> In_G2_Phase = false;
			CE -> Tumour -> back() -> In_M_Phase = false;
			
			CE -> Tumour -> back() -> Time_to_Division = r.G1() + r.S() + 1 + r.G2();
			CE -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

			CE -> Tumour -> back() -> Clone_Fitness = 0;
			CE -> Tumour -> back() -> Clone_Ready_to_Compete = false;
			CE -> Tumour -> back() -> Clone_in_Recovery_Mode = true;
			CE -> Tumour -> back() -> Clone_Induction_of_Heterogeneity = false; 
			CE -> Tumour -> back() -> Generation_ID_Counter = 0; 
			CE -> Tumour -> back() -> Generation_ID = "P-0:0";

			CE -> Population_Size++;
			
		}
		else
			throwX( "WARNING: You send to carcinogenesis() a NULL pointer, send a Clonal_Expansion pointer " );

	}// end function
	


	


	/**	FUNCTION print_Clone_DS(..) ------- # 3
		
		A basic Print of the DS.

		Input: unique constat pointer reference of type Clone.
		Ouput: Void

    ************************************************************************/
	void print_Clone_DS(unique_ptr<Clone> const  & CL)
	{
		if(CL)// is the pointer valid?
		{
			cout << "\n\n##################################\n\n CLONE DS \n{ " << endl;
			cout << "\tP_Beneficial: " << CL -> P_Beneficial << endl;
			cout << "\tP_Deleterious: " << CL -> P_Deleterious << endl;
			cout << "\tP_Increment_Mutation_Rate: " << CL -> P_Increment_Mutation_Rate << endl;
			cout << "\tP_Decrement_Mutation_Rate: " << CL -> P_Decrement_Mutation_Rate << endl;
			cout << "\tMutatation_Rate: " << CL -> Mutation_Rate << endl;
			cout << "\tDeath_Rate: " << CL -> Death_Rate << endl;
			cout << "\tDivision_Rate: " << CL -> Division_Rate << endl;
			cout << "\tProliferation_Rate: " << CL -> Proliferation_Rate << endl;
			cout << "\tMutatation_Effect: " << CL -> Mutation_Effect << endl;
			cout << "\tNumber_of_Memebers_to_Start_Heterrogeneity: " << CL -> Number_of_Memebers_to_Start_Heterogeneity << endl;
			cout << "\tRemaining_Time_in_Recovery_Phase: " << CL -> Remaining_Time_in_Recovery_Phase << endl;
			cout << "\tClone_Size: " << CL -> Clone_Size << endl;
			cout << "\tClone_Fitness: " << CL -> Clone_Fitness << endl;
			cout << "\tInitiall_Expasion_Period: " << CL -> Initiall_Expasion_Period << endl;
			cout << "\tClone_Ready_to_Compete: " << CL -> Clone_Ready_to_Compete << endl;
			cout << "\tClone_in_Recovery_Mode: " << CL  -> Clone_in_Recovery_Mode << endl;
			cout << "\tProgenitor_Ready_to_Reproduce: " << CL -> Progenitor_Ready_to_Reproduce << endl;
			cout << "\tClone_Induction_of_Heterogeneity: " << CL -> Clone_Induction_of_Heterogeneity << endl;
			cout << "\tTime of Recovery: " << CL -> Time_of_Recovery << endl;
			cout << "\tTime to Division: " << CL -> Time_to_Division << endl;
			cout << "\tName: " << CL -> Generation_ID << endl;
			cout << "\tMuts at cell div: " << CL -> Number_of_Driver_Mutations_in_Cell_Division << endl;
			cout << "}\n#################################" << endl;
		}
		else
			throwX( "WARNING: You send to print_Clone_DS() a NULL pointer, send a Clone pointer " );
	}// function end

	/**	FUNCTION Division_Model(..) ------- # 4
		
		Basic equation for updating population size.
		This funtion is used only when the clone is in an initial stage of
		expansion

		Input: unsigned int of clone Size
		Ouput: unsigened long long int of the updated clone size

    ************************************************************************/
	unsigned long long int Division_Model(unsigned int Clone_Size) 
	{
		//Clone_Size +=  Clone_Size * 2;
		return (unsigned long long int ) (Clone_Size =  Clone_Size * 2);
	}

	/**	FUNCTION Basic_Clonal_Expansion(..) ------- # 7
		
		This model computes the basic model for a clonal expansion.
		(1) P(New) ~ Binom(Clone Size, p_new)
		(2) P(Die) ~ Binom(Clone Size, p_die)
       	(3) P(Ben) ~ Binom(New Clone, p_ben)
       	(4) P(Del) ~ Binom(New Clone, p_del)


		Input: (1) unique ptr of clonal_expasion --> main DS
			   (2) int layer, which key are we using
			   (3) int ith_clone, key of which clone are we handling.
			   (4) int hours, the actual hour of the model
		Ouput: int The number of beneficial mutations.

		NOTE: This version can be modified for different schemes.

    ****************************************************************************************************************/
	tuple<unsigned int,unsigned int, unsigned int, bool, bool> Basic_Clonal_Expansion(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int hours)
	{
		tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

		get<3>(Dying_and_Newborn) =  false;
		get<4>(Dying_and_Newborn) =  true;
		unsigned int n[3];
		unsigned int m = 0;
		
		double P_DR =  CE -> Tumour -> at(Generation_ID) -> P_Expansion[0];
		double P_NB =  CE -> Tumour -> at(Generation_ID) -> P_Expansion[1]  - CE -> feedback;
		double P_NT =  1.0 - (P_DR + P_NB );

		double p[] = {P_DR, P_NB, P_NT};


		
		
		gsl_ran_multinomial ( r_global, K, CE -> Tumour -> at (Generation_ID) -> Clone_Size, p, n);
		//cout << "DY: " << n[0] << " NB: " << n[1] << " NT: "<< n[2] << endl; 
	


			get<0>(Dying_and_Newborn) = n[0];
			// For Newborn Clones 
			get<1>(Dying_and_Newborn) = n[1];

			if(get<1>(Dying_and_Newborn) > 0)
			{	
				get<2>(Dying_and_Newborn) = r.Binomial_Mutants(
																get<1>(Dying_and_Newborn),
																CE -> Tumour -> at(Generation_ID) -> Mutation_Rate
			
															   );
				m = get<2>(Dying_and_Newborn);
				get<3>(Dying_and_Newborn) =  true;
			}

			
			int o = (int) ( CE -> Tumour -> at(Generation_ID) -> Clone_Size + n[1]) - (int) (n[0] + m);
			
			if(o <= 0)
			{
				
				CE -> Tumour -> at(Generation_ID) -> clone_extinct = true;

				//cout << "CS: " <<  CE -> Tumour -> at(Generation_ID) -> Clone_Size << " D: "  
				//	 << n[0]   << " B: " << n[1] << " M: " << m << " o: "  << o   << " ID extinct "  << 
				//	 CE -> Tumour -> at(Generation_ID) -> Generation_ID << " E: " << 
				//	 CE -> Tumour -> at(Generation_ID) -> clone_extinct  << " PS: " << CE -> Population_Size << endl;

				get<4>(Dying_and_Newborn) =  false;

				CE -> Population_Size -=   (unsigned long long int) CE -> Tumour  -> at(Generation_ID) -> Clone_Size ;
				CE -> Tumour -> at(Generation_ID) -> Clone_Size = 0;
				//cout << " PS: " << CE -> Tumour -> at(Generation_ID) -> Clone_Size << endl;
				
				//getchar();

			}

			
		//}//
		//cout << "Model " << " CS: " << CE -> Tumour -> at(Generation_ID) -> Clone_Size << " DR: "  << CE -> Tumour -> at(Generation_ID) -> Death_Rate << " PR: " << CE -> Tumour -> at(Generation_ID) -> Proliferation_Rate << " Dy: " << get<0>(Dying_and_Newborn) << " NB: " << get<1>(Dying_and_Newborn) << " Muts: " << get<2>(Dying_and_Newborn) <<endl;
		return Dying_and_Newborn;
	}//Function




	string Generate_Clone_Name_ID(int Parent_Generation_ID_Counter, string Parent_ID, unsigned int years, unsigned int hours)
	{
		string heredity_pattern = Parent_ID.substr( 0, Parent_ID.find("-") );
		string Clone_ID ("");
		//cout << "I GET " << heredity_pattern << endl;

		if(heredity_pattern.compare("P") == 0)
		{
			//cout << "I GET 1 " << heredity_pattern << endl;
			string Clone_ID_tag ("PC");
			Clone_ID_tag.append(to_string(Parent_Generation_ID_Counter));Clone_ID_tag.append("-");
			Clone_ID_tag.append(to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(to_string(hours));
			Clone_ID = Clone_ID_tag;
		}
		else
		{
			
			string Clone_ID_tag ("");
			Clone_ID_tag.append(heredity_pattern); Clone_ID_tag.append(","); Clone_ID_tag.append(to_string(Parent_Generation_ID_Counter));
			Clone_ID_tag.append("-");
			Clone_ID_tag.append(to_string(years));Clone_ID_tag.append(":");Clone_ID_tag.append(to_string(hours));
			Clone_ID = Clone_ID_tag;

		}
		return Clone_ID;
	}


	void carcinogenesis_from_Driver(unique_ptr<Clonal_Expansion> const  & CE,  int Generation_ID , double driver_effect, unsigned int years, unsigned int hours)
	{

		//cout << "HERE A " << CE -> Population_Size <<endl;
		if(CE)// is the pointer NUL?
		{	
	
			double mr = CE -> Tumour -> at(Generation_ID) -> Mutation_Rate;
			double pr = CE -> Tumour -> back() -> P_Expansion[1];

			CE -> Tumour -> at(Generation_ID) -> Number_of_Driver_Mutations_in_Cell_Division++;
			string Clone_Name = Generate_Clone_Name_ID(CE -> Tumour -> at(Generation_ID) -> Generation_ID_Counter, CE -> Tumour -> at(Generation_ID)  -> Generation_ID, years, hours);
			CE -> Tumour -> at(Generation_ID) -> Generation_ID_Counter++;


			CE -> Tumour -> push_back( get_Clone_DS(CE->params) );// We update the size that way
			CE -> Tumour -> back() -> Generation_ID = Clone_Name;
			CE -> Tumour -> back() -> Progenitor_Ready_to_Reproduce = false;
			CE -> Tumour -> back() -> Initiall_Expasion_Period = true;
			CE -> Tumour -> back() -> Clone_Size = 1;
			CE -> Tumour -> back() -> Mutation_Effect = driver_effect;
			CE -> Tumour -> back() -> Time_of_Recovery = r.Recovery_After_Replication() + 48; // this should be a function
			CE -> Tumour -> back() -> Time_to_Division = r.G1() + r.S() + 1 + r.G2();

			CE -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

			CE -> Tumour -> back() -> In_G0_Phase = false;
			CE -> Tumour -> back() -> In_G1_Phase = true;
			CE -> Tumour -> back() -> In_S_Phase = false;
			CE -> Tumour -> back() -> In_G2_Phase = false;
			CE -> Tumour -> back() -> In_M_Phase = false;

			CE -> Tumour -> back() -> Remaining_Time_in_Recovery_Phase = 0;
			CE -> Tumour -> back() -> Clone_Fitness = 0;
			CE -> Tumour -> back() -> Clone_Ready_to_Compete = false;
			CE -> Tumour -> back() -> Clone_in_Recovery_Mode = true;
			CE -> Tumour -> back() -> clone_extinct = false;
			CE -> Tumour -> back() -> Clone_Induction_of_Heterogeneity = false; 
			CE -> Tumour -> back() -> Generation_ID_Counter = 0; 


			CE -> Tumour -> back() -> Mutation_Rate = r.Uniform_Mutation_Rate_2(mr, mut_rate);

			
			CE -> Tumour -> back() -> P_Expansion[1] =  r.Update_Proliferation_Rate(pr);
			CE -> Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = 1;
			CE -> Population_Size++;
		

		}
		else
			throwX( "WARNING: You send to carcinogenesis() a NULL pointer, send a Clonal_Expansion pointer " );

	}// end function

	double map_Model( double Current_Population_Size )
	{
		//cout << "DR : "<< death_rate << " Pop: " << Pop_Size << " CV: " << Current_Population_Size << endl;
		//variable2 = min2+(max2-min2)*((variable1-min1)/(max1-min1))
		return 0.0 + (1.0 - 0.0) * ((1.0 - 0.0) / (Current_Population_Size - 0.0));	
	}
	double map_Feedback( double Current_Population_Size )
	{
		//cout << "DR : "<< death_rate << " Pop: " << Pop_Size << " CV: " << Current_Population_Size << endl;

		return 0.0 + (diff - 0.0) * ((Current_Population_Size - 0.0) / ((double) PS - 0.0));	
	}


	void Modify_Proliferation_Effect(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, bool Beneficial, double muatational_effect)
	{
		if(CE)
		{
			if(Beneficial)
			{
				double eff = map_Model( CE -> Tumour -> at(Generation_ID) -> Clone_Size  ) ;
				//cout << "eff " << eff  << " P_Exp (1*) " << CE -> Tumour -> at(Generation_ID) -> P_Expansion[1] << endl;
				CE -> Tumour -> at(Generation_ID) -> P_Expansion[1] += eff;
				//cout << " P_Exp(1) " << CE -> Tumour -> at(Generation_ID) -> P_Expansion[1] << endl;
				//cout << "PR N: " << CE -> Tumour -> at(Generation_ID) -> Proliferation_Rate << endl;
				//CE -> Tumour -> at(Generation_ID) -> Mutation_Rate -= val;
				if(CE -> Tumour -> at(Generation_ID) -> Time_of_Recovery >= 48)
				{
					CE -> Tumour -> at(Generation_ID) -> Time_of_Recovery--;
				}

			}
			else // is deletrious
			{

				double eff = map_Model( CE -> Tumour -> at(Generation_ID) -> Clone_Size  ) ;
		
				CE -> Tumour -> at(Generation_ID) -> P_Expansion[1] -= eff;

				if(CE -> Tumour -> at(Generation_ID) -> Counter_of_Del_Time == CE -> Tumour -> at(Generation_ID) -> Number_Dec_Time)
				{
					CE -> Tumour -> at(Generation_ID) -> Counter_of_Del_Time = 0;
					CE -> Tumour -> at(Generation_ID) -> Time_of_Recovery++;
				}
			}
		}
	}

	void Mutational_Effect_From_Mutants(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours  )
	{
		if(Number_of_Mutants > 0)
		{
			double muatational_effect = 0.0;
			
			for(unsigned int i = 0; i < Number_of_Mutants; i++)
				carcinogenesis_from_Driver(CE,  Generation_ID,  muatational_effect,  years,  hours);
				
		}//if
	}

	void compute_Mutations(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID ,tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn, unsigned int years, unsigned int hours )
	{
		if( get<3>(Dying_and_Newborn) )
		{
			Mutational_Effect_From_Mutants(CE, Generation_ID, get<2>(Dying_and_Newborn), years, hours );

		}
		unsigned int New_poulation_Size = (
											(		
												CE -> Tumour -> at(Generation_ID) -> Clone_Size 
											 	-
												get<0>(Dying_and_Newborn)
											) 
										    -	
											CE -> Tumour  -> at(Generation_ID) -> Number_of_Driver_Mutations_in_Cell_Division 
										) 
										+ 	get<1>(Dying_and_Newborn) ;

		CE -> Population_Size +=  ( (unsigned long long int) New_poulation_Size - (unsigned long long int) CE -> Tumour  -> at(Generation_ID) -> Clone_Size) ;
		CE -> Tumour -> at(Generation_ID) -> Clone_Size = (unsigned int) New_poulation_Size;
		CE -> Tumour -> at(Generation_ID) -> Number_of_Driver_Mutations_in_Cell_Division = 0;
	}

	void Initial_Expansion_Mitosis(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int hours, unsigned int years)
	{
		//(1) Go first through G1 Phase
		if (
				CE -> Tumour -> at (Generation_ID) -> In_G1_Phase 					&&
				CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase != 0
			)
		{
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase--;	

		}// G1 phase
		//(1) Reconfiguration G1 phase new parameters
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G1_Phase 					&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase  = false; 
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> at (Generation_ID) -> In_S_Phase = true;
		}// Finishig G1 Phase
		// (2) Going to S phase
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_S_Phase 			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase != 0
				)
		{
			CE -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase--;	
		}
		//(2.A) Exiting S phase
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_S_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_S_Phase  = false;
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> at (Generation_ID) -> In_G2_Phase = true;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase != 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase--;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  = false;
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = true;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_M_Phase 
				)
		{
			unsigned long long int cells_after_division = Division_Model( CE -> Tumour -> at(Generation_ID) -> Clone_Size );
			CE -> Population_Size =  (CE -> Population_Size - (unsigned long long int) CE -> Tumour -> at(Generation_ID) -> Clone_Size ) + cells_after_division;
			CE -> Tumour -> at(Generation_ID) -> Clone_Size = (unsigned int) cells_after_division;	

			if(CE -> Tumour -> at(Generation_ID) -> Clone_Size >= CE -> Tumour -> at(Generation_ID) -> Number_of_Memebers_to_Start_Heterogeneity)
			{
				CE -> Tumour -> at(Generation_ID) ->  Initiall_Expasion_Period = false;		
			}
			
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = false;
			//cout << "Going to be quiescent? " << endl;
			//cout << "If not going to G1 " << endl;
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase = true;
			
		}

	}//End

	void Mitosis(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int hours, unsigned int years )
	{

		if (
				CE -> Tumour -> at (Generation_ID) -> In_G1_Phase 					&&
				CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase != 0
			)
		{
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase--;
			//cout << "Initial Expansion in G1 Phase of clone: " << Generation_ID << endl;	

		}// G1 phase
		//(1) Reconfiguration G1 phase new parameters
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G1_Phase 					&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase  = false; 
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> at (Generation_ID) -> In_S_Phase = true;
		}// Finishig G1 Phase
		// (2) Going to S phase
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_S_Phase 			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase != 0
				)
		{
			CE -> Tumour -> at(Generation_ID) -> Remaining_Time_in_S_Phase--;	
		}
		//(2.A) Exiting S phase
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_S_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_S_Phase  = false;
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> at (Generation_ID) -> In_G2_Phase = true;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase != 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase--;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  			&&
					CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase == 0
				)
		{
			CE -> Tumour -> at (Generation_ID) -> In_G2_Phase  = false;
			CE -> Tumour -> at (Generation_ID) -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = true;
		}
		else if (
					CE -> Tumour -> at (Generation_ID) -> In_M_Phase 			
				)
		{

			
			tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn = Basic_Clonal_Expansion( CE, Generation_ID, hours ); 
			if(get<4>(Dying_and_Newborn))
			{
				compute_Mutations(CE, Generation_ID, Dying_and_Newborn, years, hours );
			}
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = false;
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase = true;
			
		}

	}
		
    /*
     * Helper function for writing some top-level sim results to a file.
     */
    void write_results_file(std::unique_ptr<Clonal_Expansion> const& CE, std::string fpath, double runtime_secs,
                            unsigned int years, unsigned int hours, unsigned int seconds)
    {
        std::ofstream results_fstream;
        results_fstream.open(fpath);

        std::string header;
        header = "pop_size\tnum_clones\truntime\telapsed_sim_years\telapsed_sim_hours\telapsed_sim_seconds\n";
        results_fstream << header;

        int rounded_runtime_secs = (int)runtime_secs;
        std::ostringstream runtime;
        int runtime_mins = rounded_runtime_secs/60;
        runtime_secs -= runtime_mins * 60;
        int runtime_hrs = runtime_mins/60;
        runtime_mins = runtime_mins % 60;
        runtime << boost::format("%1%:%|2$02|:%3$04.1f") % runtime_hrs % runtime_mins % runtime_secs;

        results_fstream << CE->Population_Size << "\t"
                        << CE->Tumour->size() << "\t"
                        << runtime.str() << "\t"
                        << years << "\t"
                        << hours << "\t"
                        << seconds << "\n";

        results_fstream.close();
    }

    /*
     * Helper function for writing clone statistics to file.
     */
    void write_stats_file(std::unique_ptr<Clonal_Expansion> const& CE, std::string fpath)
    {
        std::ofstream stats_fstream;
        stats_fstream.open(fpath);

        std::string header;
        header = "id\tClone_size\tProliferation_Rate\tMutation_Rate\tExtinct\tG_ID\n";
        stats_fstream << header;

        for(unsigned int i=0; i < CE->Tumour->size(); i++)
        {
            std::unique_ptr<Clone>& curr_clone = CE->Tumour->at(i);

            stats_fstream << i << "\t"
                          << curr_clone->Clone_Size << "\t"
                          << curr_clone->P_Expansion[1] << "\t"
                          << curr_clone->Mutation_Rate << "\t"
                          << curr_clone->clone_extinct << "\t"
                          << curr_clone->Generation_ID  << "\n" ;
        }

        stats_fstream.close();
    }

	/**	FUNCTION compute_Tumour_Evolution(..) ------- # 10
		
		This function is the CORE of the Algorithm.
		(1) This function generates a time simulator model
		... Inside that 
		... (2) Waiting_until_next_division_round(...)
		... (3) Local_Clone_Expansion(...)
		... (4) Basic_Clonal_Expansion(...)
		... (5) Feedback

		Input: Void
		Ouput: Void

    ***********************************************************/
	void compute_Tumour_Evolution(po::variables_map params)
	{
		using namespace std;
        clock_t begin = clock(); // start timing simulation


		unique_ptr<Clonal_Expansion> const CE = get_Clonnal_Expansion_DS(params);
		print_Clonal_Expansion_DS(CE);
		// move to pass reference

		/* Time Variables */
		unsigned int seconds = 0;
		unsigned int hours = 0;
		unsigned int years = 0;

		string path = "./V1_Clone_Data/";

  		unsigned int elapsed_hours = 0; 

        ofstream tumour_size_file;
        ofstream debugging;

        // construct paths for output files
        fs::path run_dir = params["run_dir"].as<string>();

        if (!fs::is_directory(run_dir)) {
            fs::create_directory(run_dir);
        }

        fs::path ts_path = run_dir / "tumour_size.txt";
        fs::path debugging_path = run_dir / "debugging.txt";

        tumour_size_file.open(ts_path.string());
        debugging.open(debugging_path.string());

  		//bool exit_flag = true;
  		unsigned int times_to_wait = 0;

		unsigned int ith_clone = 0;


		carcinogenesis(CE);
		

		print_Clone_DS( CE -> Tumour -> at(0) );


		while(times_to_wait < stop_simulation_when_plateaus )
		{
			seconds += dt;
			if(seconds == 3600)
			{
				seconds = 0; hours ++;

				unsigned int size = CE -> Tumour  -> size() ;

				
				for (ith_clone = 0; ith_clone < size ; ith_clone++)
				{
					//(1) Is my clone extinct?
					if( !(CE -> Tumour -> at (ith_clone) -> clone_extinct) )
					{
						// (1.a) if not, is my clone in free division?
						if(CE -> Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
						{
							Initial_Expansion_Mitosis(CE, ith_clone, hours, years);
						}
						else
						{
							Mitosis( CE, ith_clone, hours, years);
						}
					}
				}//for

				if(Save_Clones)
				{
					// Saving massive files
					if(hours % 100 == 0)
					{
						for (ith_clone = 0; ith_clone < size ; ith_clone++)
						{
							ofstream myfile;
							string f = path;
							f.append(CE -> Tumour -> at (ith_clone) -> Generation_ID);
							f.append(".txt");
							myfile.open( f, ios_base::app );
    						myfile << elapsed_hours + 1 << "\t" << CE -> Tumour -> at (ith_clone) -> Clone_Size << "\n";	
						}
					}
				}// Save clones
				
				CE -> feedback =  map_Feedback( (double) CE -> Population_Size );
			
				
				//CE -> feedback =  map_Feedback( (double) CE -> Population_Size  );
				//if(hours % 10 == 0)
				//{
  					debugging << " \n\n ACTIVE CELLS " <<  CE -> Population_Size  << " CLONES " << CE -> Tumour -> size() << "   H: " << hours  << " Y: " << years << " FD: " << CE -> feedback<< " PB: "<< prolif_rate - CE -> feedback << endl;
  					
  				//	getchar();
                    tumour_size_file << CE -> Population_Size << "\t" << elapsed_hours << "\n";
					elapsed_hours += 1;
					
					if( CE -> Population_Size >  Pop_Size_p )
					{
						debugging << " \n\n ACTIVE CELLS " <<  CE -> Population_Size  << " CLONES " << CE -> Tumour -> size() << "   H: " << hours  << " Y: " << years << " FD: " << CE -> feedback<< " PB: "<< prolif_rate - CE -> feedback << endl;
  						times_to_wait++;

  						break;
  					}
				//}// if print
				
			}//seconds
			if(hours == 8764)//8764
			{
				hours = 0; years ++;
			}
		}//while

		cout << "FINISH MAIN LOOP ..... SAVING stats" << endl;
		
        tumour_size_file.close();
        debugging.close();

        clock_t end = clock(); // finish timing simulation

        double cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
        fs::path results_path = run_dir / "results.txt";
        write_results_file(CE, results_path.string(), cpu_secs, years, hours, seconds);
		
        fs::path stats_path = run_dir / "stats.txt";
        write_stats_file(CE, stats_path.string());
	}// end of function

}  // namespace blah


int main(int argc, char *argv[])
{
	using namespace std;
	using namespace core;
	
	long seed;


  r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  seed = time (NULL) * getpid();    
  gsl_rng_set (r_global, seed);  

    try
    {
        po::variables_map params = parser::parse_options(argc, argv, cout);
        compute_Tumour_Evolution(params);
        
        return EXIT_SUCCESS;
    }
    catch( exception const& error )
    {
        cerr << "!" << error.what() << endl;
    }
    gsl_rng_free (r_global);

	return EXIT_SUCCESS;

}


