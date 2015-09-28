/***************************************************************************************************

                       
    				The University of Melbourne & Peter MacCallum Cancer Centre

                    			Clonala Expansion Model Proposal
 
 Model:

			(1) P(New) ~ Binom(Clone Size, p_new)
			(2) P(Die) ~ Binom(Clone Size, p_die)
       		(3) P(Ben) ~ Binom(New Clone, p_ben)
       		(4) P(Del) ~ Binom(New Clone, p_del)
       
      
 Author:       Luis E Lara, ID: 657070, llara@student.unimelb.edu.au
Supervisor:    David Goode, PhD
 

	Compile as: 
				/usr/local/bin/g++ -Wall -I/usr/local/include -c TEM_V1.cpp --std=c++11
				/usr/local/bin/g++ -L/usr/local/lib TEM_V1.o -lgsl -lgslcblas -lm -lmpi_cxx -lmpi 
				rm TEM_V1.o
				mv a.out TEM_V1
	
	Run  by:
				mpirun -np 2 ./TEM_V1


****************************************************************************************************/
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
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
//#include <omp.h>			// Hyperthreading commands
#include <tuple>
 #include <ctime>
 #include <unistd.h>
#include <time.h>
#include <ctime>
#include <random>
#include <cmath>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iterator>
#include <utility>
#include <initializer_list>
#include <mpi.h>	// Needs to be included in order to use MPI



/******************************* Constant Varaible Definitions *******************************/

#define NUM_THREADS 4
#define CPP_NO_COPYING_OF( Clazz )      \
    Clazz( Clazz const& );              \
    Clazz& operator=( Clazz const& )

/************ MODEL PARAMETERS ************/

const double MUT_RATE = 0.0002002;//0.0002002
const double DEATH_RATE = 0.02;
const double PROLIFERATION_RATE = 0.03;

const unsigned long long MAXIMUM_POPULATION_SIZE =   3000000000; //1000000000;
const double  PS = (double) MAXIMUM_POPULATION_SIZE;

const unsigned long long DETECTABLE_POPULATION_SIZE =    400000; //1000000000;
const unsigned int STOP_GROWTH_AFTER_DIAGNOSIS = 10; //[hours]

const bool SAVE_CLONAL_EVOLUTION = false;

/* Main path of the model */
//const std::string  CLONAL_EVOLUTION_PATH = "./V1_Clone_Data/";
const std::string DEFAULT_TREATMENT_FILE = "./DRUG/V1_CTX_Scheme.drug";

const std::string VERSION = "/V1/";

const size_t K = 4;

const bool PRINT = true;

/************ MODEL PARAMETERS Don't modify these ************/

const double DIFF = PROLIFERATION_RATE - DEATH_RATE ;
const unsigned int dt = 3600;

const bool MULTIPLE_TUMOURS = false;
//const bool MULTIPLE_TREATMENTS = false;

/*
	Types of Penaltys
	1 - Mean
	2 - 25% Quartile
	3 - Static to 0.03
*/
const int PENALTY = 3;

//extern const std::string SIMULATION_PATH;

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

		unsigned int  Number_of_Memebers_to_Start_Heterogeneity;
		unsigned int Generation_ID_Counter;
		bool clone_extinct;
		double  Mutation_Rate;
       
       	/**
			The Following variables can be modified at runtime.	

       	********************************************************/

       //	double  Mutation_Effect;
       	unsigned long long int Number_of_Mutations;

		//unsigned int Remaining_Time_in_Recovery_Phase; // Hours to wait until next division
		unsigned int Clone_Size; 					  // Number of cells in teh clone 
		

		bool Initiall_Expasion_Period;
		
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
        Clone(void):
			
				Number_of_Memebers_to_Start_Heterogeneity(27),
				Generation_ID_Counter(0),
				clone_extinct(0),
				Mutation_Rate(MUT_RATE),
				Number_of_Mutations(0),
				Clone_Size(0),
				Initiall_Expasion_Period(0),
				In_G0_Phase(0),
				In_G1_Phase(0),
				In_S_Phase(0),
				In_G2_Phase(0),
				In_M_Phase(0),
				P_Expansion{DEATH_RATE, PROLIFERATION_RATE, 1.0-(DEATH_RATE + PROLIFERATION_RATE)},
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
    unique_ptr<Clone> get_Clone_DS()
    {
        return unique_ptr<Clone>( new Clone( ) );
    } // end function

    /* DS Clonal_Expansion */
	struct Clonal_Expansion
	/**********************/
	{
		
		unsigned long long int Population_Size;
		//unsigned long long int Quiescent_Size;
		
		
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
		Clonal_Expansion(void) :
				 				
				 				Population_Size 			( 0   ),
			 					feedback 					( 0.0 )
			 					{ 									}
		~Clonal_Expansion() 
		{}
			 				
	};// end DS

	/**	FUNCTION get_Clonnal_Expansion_DS()
		
		This function returns a 
		unique pointer of type Clonal_Expansion. 

    ******************************************************/
	unique_ptr<Clonal_Expansion> get_Clonnal_Expansion_DS()
	{
		return unique_ptr<Clonal_Expansion>( new Clonal_Expansion() );
	}// end function

	struct Drug
	{
		double maximum_tolareted_dose;
		double minimum_dose;
		double half_life;
		double plasma_volume;
		double dosage;
		double absorption_Fraction;
		unsigned int  interval;

		double elimintation_constant;

		double eliminated;
		double drug_in_system;
		double concentration;

		//double dx;

		unsigned int hour_counter;

		double ingestion;
		double bioavailavility;
		unsigned int cycles;
		unsigned int repetition;
		vector<unsigned int> *intake_days_in_hours = new vector<unsigned int>;
		unsigned int time_to_cell_death;
		unsigned int time_to_mutation;

		/*
			<vector> of drug concentration to save as delay
			unsigned int this is tghe index counter
			unsigned int this is the max index to reset 
			bool to print
			bool to add penalty
		*/

		tuple< vector<double> , unsigned int, unsigned int, bool> MR_Register;
		tuple< vector<double> , unsigned int, unsigned int, bool> DR_Register;

		bool done;
		//unsigned int final_counter;
		bool exit;
		string Family;
		string DrugName;
		bool Penalty_PR;
		bool Penalty_DR;
		bool Penalty_MR;


		Drug(void):
			maximum_tolareted_dose(0.0),
			minimum_dose(0.0),
			half_life(0.0),
			//plasma_volume(0.0),
			dosage(0.0),
			absorption_Fraction(0.0),
			interval(0),
			elimintation_constant(0.0),
			eliminated(0.0),
			drug_in_system(0.0),
			concentration(0.0),
			//dx(0.0),
			hour_counter(0),
			ingestion(0.0),
			bioavailavility(0.0),
			//cycles(0),
			//repetition(0),
			time_to_cell_death(0),
			time_to_mutation(0),
			MR_Register{{},0,0,0},
			DR_Register{{},0,0,0},
			done(false),
			//final_counter(1000),
			exit(false),
			Family(""),
			DrugName(""),
			Penalty_PR(false),
			Penalty_DR(false),
			Penalty_MR(false)
			{}
		
		~Drug()
			{}

	};

	unique_ptr<Drug> get_Drug()
	{
		return unique_ptr<Drug>( new Drug() );
	}

	/*
		Treatment DS
	*/
	struct Treatment
	{
		// Whcih type of shcme we need: Main, Adjuvant, Both, Metro, etc
		string Scheme;
		unsigned int Number_of_Drugs;
		string Type_of_Treatment;
		unsigned int Repetition_Day_in_Hours;
		unsigned int Cycles;
		double Maximum_Tolareted_Toxicity; //ug/ml
		double Toxicity;
		double Plasma_Volume;
		double Total_Drug_Ingestion;
		double Total_Drug_Eliminated;
		double Total_Drug_in_System;
		double Total_Drug_Concentration;
		double dx;
		unsigned int hours_counter;
		unsigned int hours_to_Exit;
		bool Trial_Done;
		vector<unique_ptr<Drug> > *AntiCancer_Agents = new vector<unique_ptr<Drug> >;

		/* CONSTRUCTOR */
		Treatment(void):
						Scheme (""),
						Number_of_Drugs(0),
						Type_of_Treatment(""),
						Repetition_Day_in_Hours(0),
						Cycles(0),
						Maximum_Tolareted_Toxicity(0.0),
						Toxicity(0.0),
						Plasma_Volume(0.0),
						Total_Drug_Ingestion(0.0),
						Total_Drug_Eliminated(0.0),
						Total_Drug_in_System(0.0),
						Total_Drug_Concentration(0.0),
						dx(0.0),
						hours_counter(0),
						hours_to_Exit(300),
						Trial_Done(0)
						

						{}
		~Treatment()
						{}
	};

	unique_ptr<Treatment> get_Treatment_DS()
	{
		return unique_ptr<Treatment>( new Treatment() );
	}



	class Random 
	{
		public:

    	Random() = default;
    	Random(std::mt19937::result_type seed) : eng(seed) {}
    	unsigned int Recovery_After_Replication();
    	unsigned int G1();
    	unsigned int S();
    	unsigned int G2();
    	unsigned int Poisson();
    	double Z();
    	double Binomial_dying(unsigned int Clone_Size, double Death_Rate);
    	double Binomial_newborn(unsigned int Clone_Size, double Adjusted_Proliferation_Rate);
    	unsigned int Binomial_Mutants(unsigned int NewBorn_Size, double Mutation_Rate);
    	//unsigned int* Mitosis_Multinomial(unsigned int Clone_Size);
    	//unsigned int Uniform_Mutation_Division_Memebers(unsigned int Clone_Size);
    	double Uniform_Mutation_Rate(double mu_rate);
    	double Update_Proliferation_Rate(double Proliferation_Rate);
    	double Uniform_Mutation_Rate_2(double Parent_mu_rate);


		private:        
    	std::mt19937 eng{std::random_device{}()};
    	

  
 		
	};


	unsigned int Random::Recovery_After_Replication()
	{
    	return uniform_int_distribution<unsigned int>{12, 15}(eng);
	}

	unsigned int Random::G1()
	{
		return uniform_int_distribution<unsigned int>{12, 15}(eng);
	}

	unsigned int Random::S()
	{
		return uniform_int_distribution<unsigned int>{5, 10}(eng);
	}

	unsigned int Random::G2()
	{
		return uniform_int_distribution<unsigned int>{3, 5}(eng);
	}

	unsigned int Random::Poisson()
	{
		return poisson_distribution<unsigned int>{1}(eng);
	}

	double Random::Z()
	{
		return normal_distribution<double>{0.0,1.0}(eng);
	}

	double Random::Binomial_dying(unsigned int Clone_Size, double Death_Rate)
	{
		return binomial_distribution<unsigned int>{Clone_Size, Death_Rate}(eng);
	}
	
	double Random::Binomial_newborn(unsigned int Clone_Size, double Adjusted_Proliferation_Rate)
	{
		return binomial_distribution<unsigned int>{Clone_Size, Adjusted_Proliferation_Rate }(eng);
	}

	unsigned int Random::Binomial_Mutants(unsigned int NewBorn_Size, double Mutation_Rate)
	{
		return binomial_distribution<unsigned int>{NewBorn_Size, Mutation_Rate}(eng);
	}

	double Random::Uniform_Mutation_Rate(double mu_rate)
	{
		return uniform_real_distribution<double>{mu_rate/2.0, mu_rate*2.0}(eng);
	}
	double Random::Update_Proliferation_Rate(double Parent_Proliferation_Rate)
	{
		double U = normal_distribution<double>{Parent_Proliferation_Rate, 0.000001 }(eng);
		if(U > 1.0)
		{
			U = 1.0;
		}
		if (U < 0.0)
		{
			U = 0.0;
		}
		
		return U;
	}

	double Random::Uniform_Mutation_Rate_2(double Parent_mu_rate)
	{

		double U;// = normal_distribution<double>{Parent_mu_rate, mut_rate }(eng);
		U = gsl_ran_beta (r_global, 0.1, 100.0);
		//U = gsl_cdf_beta_P (Parent_mu_rate, 2.0,  2.0);

	
			//cout << " parent MR " << Parent_mu_rate << " Child  "<< U << endl;
			//getchar();
		
		
		return U;
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
			cout << "\tSelective_Pressure: " << CE -> feedback << endl;
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
		
			CE -> Tumour -> push_back( get_Clone_DS() );// We update the size that way
			//CE -> Tumour -> back() -> Progenitor_Ready_to_Reproduce = false;
			CE -> Tumour -> back() -> Initiall_Expasion_Period = true;
			CE -> Tumour -> back() -> Clone_Size = 1;
			//CE -> Tumour -> back() -> Time_of_Recovery = r.Recovery_After_Replication() + 48; // this should be a function
			
			CE -> Tumour -> back() -> In_G0_Phase = false;
			CE -> Tumour -> back() -> In_G1_Phase = true;
			CE -> Tumour -> back() -> In_S_Phase = false;
			CE -> Tumour -> back() -> In_G2_Phase = false;
			CE -> Tumour -> back() -> In_M_Phase = false;
			
			CE -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

 
			CE -> Tumour -> back() -> Generation_ID_Counter = 0; 
			CE -> Tumour -> back() -> Generation_ID = "P-0:0";

			CE -> Population_Size++;
			
		}
		else
			throwX( "WARNING: You send to carcinogenesis() a NULL pointer, send a Clonal_Expansion pointer " );

	}// end function

	void Default_Drug_Parameters(unique_ptr<Drug>  const & CTX)
	{


		CTX -> maximum_tolareted_dose = 20.0; // ug/ml
		CTX -> minimum_dose = 10.0; // ug/ml
		CTX -> half_life = 22.0; //hr
		//CTX -> plasma_volume = 3000.0; //ml
		CTX -> dosage = 300.0 * 1000.0 ; //ug
		CTX -> absorption_Fraction = 0.12; //(unitless)
		CTX -> interval = 168;
		
		CTX -> elimintation_constant = -log(0.5)/(CTX -> half_life); // 1/hr

		CTX -> drug_in_system = CTX -> absorption_Fraction * CTX -> dosage;  // ug
		CTX -> concentration = CTX -> drug_in_system/ CTX -> plasma_volume;
		//CTX -> dx = 1.0; //hr
		CTX -> hour_counter = 0;
		CTX -> ingestion = 0.0;

		CTX -> bioavailavility = 90.0;
		
		CTX -> cycles = 100;
		CTX -> repetition = 12; //hours 672
		CTX -> intake_days_in_hours -> push_back(1);
		CTX -> intake_days_in_hours -> push_back(8);//336
		CTX -> time_to_cell_death = 8; //336
		CTX -> time_to_mutation = 4; //168

		get<0>(CTX -> MR_Register).resize(CTX -> time_to_mutation);
		get<0>(CTX -> DR_Register).resize(CTX -> time_to_cell_death);

		
		get<0>(CTX -> MR_Register)[0]= (CTX -> concentration);
		get<1>(CTX -> MR_Register) = 1;
		get<2>(CTX -> MR_Register) = CTX -> time_to_mutation;

		get<0>(CTX -> DR_Register)[0] = (CTX -> concentration);
		get<1>(CTX -> DR_Register) = 1;
		get<2>(CTX -> DR_Register) = CTX -> time_to_cell_death;
		// -> push_back( CTX -> concentration) ;
		//get<2>(CTX -> MR_Register) = 1;
		//get<3>(CTX -> MR_Register) = CTX -> time_to_mutation;


	}
	
	void print_Drug_Parameters(unique_ptr<Drug>  const & CTX)
	{
		cout << "\n\n##################################\n\n DRUG DS \n{ " << endl;
		cout << "\tMax Tolerated Dose: " << CTX -> maximum_tolareted_dose << endl;
		cout << "\tMin Dose To Clinical Benefit: " << CTX -> minimum_dose << endl;
		cout << "\tDrug Half Life: " << CTX -> half_life << endl;
		cout << "\tPlasma Volume: " << CTX -> plasma_volume << endl;
		cout << "\tDosage: " << CTX -> dosage << endl;
		cout << "\tAbsortion Fraction: " << CTX -> absorption_Fraction << endl;
		cout << "\tInterval: " << CTX -> interval << endl;
		cout << "\tElimination Constant: " << CTX -> elimintation_constant << endl;
		cout << "\tDrug in System: " << CTX -> drug_in_system << endl;
		cout << "\tConcentration: " << CTX -> concentration << endl;
		//cout << "\tdx: " << CTX -> dx << endl;
		cout << "\thour counter: " << CTX -> hour_counter << endl;
		cout << "\tingestion: " << CTX -> ingestion << endl;
		cout << "\tbioavailavility: " << CTX -> bioavailavility << endl;
		cout << "\tcycles: " << CTX -> cycles << endl;
		cout << "\trepetition: " << CTX -> repetition << endl;
		cout << "\tIntake Days [in hours]: [" ;
		for (vector<unsigned int>::const_iterator i = CTX -> intake_days_in_hours -> begin(); i != CTX -> intake_days_in_hours -> end(); ++i)
    		cout << *i << ' ';
    	cout << "]" <<endl;
		cout << "\tTime to Cell Death: " << CTX -> time_to_cell_death << endl;
		cout << "\tTime to Mutation: " << CTX -> time_to_mutation << endl;

		cout << "\tMR Register: [ " << get<0>(CTX -> MR_Register).at(0) 
			 << " , " << get<1>(CTX -> MR_Register) 
			 << ", " << get<2>(CTX -> MR_Register)  << "]" << endl;

		cout << "\tDR Register: [ " << get<0>(CTX -> DR_Register).at(0) 
			 << " , " << get<1>(CTX -> DR_Register) 
			 << ", " << get<2>(CTX -> DR_Register)  << "]" << endl;

		cout << "}\n#################################" << endl;

	}


	


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

			cout << "\tNumber_of_Memebers_to_Start_Heterrogeneity: " << CL -> Number_of_Memebers_to_Start_Heterogeneity << endl;
			cout << "\tGeneration ID Counter " << CL -> Generation_ID_Counter << endl;
			cout << "\tExtinct " << CL -> clone_extinct << endl;
			cout << "\tMutation Rate " << CL -> Mutation_Rate << endl;
			cout << "\tNumber of Mutations " << CL -> Number_of_Mutations << endl;
			cout << "\tClone_Size: " << CL -> Clone_Size << endl;
			cout << "\tInitiall_Expasion_Period: " << CL -> Initiall_Expasion_Period << endl;
			cout << "\tIn_G0_Phase: " << CL -> In_G0_Phase << endl;
			cout << "\tIn_G1_Phase: " << CL -> In_G1_Phase << endl;
			cout << "\tIn_S_Phase: " << CL -> In_S_Phase << endl;
			cout << "\tIn_G2_Phase: " << CL -> In_G2_Phase << endl;
			cout << "\tIn_M_Phase: " << CL -> In_M_Phase << endl;
			cout << "\tDEATH RATE: " << CL -> P_Expansion[0] << endl;
			cout << "\tPROL RATE: " << CL -> P_Expansion[1] << endl;
			cout << "\tNT RATE: " << CL -> P_Expansion[2] << endl;
			cout << "\tRemaining_Time_in_G1_Phase: " << CL -> Remaining_Time_in_G1_Phase << endl;
			cout << "\tRemaining_Time_in_S_Phase: " << CL -> Remaining_Time_in_S_Phase << endl;
			cout << "\tRemaining_Time_in_G2_Phase: " << CL -> Remaining_Time_in_G2_Phase << endl;
			cout << "\tRemaining_Time_in_M_Phase: " << CL -> Remaining_Time_in_M_Phase << endl;
			cout << "\tName: " << CL -> Generation_ID << endl;
			
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
		unsigned int mutant_cells = 0;


		/*
			1) Let's get basline parameters
		*/
		double P_DR =  CE -> Tumour -> at(Generation_ID) -> P_Expansion[0];
		double P_NB =  CE -> Tumour -> at(Generation_ID) -> P_Expansion[1];
		
		/*
			2) Apply enviromental penalty to growth variables
		*/

		// if the penalty greater to PR lets cap it to 0 
		if(CE -> feedback >= P_NB)
		{
			P_NB = 0;
		}
		else
		{
			P_NB -= CE -> feedback;
		}

		/*
			3) Adjust Probability Mass
		*/
		double P_NT =  1.0 - (P_DR + P_NB );
		//structure that holds potential values
		double p[] = {P_DR, P_NB, P_NT};

		/*
			4) Sample from multinomial as
			[Died, Newborn, Other] ~ multinom(p{p_dr,p_nb,p_nt}, Clone_Size)
		*/
		// n will hold the numbers assoiated in the sample 
		gsl_ran_multinomial ( r_global, K, CE -> Tumour -> at (Generation_ID) -> Clone_Size, p, n);
		//cout << "DY: " << n[0] << " NB: " << n[1] << " NT: "<< n[2] << endl; 
	

		get<0>(Dying_and_Newborn) = n[0];// For Dying Clones 		
		get<1>(Dying_and_Newborn) = n[1]; // For Newborn Clones 
		
		/*
			5) Check for Mutants
		*/
		if(get<1>(Dying_and_Newborn) > 0)
		{	
			get<2>(Dying_and_Newborn) = r.Binomial_Mutants(
															get<1>(Dying_and_Newborn),
															CE -> Tumour -> at(Generation_ID) -> Mutation_Rate
														   );
			// Announce that there are mutants
			mutant_cells = get<2>(Dying_and_Newborn); // Mutant cells
			get<3>(Dying_and_Newborn) =  true;// flag that indicates that there are mutant cells

			// Adjust Newborn gain as
			// Newborn - Mutant cells
			if(mutant_cells > get<1>(Dying_and_Newborn))
			{
				cout << "More mutant cells than newborn? " << " MC: " << mutant_cells << " NB: " << get<1>(Dying_and_Newborn) << endl;
			}
			else
			{
				get<1>(Dying_and_Newborn) -= mutant_cells;
			}

		}
		

		// Check for clonal extinction
		int o = (int) ( CE -> Tumour -> at(Generation_ID) -> Clone_Size + n[1]) - (int) (n[0] + mutant_cells);
			
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

	/*
		Oeverla Penalty to the population
	*/
	double Generate_PR_Drug_Penalty(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		double chemo_PR_penalty;
		for(ith_drug=0; ith_drug < size; ith_drug++)
			if( Therapy -> AntiCancer_Agents -> at(ith_drug) ->  Penalty_PR )
				chemo_PR_penalty += Therapy -> AntiCancer_Agents -> at(ith_drug) -> concentration/100.0;

		Therapy -> Total_Drug_Concentration = chemo_PR_penalty;
		//cout << "PR penalty: " << Therapy -> Total_Drug_Concentration << endl;
	
		return Therapy -> Total_Drug_Concentration;
	}

	double Generate_MR_Drug_Penalty(unique_ptr<Treatment> const & Therapy, double baseline_MR, string Generation_ID )
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		double MR_Penalty = 0.0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
				

			//if(Generation_ID.compare ("P-0:0") == 0)
			//	cout << "INITIAL VALUES " << "ID" <<  Generation_ID 
			//		 << "  drug " << ith_drug
			//		 << " baseline_MR " << baseline_MR 
			//		 << "  Value <3> " <<  get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register) 
			//		 << "  FLAG  " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> Penalty_MR 
			//		 << endl;

			if( get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register) && (Therapy -> AntiCancer_Agents -> at(ith_drug) -> Penalty_MR ) )
			{
				
				MR_Penalty += get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)[get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)]/100.0;
				//if(Generation_ID.compare ("P-0:0") == 0){
				//	cout << " ID" <<  Generation_ID 
				//		 << " drug " << ith_drug 
				//		 << " Val " << get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)[get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)] 
				//		 << "idx " << get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)
				//		 << " MR PENALTY " << MR_Penalty << endl; getchar();}
			}

			
		}


		//if(Generation_ID.compare ("P-0:0") == 0)
		//	cout<< " Total penalty " << baseline_MR + MR_Penalty << endl;

		return (baseline_MR + MR_Penalty);
	}

	double Generate_DR_Drug_Penalty(unique_ptr<Treatment> const & Therapy, double baseline_DR, string Generation_ID )
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		double DR_Penalty = 0.0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
				

			//if(Generation_ID.compare ("P-0:0") == 0)
			//	cout << "INITIAL VALUES " << "ID" <<  Generation_ID 
			//		 << "  drug " << ith_drug
			//		 << " baseline_MR " << baseline_MR 
			//		 << "  Value <3> " <<  get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register) 
			//		 << "  FLAG  " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> Penalty_MR 
			//		 << endl;

			if( get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register) && (Therapy -> AntiCancer_Agents -> at(ith_drug) -> Penalty_DR ) )
			{
				
				DR_Penalty += get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register)[get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register)]/100.0;
				//if(Generation_ID.compare ("P-0:0") == 0){
				//	cout << " ID" <<  Generation_ID 
				//		 << " drug " << ith_drug 
				//		 << " Val " << get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)[get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)] 
				//		 << "idx " << get<1>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)
				//		 << " MR PENALTY " << MR_Penalty << endl; getchar();}
			}

			
		}


		//if(Generation_ID.compare ("P-0:0") == 0)
		//	cout<< " Total penalty " << baseline_DR + DR_Penalty << endl;

		return (baseline_DR + DR_Penalty);
	}

	//DR penalty


	double apply_chemo_penalty(double chemo_PR_penalty, double PR)
	{
		if(chemo_PR_penalty >= PR)
			PR = 0.0;
		else
			PR -= chemo_PR_penalty;
		return PR;
	}

	double apply_size_penalty(double size_penalty, double PR)
	{
		if(size_penalty >= PR)
			PR = 0.0;
		else
			PR -= size_penalty;
		return PR;
	}


	tuple<unsigned int,unsigned int, unsigned int, bool, bool> Basic_Clonal_Expansion_Treatment(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int hours, unique_ptr<Treatment> const & Therapy )
	{
		tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn (0, 0, 0, false, true);

		get<3>(Dying_and_Newborn) =  false;
		get<4>(Dying_and_Newborn) =  true;
		unsigned int n[3];
		unsigned int mutant_cells = 0;

		//if((CE -> Tumour -> at(Generation_ID) -> Generation_ID).compare ("P-0:0") == 0)
		//	cout << "Hours " << hours << " MR " <<  CE -> Tumour -> at(Generation_ID) -> Mutation_Rate ;

		double MR = CE -> Tumour -> at(Generation_ID) -> Mutation_Rate;
		//cout << "MR -> " << MR <<endl;
		/*
			1) Let's get basline parameters
		*/
		double P_DR =  CE -> Tumour -> at(Generation_ID) -> P_Expansion[0];
		double P_NB = CE -> Tumour -> at(Generation_ID) -> P_Expansion[1];
		//cout << "ID: " << Generation_ID << " PB(t-1): " << P_NB ;


		//if((CE -> Tumour -> at(Generation_ID) -> Generation_ID).compare ("P-0:0") == 0)
		//	cout << "Hours " << hours << " DR " <<  P_DR ;


		/*
			2) Penlaty on PR
		*/
		//double chemo_PR_penalty =  Generate_PR_Drug_Penalty( Therapy );

		//cout << " Conc: " << chemo_PR_penalty ;
		if( !(Therapy -> Trial_Done) )
			P_NB = apply_chemo_penalty( Generate_PR_Drug_Penalty( Therapy ),  P_NB );

		P_NB = apply_size_penalty(CE -> feedback, P_NB );

		if( !(Therapy -> Trial_Done) )
			MR =  Generate_MR_Drug_Penalty( Therapy, MR, CE -> Tumour -> at(Generation_ID) -> Generation_ID );

		if( !(Therapy -> Trial_Done) )
			P_DR =  Generate_DR_Drug_Penalty( Therapy, P_DR, CE -> Tumour -> at(Generation_ID) -> Generation_ID );

	//	if((CE -> Tumour -> at(Generation_ID) -> Generation_ID).compare ("P-0:0") == 0)
	//		cout  << " NEW DR " <<  P_DR << endl;
		//cout << " with  penlaty: >> " << MR << " Ex: " <<CE -> Tumour -> at(Generation_ID) -> clone_extinct << endl;

		double P_NT =  1.0 - (P_DR + P_NB );
		double p[] = {P_DR, P_NB, P_NT};


		
		gsl_ran_multinomial ( r_global, K, CE -> Tumour -> at (Generation_ID) -> Clone_Size, p, n);
		//cout << "DY: " << n[0] << " NB: " << n[1] << " NT: "<< n[2] << endl; 
	


			get<0>(Dying_and_Newborn) = n[0];
			// For Newborn Clones 
			get<1>(Dying_and_Newborn) = n[1];

			if(get<1>(Dying_and_Newborn) > 0)
			{	
				get<2>(Dying_and_Newborn) = r.Binomial_Mutants( get<1>(Dying_and_Newborn), MR );
				mutant_cells = get<2>(Dying_and_Newborn);
				get<3>(Dying_and_Newborn) =  true;
				// Adjust Newborn gain as
				// Newborn - Mutant cells
				if(mutant_cells > get<1>(Dying_and_Newborn))
				{
					cout << "More mutant cells than newborn? " << " MC: " << mutant_cells << " NB: " << get<1>(Dying_and_Newborn) << endl;
				}
				else
				{
					get<1>(Dying_and_Newborn) -= mutant_cells;
				}
			}

			
			int o = (int) ( CE -> Tumour -> at(Generation_ID) -> Clone_Size + n[1]) - (int) (n[0] + mutant_cells);
			
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


	void carcinogenesis_from_Driver(unique_ptr<Clonal_Expansion> const  & CE,  int Generation_ID , unsigned int years, unsigned int hours)
	{

		//cout << "HERE A " << CE -> Population_Size <<endl;
		if(CE)// is the pointer NUL?
		{	
		

			double mr = CE -> Tumour -> at(Generation_ID) -> Mutation_Rate;
			unsigned long long int NOM = CE -> Tumour -> at (Generation_ID) -> Number_of_Mutations;
			double pr = CE -> Tumour -> back() -> P_Expansion[1];


			//CE -> Tumour -> at(Generation_ID) -> Number_of_Driver_Mutations_in_Cell_Division++;
			string Clone_Name = Generate_Clone_Name_ID(CE -> Tumour -> at(Generation_ID) -> Generation_ID_Counter, CE -> Tumour -> at(Generation_ID)  -> Generation_ID, years, hours);
			CE -> Tumour -> at(Generation_ID) -> Generation_ID_Counter++;


			CE -> Tumour -> push_back( get_Clone_DS() );// We update the size that way
			CE -> Tumour -> back() -> Generation_ID = Clone_Name;
			
			CE -> Tumour -> back() -> Initiall_Expasion_Period = false;//true
			CE -> Tumour -> back() -> Clone_Size = 1;
			
			//CE -> Tumour -> back() -> Time_of_Recovery = r.Recovery_After_Replication() + 48; // this should be a function
			CE -> Tumour -> back() -> Number_of_Mutations = NOM + 1;
			//CE -> Tumour -> back() -> Time_to_Division = r.G1() + r.S() + 1 + r.G2();

			CE -> Tumour -> back() -> Remaining_Time_in_G1_Phase = r.G1();
			CE -> Tumour -> back() -> Remaining_Time_in_S_Phase = r.S();
			CE -> Tumour -> back() -> Remaining_Time_in_G2_Phase = r.G2();
			CE -> Tumour -> back() -> Remaining_Time_in_M_Phase = 1;

			CE -> Tumour -> back() -> In_G0_Phase = false;
			CE -> Tumour -> back() -> In_G1_Phase = true;
			CE -> Tumour -> back() -> In_S_Phase = false;
			CE -> Tumour -> back() -> In_G2_Phase = false;
			CE -> Tumour -> back() -> In_M_Phase = false;

			//CE -> Tumour -> back() -> Remaining_Time_in_Recovery_Phase = 0;
			

			CE -> Tumour -> back() -> clone_extinct = false;
			
			CE -> Tumour -> back() -> Generation_ID_Counter = 0; 


			CE -> Tumour -> back() -> Mutation_Rate = r.Uniform_Mutation_Rate_2(mr);

			
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

	void map_Feedback_orig( unique_ptr<Clonal_Expansion> const & CE )
	{
		//cout << "DR : "<< death_rate << " Pop: " << Pop_Size << " CV: " << Current_Population_Size << endl;

		CE -> feedback =  0.0 + (DIFF - 0.0) * (( (double) CE -> Population_Size - 0.0) / ((double) PS - 0.0));	
	}

	/*
		Penlaty Model

	*/
	void map_Feedback( unique_ptr<Clonal_Expansion> const & CE, double scale  )
	{
		unsigned int ith_clone = 0;
		unsigned int size = CE -> Tumour  -> size();
		double avg = 0.0;
		unsigned long long int k = 0;
		for (ith_clone = 0; ith_clone < size ; ith_clone++)
			if( !(CE -> Tumour -> at (ith_clone) -> clone_extinct) )
			{
				avg += CE -> Tumour -> at(ith_clone) -> P_Expansion[1] ;
				k++;
			}

		double diff = avg/((double) k * scale);

		CE -> feedback = 0.0 + (diff - 0.0) * (( (double) CE -> Population_Size - 0.0) / ((double) PS - 0.0));	
	}//map_feedback

	void Size_Dependent_Penalty(unique_ptr<Clonal_Expansion> const & CE)
	{
		
		switch (PENALTY)
		{
			case 1:
			{
				map_Feedback(CE, 1.0);
				break;
			}

			case 2:
			{
				map_Feedback(CE, 2.0);
				break;
			}

			case 3:
			{
				map_Feedback_orig( CE );
				break;
			}

		}	
			
	}//funtion
				

	void Mutational_Effect_From_Mutants(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int Number_of_Mutants, unsigned int years, unsigned int hours  )
	{
		if(Number_of_Mutants > 0) // redundant
		{
			//double muatational_effect = 0.0;
			
			for(unsigned int i = 0; i < Number_of_Mutants; i++)
				carcinogenesis_from_Driver(CE,  Generation_ID,  years,  hours);
				
		}//if
	}

	void compute_Mutations(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID ,tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn, unsigned int years, unsigned int hours )
	{
		/*
			Version 1 => requires that all new mutants are considered new clones
		*/
		if( get<3>(Dying_and_Newborn) )
		{
			Mutational_Effect_From_Mutants(CE, Generation_ID, get<2>(Dying_and_Newborn), years, hours );

		}
		/*
			From Basic Clonal expansion we updated Newbrn - mutants = Gain
			Therfore
			CS = (CS - Death) Gain 

		*/
		unsigned int New_poulation_Size = ( CE -> Tumour -> at(Generation_ID) -> Clone_Size  - get<0>(Dying_and_Newborn) ) 
										+ 	get<1>(Dying_and_Newborn) ;

		// Update population size by
		// Pop size =Pop_Size + (CS[t] - CS[t-1]) 
		CE -> Population_Size +=  ( (unsigned long long int) New_poulation_Size - (unsigned long long int) CE -> Tumour  -> at(Generation_ID) -> Clone_Size) ;
		CE -> Tumour -> at(Generation_ID) -> Clone_Size = (unsigned int) New_poulation_Size;
		
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
			
			// Only if there are new mutants
			if(get<4>(Dying_and_Newborn))
			{
				compute_Mutations(CE, Generation_ID, Dying_and_Newborn, years, hours );
			}
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = false;
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase = true;
			
		}

	}




	void Mitosis_Treatment(unique_ptr<Clonal_Expansion> const & CE, int Generation_ID, unsigned int hours, unsigned int years, unique_ptr<Treatment> const & Therapy )
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

			
			tuple<unsigned int, unsigned int, unsigned int, bool, bool> Dying_and_Newborn = Basic_Clonal_Expansion_Treatment( CE, Generation_ID, hours, Therapy ); 
			if(get<4>(Dying_and_Newborn))
			{
				compute_Mutations(CE, Generation_ID, Dying_and_Newborn, years, hours );
			}
			CE -> Tumour -> at (Generation_ID) -> In_M_Phase = false;
			CE -> Tumour -> at (Generation_ID) -> In_G1_Phase = true;
			
		}

	}


	/*
		This funtion stores all the evolution of the clones
	*/
	void store_clonal_evolution(unique_ptr<Clonal_Expansion> const & CE, unsigned int elapsed_hours,string BasePath , unsigned int hours = 0 )
	{
		unsigned int ith_clone = 0;
		unsigned int size = CE -> Tumour  -> size() ;


		if(hours % 100 == 0)
		{	//Create a file for each clone 
			for (ith_clone = 0; ith_clone < size ; ith_clone++)
			{
				ofstream myfile;
				string file_path = BasePath +"/Clonal_Evolution/";
				file_path.append(CE -> Tumour -> at (ith_clone) -> Generation_ID);
				file_path.append(".txt");
				myfile.open( file_path, ios_base::app );
    			myfile << elapsed_hours + 1 << "\t" << CE -> Tumour -> at (ith_clone) -> Clone_Size << "\n";	
			}//for
		}//if hours wait
	}//function

	void update_population(unique_ptr<Clonal_Expansion> const & CE, unsigned int hours, unsigned int years)
	{
		unsigned int size = CE -> Tumour  -> size() ;
		unsigned int ith_clone = 0;
				
		for (ith_clone = 0; ith_clone < size ; ith_clone++) 			//(1) Is my clone not extinct?	
			if( !(CE -> Tumour -> at (ith_clone) -> clone_extinct) )	// (1.a) if not, is my clone in free division?
				if(CE -> Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
					Initial_Expansion_Mitosis(CE, ith_clone, hours, years);
				else
					Mitosis( CE, ith_clone, hours, years );
			else
				continue;												// Nothin to do with this clone, is dead
			
	}//function

	/*
		Print the state of the simulation
	*/
	void print_Status(unique_ptr<Clonal_Expansion> const & CE, unsigned int hours, unsigned int years, int myID,  bool each_100 = false )
	{
		
		if(each_100 && (hours % 100 == 0))
		{	
			usleep( (myID*1) + 10);
			cout << " \n\n ACTIVE CELLS " <<  CE -> Population_Size  
				 << " CLONES " << CE -> Tumour -> size() 
				 << "   H: " << hours  
				 << " Y: " << years 
			 	<< " FD: " << CE -> feedback 
			 	<< " ID: " << myID  
			 	<< endl;
			 }
		 else if(each_100 == false)
			cout << " \n\n ACTIVE CELLS " <<  CE -> Population_Size  
				 << " CLONES " << CE -> Tumour -> size() 
				 << "   H: " << hours  
				 << " Y: " << years 
			 	<< " FD: " << CE -> feedback
			 	<< " ID: " << myID  
			 	<< endl;
		
			
	}// print command

	unsigned int Save_Growth_Curve(unique_ptr<Clonal_Expansion> const & CE, unsigned int elapsed_hours, ofstream &Tumour_Evolution )
	{
		Tumour_Evolution << CE -> Population_Size << "\t" << elapsed_hours << "\n";
		return elapsed_hours + 1;
	}

	

	void Store_Population_Stats(unique_ptr<Clonal_Expansion> const & CE, string path)
	{
		unsigned int ith_clone = 0;
		ofstream Pop_Stats;
		Pop_Stats.open (path);
		Pop_Stats << "id\t" << "Clone_Size\t" << "Proliferation_Rate\t" << "Mutation_Rate\t" << "Extinct\t" << "G_ID\t"<<"Number_of_Mutations" <<"\n";
  		
  		for( ith_clone = 0; ith_clone < CE-> Tumour -> size() ; ith_clone ++)
  		{
  			Pop_Stats << ith_clone << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> Clone_Size << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> P_Expansion[1] << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Mutation_Rate << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> clone_extinct << "\t"
  					  << CE -> Tumour -> at (ith_clone) -> Generation_ID << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Number_of_Mutations  
  					  << "\n" ;


  		}

  		Pop_Stats.close();
	}

	unsigned int Abort_Condition(unique_ptr<Clonal_Expansion> const & CE, unsigned int times_to_wait)
	{
		if( CE -> Population_Size >  DETECTABLE_POPULATION_SIZE )
  			times_to_wait++;

  		return times_to_wait;
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
	void compute_Tumour_Evolution(unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID)
	{
		using namespace std;

		/* 
			Time Variables 
		*/
		unsigned int seconds = 0;
		unsigned int hours = 0;
		unsigned int years = 0;
		unsigned int elapsed_hours = 0;
		unsigned int times_to_wait = 0;

		
		/*
			File Stream	
		*/ 
  		ofstream Tumour_Evolution;
  		string Growth = BasePath + "/Growth/" + "Initial_Growth.txt"; 
  		string Stats = BasePath + "/Stats/" + "All_Clones_Prior.txt";
  		// cconsider havinf as iput the file names for this and other files
  		Tumour_Evolution.open (Growth);

		carcinogenesis(CE);
		print_Clone_DS( CE -> Tumour -> at(0) );


		while( ((times_to_wait < STOP_GROWTH_AFTER_DIAGNOSIS) && (years < 50) ) && (CE -> Population_Size > 0) )
		{
			seconds += dt;
			/*
				Each hour compute...
			*/
			if(seconds == 3600)
			{
				seconds = 0; hours ++;
				
				update_population(CE, hours, years);

				if( SAVE_CLONAL_EVOLUTION )
					store_clonal_evolution( CE, elapsed_hours, BasePath ); // to save the clone
				// Consider higher valuesof clones to save only
			
				Size_Dependent_Penalty( CE );
			
				if(PRINT)
					print_Status(CE, hours, years, myID, true );

				elapsed_hours = Save_Growth_Curve(CE, elapsed_hours, Tumour_Evolution ); // recursive function

				times_to_wait = Abort_Condition(CE, times_to_wait);
				
			}//seconds
			if(hours == 8764)//8764
			{
				hours = 0; years ++;
			}
		}//while

		cout << "FINISH Growth Phase ID: " << myID << endl;
		Tumour_Evolution.close();

		Store_Population_Stats( CE, Stats );
		
	}// end of function

/*
	bool check_intake_day(unique_ptr<Drug> const & CTX)
	{

		bool intake = false;
		for (vector<unsigned int>::const_iterator i = CTX -> intake_days -> begin(); i != CTX -> intake_days -> end(); ++i)
		{
			if(CTX -> hour_counter == *i)
			{
				intake = true;
				break;
			}
		}
		return intake;
	}// End model
	*/

	
	void update_population_treatment(unique_ptr<Clonal_Expansion> const & CE,  unique_ptr<Treatment> const & Therapy, unsigned int hours, unsigned int years)
	{
		unsigned int size = CE -> Tumour  -> size() ;
		unsigned int ith_clone = 0;
				
		for (ith_clone = 0; ith_clone < size ; ith_clone++) 			//(1) Is my clone not extinct?	
			if( !(CE -> Tumour -> at (ith_clone) -> clone_extinct) )	// (1.a) if not, is my clone in free division
				if(CE -> Tumour -> at (ith_clone) -> Initiall_Expasion_Period )
					Initial_Expansion_Mitosis(CE, ith_clone, hours, years);
				else
					Mitosis_Treatment( CE, ith_clone, hours, years,  Therapy);
			else
				continue;

	}

	

	void update_drug_registers(unique_ptr<Drug> const & CTX)
	{
		//Update Registers
		if(get<1>(CTX -> MR_Register) == get<2>(CTX -> MR_Register))
		{
			get<1>(CTX -> MR_Register) = 0;
			get<3>(CTX -> MR_Register) = true;
			/*
			cout << "MR: " << get<0>(CTX -> MR_Register).at(0) << " " << get<0>(CTX -> MR_Register).at(1) << endl;
			vector<double>::size_type sz = get<0>(CTX -> MR_Register).size();
			cout << "\n REGISTER MR FULL size : "  << sz  << " [ ";
			for (vector<double>::const_iterator i = get<0>(CTX -> MR_Register).begin(); i != get<0>(CTX -> MR_Register).end(); ++i)
    		cout << *i << ' ';
    		cout << " \n"<< endl;
    		getchar();
    		*/

		}
		get<0>(CTX -> MR_Register)[get<1>(CTX -> MR_Register)] = CTX -> concentration;
		get<1>(CTX -> MR_Register)++;
		//cout <<"MR C:  "<< get<1>(CTX -> MR_Register) << endl;

		if(get<1>(CTX -> DR_Register) == get<2>(CTX -> DR_Register))
		{
			get<1>(CTX -> DR_Register) = 0;
			get<3>(CTX -> DR_Register) = true;
			//vector<double>::size_type sz = get<0>(CTX -> DR_Register).size();
			//cout << "\n REGISTER DR FULL size: " << sz << " [ " ;
			//for (vector<double>::const_iterator i = get<0>(CTX -> DR_Register).begin(); i != get<0>(CTX -> DR_Register).end(); ++i)
    		//	cout << *i << ' ';
    		//cout << " \n"<< endl;
    		//getchar();
		}
		get<0>(CTX -> DR_Register)[get<1>(CTX -> DR_Register)] = CTX -> concentration;
		get<1>(CTX -> DR_Register)++;
		//cout <<"DR C:  "<< get<1>(CTX -> MR_Register) << endl;
	}

	bool is_extinct(unique_ptr<Clone> const  & Clone)
	{
		return (Clone -> clone_extinct);
	}


	void Purge_Death_Clones(unique_ptr<Clonal_Expansion> const & CE)
	{
		
		cout << "REMOVING DEATH CLONES" << endl;

		CE ->Tumour -> erase ( remove_if(CE ->Tumour -> begin(), CE ->Tumour -> end (), is_extinct ), CE ->Tumour -> end ()  );



		//vector<unique_ptr<Clone> >::iterator it;
		//for(it = CE ->Tumour -> begin(); it != CE ->Tumour ->  end();)
    	//	if ((*it)->clone_extinct )
        //		it = CE -> Tumour -> erase(it);  // Returns the new iterator to continue from.
    	//	else
        //		++it;

        cout << "Effective Number of Clones: " << CE -> Tumour -> size() << endl;
	}//



	vector<string> &split(const string &s, char delim, vector<string> &elems) 
	{
    	stringstream ss(s);
    	string item;
    	while (getline(ss, item, delim)) 
        	elems.push_back(item);
   	
    	return elems;
	}


	vector<string> split(const string &s, char delim) 
	{
    	vector<string> elems;
    	split(s, delim, elems);
    	return elems;
	}

	unsigned int intake_days_to_hours(unsigned int day)
	{
		unsigned int hour = 0;
		if(day != 1)
		{
			hour = day * 24;
		}
		return hour;

	}

	void StringSeq_to_int_Vec(string sequence, unique_ptr<Treatment> const & Therapy)
	{

		vector<string> Sseq = split(sequence, ',');
		vector<unsigned int> Usequence;

		for (auto &s : Sseq) 
		{
    		stringstream parser(s);
    		unsigned int x = 0;

    		parser >> x;

    		Usequence.push_back(x);
		}
    	
    	for(vector<unsigned int>::size_type i = 0; i != Usequence.size(); i++) 
    		Therapy -> AntiCancer_Agents -> back() -> intake_days_in_hours -> push_back( intake_days_to_hours( Usequence[i] ) );
		

		/*
		for(vector<unsigned int>::iterator it = 
			Therapy -> AntiCancer_Agents -> back() -> intake_days -> begin(); 
			it != Therapy -> AntiCancer_Agents -> back() -> intake_days -> end(); 
			++it) 
		{
    		cout << *it << " ";
		}cout << "\n";
		*/

	}//

	void Drug_Penalty_Type( unique_ptr<Drug> const & AntiCancer_Agent, string type)
	{
		if(type.compare("PR") == 0)
			AntiCancer_Agent -> Penalty_PR = true;
	
		if(type.compare("MR") == 0)
			AntiCancer_Agent -> Penalty_MR = true;

		if(type.compare("DR") == 0)
			AntiCancer_Agent -> Penalty_DR= true;

		if(type.compare("PR_MR") ==  0)
		{
			//cout << " HERE YRS" << endl;
			//getchar();
			AntiCancer_Agent -> Penalty_PR = true;
			AntiCancer_Agent -> Penalty_MR = true;
		}

		if(type.compare("DR_MR") == 0)
		{
			AntiCancer_Agent -> Penalty_DR = true;
			AntiCancer_Agent -> Penalty_MR = true;
		}

		if(type.compare("PR_DR") == 0 )
		{
			AntiCancer_Agent -> Penalty_PR = true;
			AntiCancer_Agent -> Penalty_DR = true;
		}

		if(type.compare("ALL") == 0)
		{
			AntiCancer_Agent -> Penalty_PR = true;
			AntiCancer_Agent -> Penalty_DR = true;
			AntiCancer_Agent -> Penalty_MR = true;
		}


	}

	void Store_Drugs(unique_ptr<Treatment> const & Therapy,  unordered_map<string,string> CTX_Scheme )
	{
		if(Therapy)
		{
			//size_type sz;
			string Drug = "Drug_";
			for(unsigned int i = 0; i < Therapy -> Number_of_Drugs; i++)
			{ 
				//Drug += to_string(i+1) + "_";
				//cout << Drug << endl;
				Therapy -> AntiCancer_Agents -> push_back( get_Drug() );
				Therapy -> AntiCancer_Agents -> back() -> Family = CTX_Scheme[Drug + to_string(i+1) + "_" + "Family"];
				Therapy -> AntiCancer_Agents -> back() -> DrugName = CTX_Scheme[Drug + to_string(i+1) + "_" + "Name"];
				Therapy -> AntiCancer_Agents -> back() -> maximum_tolareted_dose = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "MTD"] ); // ug/ml
				//cout << Therapy -> AntiCancer_Agents -> back() -> maximum_tolareted_dose << endl;

				Therapy -> AntiCancer_Agents -> back() -> minimum_dose = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "MD"] ); // ug/ml
				//cout << Therapy -> AntiCancer_Agents -> back() -> minimum_dose << endl;

				Therapy -> AntiCancer_Agents -> back() -> dosage = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "Dosage"] );
				//cout << Therapy -> AntiCancer_Agents -> back() -> dosage << endl;

				Therapy -> AntiCancer_Agents -> back() -> absorption_Fraction = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "Absortion"] );
				//cout << Therapy -> AntiCancer_Agents -> back() -> absorption_Fraction << endl;

				StringSeq_to_int_Vec(CTX_Scheme[Drug + to_string(i+1) + "_" + "Days"], Therapy);

				Therapy -> AntiCancer_Agents -> back() -> bioavailavility = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "Bioavailability"] );
				//cout << Therapy -> AntiCancer_Agents -> back() -> bioavailavility << endl;

				Therapy -> AntiCancer_Agents -> back() -> time_to_mutation = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "TTM"] );
				//get<0>( Therapy -> AntiCancer_Agents -> back() ->  MR_Register).resize(Therapy -> AntiCancer_Agents -> back() -> time_to_mutation);

				//get<0>(CTX -> MR_Register).resize(CTX -> time_to_mutation);
				//cout << Therapy -> AntiCancer_Agents -> back() -> time_to_mutation << endl;

				Therapy -> AntiCancer_Agents -> back() -> time_to_cell_death = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "TTD"] );
				Therapy -> AntiCancer_Agents -> back() -> half_life = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "HalfLife"] ); 
				//cout << Therapy -> AntiCancer_Agents -> back() -> time_to_cell_death << endl;
				//Therapy -> AntiCancer_Agents -> back() -> half_life = stod( CTX_Scheme[Drug + to_string(i+1) + "_" + "Penlaty"] );

				Drug_Penalty_Type( Therapy -> AntiCancer_Agents -> back() , CTX_Scheme[Drug + to_string(i+1) + "_" + "Penlaty"]); 

			}// for

		}//Therapy

	}//Function

	//Print Function

	void print_Treatment_DS(unique_ptr<Treatment> const & Therapy)
	{
		cout << "\n\n##################################\n\n TREATMENT DS \n{ " << endl;
		cout << "Scheme: " << Therapy -> Scheme 
  			 << "\n\t # Drugs: " << Therapy -> Number_of_Drugs//Therapy -> Number_of_Drugs
  			 << "\n\t Type_of_Treatment: " << Therapy -> Type_of_Treatment
  			 << "\n\t Repetition [hours]: " << Therapy -> Repetition_Day_in_Hours//Therapy -> Repetition_Day
  			 << "\n\t Cycles: "	   << Therapy -> Cycles//Therapy -> Cycles
  			 << "\n\t dx: " << Therapy -> dx 
  			 << "\n\t Hours Counter: " << Therapy -> hours_counter 
  			 << endl;
  		cout << "}\n#################################" << endl;

  		cout << "\n\n##################################\n\n DRUGS \n{ " << endl;
  		for(unsigned int i = 0; i < Therapy -> Number_of_Drugs; i++)
		{ 
			cout << "DRUG " + to_string(i);
			cout << "\n\t Family: " << Therapy -> AntiCancer_Agents -> at(i) -> Family
				 << "\n\t Drug Name: " << Therapy -> AntiCancer_Agents -> at(i) -> DrugName
				 << "\n\t Max Tolared Dose: " << Therapy -> AntiCancer_Agents -> at(i) -> maximum_tolareted_dose
				 << "\n\t Minimum Dose: " << Therapy -> AntiCancer_Agents -> at(i) -> minimum_dose
				 << "\n\t Dosage: " << Therapy -> AntiCancer_Agents -> at(i) -> dosage
				 << "\n\t Absorption Fraction: " << Therapy -> AntiCancer_Agents -> at(i) -> absorption_Fraction
				 <<"\n\t Intake Days: ";
			for(vector<unsigned int>::iterator it = 
				Therapy -> AntiCancer_Agents -> at(i) -> intake_days_in_hours -> begin(); 
				it != Therapy -> AntiCancer_Agents -> at(i) -> intake_days_in_hours -> end(); 
				++it) 
    				cout << *it << " ";
						
			cout << "\n\t Bioavailavility: " << Therapy -> AntiCancer_Agents -> at(i) -> bioavailavility
				 << "\n\t Time To Modify Mutation Rate: " << Therapy -> AntiCancer_Agents -> at(i) -> time_to_mutation
				 << "\n\t Time To Modify Death Rate: " << Therapy -> AntiCancer_Agents -> at(i) -> time_to_cell_death
				 << "\n\t Ingestion: " << Therapy -> AntiCancer_Agents -> at(i) -> ingestion
				 << "\n\t Half Life: " << Therapy -> AntiCancer_Agents -> at(i) -> half_life
				 << "\n\t Penalty PR: " << Therapy -> AntiCancer_Agents -> at(i) -> Penalty_PR
				  << "\n\t Penalty MR: " << Therapy -> AntiCancer_Agents -> at(i) -> Penalty_MR
				   << "\n\t Penalty DR: " << Therapy -> AntiCancer_Agents -> at(i) -> Penalty_DR
				<< endl;
			}
			 cout << "}\n#################################" << endl;
  			
	}

	/*
		This Function generates the Cehmeotherapy Reading
	*/
	void Assign_CTX_Input_Values(unordered_map<string,string> CTX_Scheme, unique_ptr<Treatment> const & Therapy)
	{
		//unique_ptr<Treatment> const Therapy = get_Treatment_DS();
		/* Requires that the Input File is correct */
		// 1) Read CTX Scheme
		Therapy -> Scheme = CTX_Scheme["Scheme"];
		// 2) Create the space of the number of drugs
		Therapy -> Number_of_Drugs  = stoul ( CTX_Scheme["Number_of_Drugs"], nullptr, 0 );
		// 3) Type of treatment
		Therapy -> Type_of_Treatment = CTX_Scheme["Type"];
		// 4) Repetition Day
		Therapy -> Repetition_Day_in_Hours = stoul (CTX_Scheme["Repetition_Day"], nullptr, 0) * 24; // to convert days to hours
		// 5) Number of Cycles
		Therapy -> Cycles = stoul ( CTX_Scheme["Number_of_Cycles"], nullptr, 0);
		// 6) set dx
		Therapy -> dx = 1.0;
		// 7) Hours counter
		Therapy -> hours_counter = 0;

		//Defualt Parameters
		Therapy -> Maximum_Tolareted_Toxicity = 50.0;
		Therapy -> Plasma_Volume =  3000.0; //ml	

  		Store_Drugs(Therapy, CTX_Scheme);
  		print_Treatment_DS(Therapy);

	}//Function+



	void Read_CTX(string path, unique_ptr<Treatment> const & Therapy)
	{
		cout << "READING FILE: " << path << endl;
		ifstream infile(path);
		string line;
		unordered_map<string,string> CTX_Scheme;
		while (getline(infile, line)) // Reading line by line
		{
    		istringstream iss(line); // define the stream to split tje content
    		vector<string> tokens{
    								istream_iterator<string>{iss}, 
                	      			istream_iterator<string>{} 
                    	  		}; // Store the string in two different positions of the vector
            CTX_Scheme.emplace (tokens[0], tokens[1]);
           // cout << tokens[0] << " " << tokens[1] << endl;
		}//end while
		Assign_CTX_Input_Values( CTX_Scheme, Therapy );

	}// end function read CTX

	

	

	

	void print_Drug_Ingestion(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;

		cout << " INGESTION OF DRUGS " << endl;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			cout << "DRUG: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> DrugName
				 << " >> Ingestion: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion
				 << endl;
		}
		cout << "TOTAL DRUG INGESTED: " << Therapy -> Total_Drug_Ingestion << endl;


	}//function

/*
	void First_Intake(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
	
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			
			Drug_Diffusion( Therapy -> AntiCancer_Agents -> at(ith_drug) );
			
			Therapy -> Total_Drug_Ingestion += Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion;
		}

		print_Drug_Ingestion( Therapy );

		//cout << "Ingestion " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion << endl;
		//cout << "Total drug Ingestion " << Therapy -> Total_Drug_Ingestion << endl;

	}//function
	*/

	/*
		Must be set per drug
	*/
	void Set_Elimination_Constants(unique_ptr<Drug> const & AntiCancer_Agent)
	{
		
		AntiCancer_Agent -> elimintation_constant = -log(0.5) / AntiCancer_Agent -> half_life;
	}

	/*
		Must be set per drug
	*/
	double Initial_Dosing(unique_ptr<Drug> const & AntiCancer_Agent)
	{
		AntiCancer_Agent -> drug_in_system = AntiCancer_Agent -> absorption_Fraction * (AntiCancer_Agent -> dosage * 1000.0); 
		return AntiCancer_Agent -> drug_in_system;
	}

	double Initial_Concentration(unique_ptr<Drug> const & AntiCancer_Agent, double Plasma_Volume)
	{
		AntiCancer_Agent -> concentration = AntiCancer_Agent -> drug_in_system/Plasma_Volume; 
		return AntiCancer_Agent -> concentration ;
	}

	void print_Drug_Intitial_Dosing(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			cout << "Elimination constant: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> elimintation_constant << endl;
			cout << "Drug: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> DrugName
				 << " In system: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> drug_in_system << endl;
		}
			cout << "Total Drug in System: " << Therapy -> Total_Drug_in_System << endl;
	}// function

	void print_Drug_Intitial_Concentration( unique_ptr<Treatment> const & Therapy )
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			cout << "Drug: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> DrugName
				 << " Concentration: " << Therapy -> AntiCancer_Agents -> at(ith_drug) -> concentration << endl;
		}
			cout << "Total Drug Concentration in System: " << Therapy -> Total_Drug_Concentration << endl;
	}

	void Inititalise_Drug_Registers(unique_ptr<Drug> const & AntiCancer_Agent)
	{
		get<0>( AntiCancer_Agent -> MR_Register ).resize( AntiCancer_Agent -> time_to_mutation );
		get<0>( AntiCancer_Agent -> DR_Register ).resize( AntiCancer_Agent -> time_to_cell_death );

		get<0>( AntiCancer_Agent -> MR_Register )[0] = (AntiCancer_Agent -> concentration ); // store intial concentration
		get<1>( AntiCancer_Agent -> MR_Register ) = 1; 										 // Start counter in 1 
		get<2>( AntiCancer_Agent -> MR_Register ) = AntiCancer_Agent -> time_to_mutation;	// Initialise comparing time

		get<0>( AntiCancer_Agent -> DR_Register )[0] = ( AntiCancer_Agent -> concentration ); // store intial concentration
		get<1>( AntiCancer_Agent -> DR_Register  ) = 1;										  // Start counter in 1
		get<2>( AntiCancer_Agent -> DR_Register ) = AntiCancer_Agent -> time_to_cell_death;	  // Initialise comparing time

	}

	void Initialise_Drugs(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		double Plasma_Volume = Therapy -> Plasma_Volume;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			Set_Elimination_Constants( Therapy -> AntiCancer_Agents -> at(ith_drug) );
		
			Therapy -> Total_Drug_in_System += Initial_Dosing( Therapy -> AntiCancer_Agents -> at(ith_drug) );

			Therapy -> Total_Drug_Concentration += Initial_Concentration( Therapy -> AntiCancer_Agents -> at(ith_drug),  Plasma_Volume);

			Inititalise_Drug_Registers( Therapy -> AntiCancer_Agents -> at(ith_drug) );
				
		}
		print_Drug_Intitial_Dosing( Therapy);
		print_Drug_Intitial_Concentration( Therapy);


	}// Function
	//Set initial parameters
	//Set initial dose
	//Modify the model in terms of multiple drugs

	void print_Treatment_Evolution(unique_ptr<Clonal_Expansion> const & CE, unique_ptr<Treatment> const & Therapy, unsigned int hours, unsigned int years, int myID, bool each_100 = false)
	{
		
		if(each_100 && (hours % 100 == 0))
		{
			usleep( (myID*1) + 10 );
			cout << " A C " <<  CE -> Population_Size  
  				 << " CLONES " << CE -> Tumour -> size() 
  			 	<< "   H: " << hours  
  			 	<< " Y: " << years 
  			 	<< " FD: " << CE -> feedback 
  			 	<< " EL: " << Therapy -> Total_Drug_Eliminated
			 	<< " DINS: " << Therapy -> Total_Drug_in_System
			 	<< " CONC: "  << Therapy -> Total_Drug_Concentration
			 	<< " ING: "  << Therapy -> Total_Drug_Ingestion
			 	<< " ID: " << myID
  			 	<< endl;
  			 }
  		else if(each_100 == false)
  			cout << " A C " <<  CE -> Population_Size  
  				 << " CLONES " << CE -> Tumour -> size() 
  			 	<< "   H: " << hours  
  			 	<< " Y: " << years 
  			 	<< " FD: " << CE -> feedback 
  			 	<< " EL: " << Therapy -> Total_Drug_Eliminated
			 	<< " DINS: " << Therapy -> Total_Drug_in_System
			 	<< " CONC: "  << Therapy -> Total_Drug_Concentration
			 	<< " ING: "  << Therapy -> Total_Drug_Ingestion
			 	<< " ID: " << myID
  			 	<< endl;
	}

	double update_drug_elimination( unique_ptr<Drug> const & AntiCancer_Agent )
	{
		AntiCancer_Agent -> eliminated = AntiCancer_Agent -> elimintation_constant * AntiCancer_Agent -> drug_in_system;
		return AntiCancer_Agent -> eliminated;

	}

	double update_drug_in_system( unique_ptr<Drug> const & AntiCancer_Agent )
	{
		AntiCancer_Agent -> drug_in_system = AntiCancer_Agent -> drug_in_system + AntiCancer_Agent -> ingestion - AntiCancer_Agent -> eliminated;
		return AntiCancer_Agent -> drug_in_system;
	}

	double update_drug_concentration(unique_ptr<Drug> const & AntiCancer_Agent, double Plasma_Volume)
	{
		AntiCancer_Agent -> concentration = AntiCancer_Agent -> drug_in_system/Plasma_Volume;
		return AntiCancer_Agent -> concentration;
	}


	void update_drug_parmeters( unique_ptr<Treatment> const & Therapy )
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;
		double Plasma_Volume = Therapy -> Plasma_Volume;

		/*
			Not sure about this
		*/
			Therapy -> Total_Drug_Eliminated = 0.0;
			Therapy -> Total_Drug_in_System = 0.0;
			Therapy -> Total_Drug_Concentration = 0.0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			Therapy -> Total_Drug_Eliminated += update_drug_elimination( Therapy -> AntiCancer_Agents -> at(ith_drug) );
			Therapy -> Total_Drug_in_System += update_drug_in_system( Therapy -> AntiCancer_Agents -> at(ith_drug) );
			Therapy -> Total_Drug_Concentration += update_drug_concentration (Therapy -> AntiCancer_Agents -> at(ith_drug), Plasma_Volume );
		}
	}

	

	bool check_intake_day(unique_ptr<Drug> const & AntiCancer_Agent, unsigned int hours_counter)
	{

		bool intake = false;
		for (vector<unsigned int>::const_iterator i = AntiCancer_Agent -> intake_days_in_hours -> begin(); i != AntiCancer_Agent -> intake_days_in_hours -> end(); ++i)
		{
			if( hours_counter == *i)
			{
				intake = true;
				break;
			}
		}
		return intake;
	}// End model



	void Drug_Diffusion(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;

		/*Check*/
		Therapy -> Total_Drug_Ingestion = 0.0;
		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			if( check_intake_day( Therapy -> AntiCancer_Agents -> at(ith_drug), Therapy -> hours_counter   ) )
			{
				Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion = Therapy -> AntiCancer_Agents -> at(ith_drug) -> absorption_Fraction * Therapy -> AntiCancer_Agents -> at(ith_drug) -> dosage * 1000.0;
			}
			else
			{
				Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion = 0;
			}
			Therapy -> Total_Drug_Ingestion += Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion;
				
		}

	}//Function

	void check_cycle_termination(unique_ptr<Treatment> const & Therapy)
	{
		
		if(Therapy -> Repetition_Day_in_Hours == Therapy -> hours_counter)
		{
			Therapy -> hours_counter = 0;
			Therapy -> Cycles--;
		}
		
	}//

	void Drug_Scheme(unique_ptr<Treatment> const & Therapy)
	{
		if( (Therapy -> Cycles) != 0 )
		{
			check_cycle_termination( Therapy );
			Drug_Diffusion( Therapy );
		}
		else
		{
			unsigned int size = Therapy -> Number_of_Drugs;
			unsigned int ith_drug = 0;
			Therapy -> Total_Drug_Ingestion = 0.0;

			for(ith_drug=0; ith_drug < size; ith_drug++)
				Therapy -> AntiCancer_Agents -> at(ith_drug) -> ingestion = 0;

			Therapy -> hours_to_Exit--;
			if(Therapy -> hours_to_Exit == 0)
			{
				Therapy -> Trial_Done = true;
			
			}
							
		}
	}

	unsigned int Save_Growth_Curve(unique_ptr<Clonal_Expansion> const & CE, unique_ptr<Treatment> const & Therapy , unsigned int elapsed_hours, ofstream &Tumour_Evolution )
	{
		Tumour_Evolution << CE -> Population_Size << "\t" << elapsed_hours << "\t" << Therapy -> Total_Drug_Concentration << "\n";
		return elapsed_hours + 1;
	}

	void Update_Drug_MR_Registers(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;

		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			if( get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ) == get<2>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ) )
			{
				get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ) = 0;
				get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ) = true;
					
				//cout << "MR: " << get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ).at(0) << " " << get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ).at(1) << endl;
				//vector<double>::size_type sz = get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register).size();
				//cout << "\n REGISTER MR FULL size : "  << sz  << " [ ";
				//for (vector<double>::const_iterator i = get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ).begin(); i != get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register ).end(); ++i)
    			//	cout << *i << ' ';
    			//cout << " ]\n"<< endl;
    			//getchar();
    				

			}
			get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register )[get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register )] = Therapy -> AntiCancer_Agents -> at(ith_drug) -> concentration;
			get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> MR_Register)++;
		}//For
	}

	void Update_Drug_DR_Registers(unique_ptr<Treatment> const & Therapy)
	{
		unsigned int size = Therapy -> Number_of_Drugs;
		unsigned int ith_drug = 0;

		for(ith_drug=0; ith_drug < size; ith_drug++)
		{
			if( get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ) == get<2>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ) )
			{
				get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ) = 0;
				get<3>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ) = true;
					
				//cout << "DR: " << get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ).at(0) << " " << get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ).at(1) << endl;
				//vector<double>::size_type sz = get<0>(Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register).size();
				//cout << "\n REGISTER DR FULL size : "  << sz  << " [ ";
				//for (vector<double>::const_iterator i = get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ).begin(); i != get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register ).end(); ++i)
    			//	cout << *i << ' ';
    			//cout << " ]\n"<< endl;
    			//getchar();
    				

			}
			get<0>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register )[get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register )] = Therapy -> AntiCancer_Agents -> at(ith_drug) -> concentration;
			get<1>( Therapy -> AntiCancer_Agents -> at(ith_drug) -> DR_Register)++;
		}//For
	}




	void treatment(unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID, string CTX_file)
	{
	
		
		print_Clonal_Expansion_DS(CE);
		
		unsigned int seconds = 0;
		unsigned int hours = 0;
		unsigned int years = 0;
		unsigned int elapsed_hours = 0;

		string Growth = BasePath + "/Growth/" + "Treatment_Growth.txt";
		string Stats = BasePath + "/Stats/" + "All_Clones_Posterior.txt";
		//unsigned int times_to_wait = 0;
	
		//unsigned int elapsed_hours = 0;
		//unsigned int elapsed_hours = 0;

		unique_ptr<Treatment> const Therapy = get_Treatment_DS();
	
		Read_CTX(CTX_file, Therapy);

		//unique_ptr<Drug> const CTX = get_Drug();


		Initialise_Drugs( Therapy );
		print_Treatment_Evolution( CE,  Therapy,  hours,  years, myID, true );

		//getchar();
		

		ofstream Tumour_Evolution;
		ofstream Pop_Stats;
		Tumour_Evolution.open (Growth);



		//while( (hours < 1600) && !(CTX -> exit) && (CE -> Population_Size < 1000000000 ) )
		while( (years < 6) && !(Therapy -> Trial_Done) )
		{
			seconds += dt;
			if(seconds == 3600)
			{
				seconds = 0; hours ++; Therapy -> hours_counter++;

				update_population_treatment( CE, Therapy, hours, years);
				
				Size_Dependent_Penalty( CE );
				
				if( !(Therapy -> Trial_Done) )
				{

					update_drug_parmeters( Therapy );
					Update_Drug_MR_Registers( Therapy );
					Update_Drug_DR_Registers( Therapy );

				}
				
				if(PRINT)
					print_Treatment_Evolution( CE,  Therapy,  hours,  years, myID, true );

				
				//update_drug_registers( CTX );
				elapsed_hours = Save_Growth_Curve(CE, Therapy ,elapsed_hours, Tumour_Evolution ); // recursive function

				
			}//seconds
			if( !(Therapy -> Trial_Done) )
				Drug_Scheme( Therapy );
			


			if(hours == 8764)//8764
			{
				hours = 0; years ++;
			}
		}//while

		Tumour_Evolution.close();
		Store_Population_Stats(CE, Stats );

	}//treatment end

	bool ToBool( const std::string & s ) 
	{
  		return s.at(0) == '1';
	}

	void add_Clone_to_DS(unique_ptr<Clonal_Expansion> const & CE, vector<string> tokens)
	{
		CE -> Tumour -> push_back( get_Clone_DS() );
		CE -> Tumour -> back() -> Number_of_Memebers_to_Start_Heterogeneity = stoul(tokens[1], nullptr,0);
		CE -> Tumour -> back() -> Generation_ID_Counter = stoul(tokens[2], nullptr,0);
		CE -> Tumour -> back() -> clone_extinct = ToBool(tokens[3]);
		CE -> Tumour -> back() -> Mutation_Rate = stod(tokens[4]);
		CE -> Tumour -> back() -> Number_of_Mutations = stoull(tokens[5]);
		CE -> Tumour -> back() -> Clone_Size = stoul(tokens[6]);
		CE -> Tumour -> back() -> Initiall_Expasion_Period = ToBool(tokens[7]);
		CE -> Tumour -> back() -> In_G0_Phase = ToBool(tokens[8]);
		CE -> Tumour -> back() -> In_G1_Phase = ToBool(tokens[9]);
		CE -> Tumour -> back() -> In_S_Phase = ToBool(tokens[10]);
		CE -> Tumour -> back() -> In_G2_Phase = ToBool(tokens[11]);
		CE -> Tumour -> back() -> In_M_Phase = ToBool(tokens[12]);
		CE -> Tumour -> back() -> P_Expansion[0] = stod(tokens[13]);
		CE -> Tumour -> back() -> P_Expansion[1] = stod(tokens[14]);
		CE -> Tumour -> back() -> P_Expansion[2] = 1 - (CE -> Tumour -> back() -> P_Expansion[0]+ CE -> Tumour -> back() -> P_Expansion[1]);
		CE -> Tumour -> back() -> Remaining_Time_in_G1_Phase = stoul(tokens[15]);
		CE -> Tumour -> back() -> Remaining_Time_in_S_Phase = stoul(tokens[16]);
		CE -> Tumour -> back() -> Remaining_Time_in_G2_Phase = stoul(tokens[17]);
		CE -> Tumour -> back() -> Remaining_Time_in_M_Phase = stoul(tokens[18]);
		CE -> Tumour -> back() -> Generation_ID = tokens[19];
		//print_Clone_DS(CE -> Tumour -> back());
	}

	void open_Tumour_Population(string path, unique_ptr<Clonal_Expansion> const & CE)
	{
		cout << "READING FILE: " << path << endl;
		ifstream infile(path);
		string line;

		int i = 0;
		
		while (getline(infile, line)) // Reading line by line
		{
    		istringstream iss(line); // define the stream to split the content
    		//cout << line << endl;
    		vector<string> tokens{
    								istream_iterator<string>{iss}, 
                	      			istream_iterator<string>{} 
                    	  		}; // Store the string in two different positions of the vector
            if(i != 0 )
            {
            	add_Clone_to_DS(CE,  tokens);
            	/*
            	cout << tokens[0] << " " 
            		 << tokens[1] << " "
            	 	<< tokens[2] << " " 
            	 	<< tokens[3] << " "
            	 	<< tokens[4] << " "
            	 	<< tokens[5] << " "
            	 	<< tokens[6] << " " 
            	 	<< tokens[7] << " "
            	 	<< tokens[8] << " " 
            	 	<< tokens[9] << " "
            	 	<< tokens[10] << " "
            	 	<< tokens[11] << " "
            		 << tokens[12] << " "
            	 	<< tokens[13] << " "
            	 	<< tokens[14] << " "
            	 	<< tokens[15] << " "
            	 	<< tokens[16] << " "
            		 << tokens[17] << " "
            		 << tokens[18] << " "
            		 << tokens[19] 
            		 << endl;
            		 */
			}
			i++;



		}//end while

	}



	void Store_ALL_Population_Stats(unique_ptr<Clonal_Expansion> const & CE, string path)
	{
		unsigned int ith_clone = 0;
		ofstream Pop_Stats;
		Pop_Stats.open (path);
		Pop_Stats << "Number_of_Memebers_to_Start_Heterogeneity\t" 
				  << "Generation_ID_Counter\t" 
				  << "clone_extinct\t" 
				  << "Mutation_Rate\t" 
				  << "Number_of_Mutations\t" 
				  << "Clone_Size\t" 
				  << "Initiall_Expasion_Period\t" 
				  << "In_G0_Phase\t"
				  << "In_G1_Phase\t"
				  << "In_S_Phase\t"
				  << "In_G2_Phase\t"
				  << "In_M_Phase\t"
				  << "DEATH_RATE\t"
				  << "PROLIFERATION_RATE\t"
				  << "Remaining_Time_in_G1_Phase\t"
				  << "Remaining_Time_in_S_Phase\t"
				  << "Remaining_Time_in_G2_Phase\t"
				  << "Remaining_Time_in_M_Phase\t"
				  << "Generation_ID"
				  <<"\n";
  		
  		for( ith_clone = 0; ith_clone < CE-> Tumour -> size() ; ith_clone ++)
  		{
  			Pop_Stats << ith_clone << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> Number_of_Memebers_to_Start_Heterogeneity << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Generation_ID_Counter << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> clone_extinct << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Mutation_Rate << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Number_of_Mutations << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Clone_Size << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> Initiall_Expasion_Period << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> In_G0_Phase << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> In_G1_Phase << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> In_S_Phase << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> In_G2_Phase << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> In_M_Phase << "\t" 
  					  << CE -> Tumour -> at(ith_clone) -> P_Expansion[0] << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> P_Expansion[1] << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Remaining_Time_in_G1_Phase << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Remaining_Time_in_S_Phase << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Remaining_Time_in_G2_Phase << "\t"
  					  << CE -> Tumour -> at(ith_clone) -> Remaining_Time_in_M_Phase << "\t"
  					  << CE -> Tumour -> at (ith_clone) -> Generation_ID << "\t"
  					  << "\n" ;

  		}

  		Pop_Stats.close();
	}

	string get_Working_Dir()
	{
		char temp[500];
   		return ( getcwd(temp, 500) ? string( temp ) : string("") );
	}

	string set_Working_Dir()
	{
		return ( get_Working_Dir() + "/results/") ;
	}

	void create_Dir(string path)
	{
		struct stat st = {0};
		if (stat(path.c_str(), &st) == -1) 
    		mkdir(path.c_str(), 0700);
	}

	void create_Subflders(string path)
	{
		path = path +"/";
		string Stats = path +"Stats";
		string Growth = path +"Growth";
		struct stat st = {0};

		if (stat(Stats.c_str(), &st) == -1) 
    		mkdir(Stats.c_str(), 0700);

    	if (stat(Growth.c_str(), &st) == -1) 
    		mkdir(Growth.c_str(), 0700);

		if(SAVE_CLONAL_EVOLUTION)
		{
			string Clonal_Evolution = path +"Clonal_Evolution";
			if (stat(Clonal_Evolution.c_str(), &st) == -1) 
    			mkdir(Clonal_Evolution.c_str(), 0700);
		}
	}//function

	string generate_TimeStamp()
	{
		time_t rawtime;
  		struct tm * timeinfo;
  		char buffer [80];

  		time (&rawtime);
  		timeinfo = localtime (&rawtime);

  		strftime (buffer,80,"%a_%F_%H_%M_%S",timeinfo);
 		puts (buffer);

 		string timeStamp(buffer);
 		return(timeStamp);
	}

	string create_Store_Directories()
	{
		
		string BasePath = set_Working_Dir() ;
		create_Dir(BasePath);

		string SimulationPath = BasePath + generate_TimeStamp();
		create_Dir(SimulationPath);

 		cout << "Global path @ " << SimulationPath << endl;
 	
 		return SimulationPath;

	}// funcrtion

	string create_myID_Folder(string SimulationPath, int myID)
	{
		string ID_Path = SimulationPath + "/ID_" + to_string(myID);
		create_Dir(ID_Path);
		return (ID_Path + "/");

	}

	void send(string const& str, int dest, int tag, MPI_Comm comm)
	{
		unsigned len = str.size();
		vector<char> nonconst_str(str.begin(), str.end());
		nonconst_str.push_back('\0');

		MPI_Send(&len, 1, MPI_UNSIGNED, dest, tag, comm);
		if (len != 0){
			MPI_Send(&nonconst_str[0], len, MPI_CHAR, dest, tag, comm);
			cout << "Message sent " << endl;
		}
		//MPI_Send(&PingPongCount, 1, MPI_INT, yourID,0, MPI_COMM_WORLD );
	}

	void recv(string& str, int src, int tag, MPI_Comm comm)
	{
		unsigned len;
		MPI_Status s;
		cout << "Sending  " <<endl;
		MPI_Recv(&len, 1, MPI_UNSIGNED, src, tag, comm, &s);
		cout << "RECV len:  " << len << endl;

		if(len != 0)
		{
			vector<char> tmp(len);
			MPI_Recv(tmp.data(), len, MPI_CHAR, src, tag, comm, &s);
			str.assign(tmp.begin(), tmp.end());
		}
		else
		{
			str.clear();
		}
		cout << "DATA RECV: " << str << endl;
	}

	string Path_Bcast_From_Master( MPI_Comm comm)
	{
		char path[1024];
		string tmp = create_Store_Directories();
    	strncpy(path, tmp.c_str(), sizeof(path));
    	path[sizeof(path) - 1] = 0;

		int pathLength = sizeof(path);
		MPI_Bcast (&pathLength, 1, MPI_INT, 0, comm);
		MPI_Bcast (path, pathLength, MPI_CHAR, 0, comm);
		cout << " Path From Master sent" << endl;
		return tmp; 
	}

	string Path_Bcast_From_Salves(int myID, MPI_Comm comm)
	{
		unsigned pathLength;
		char * Path;
		MPI_Bcast (&pathLength, 1, MPI_INT, 0, comm);
   		Path = (char *) malloc (pathLength);
   		MPI_Bcast (Path, pathLength, MPI_CHAR, 0, comm);
    	printf ("Process %d: %s\n", myID, Path);
    	string str(Path);
    	return Path;
	}

	string generate_subfolders(string BasePath, int myID)
	{
		BasePath = BasePath +"/ID_" + to_string(myID);
        create_Dir(BasePath);
        create_Subflders(BasePath);
        return BasePath;
	}

	string ReplaceAll(string str, const string& from, const string& to) 
	{
    	size_t start_pos = 0;
    	while((start_pos = str.find(from, start_pos)) != string::npos) 
    	{
        	str.replace(start_pos, from.length(), to);
        	start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    	}
    	return str;
	}

	void Multiple_Tumour_Scheme(unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID)
	{
			compute_Tumour_Evolution( CE,  BasePath, myID );
        	//MPI_Barrier(MPI_COMM_WORLD);
        	Purge_Death_Clones( CE );
        	//MPI_Barrier(MPI_COMM_WORLD);
        	treatment( CE, BasePath, myID, DEFAULT_TREATMENT_FILE);
        	cout << "Porcess " << myID << " DONE " <<endl;
        	//MPI_Barrier(MPI_COMM_WORLD);
        	MPI_Finalize();
        	//Decide to save alive population
	}

	void set_Expansion_Struct_Parameters(unique_ptr<Clonal_Expansion> const & CE)
	{
		unsigned int ith_clone = 0;
		unsigned int size = CE -> Tumour  -> size();
		unsigned long long int Population_Size = 0;
		
		for (ith_clone = 0; ith_clone < size ; ith_clone++)
			if( !(CE -> Tumour -> at (ith_clone) -> clone_extinct) )
				Population_Size += (unsigned long long int) CE -> Tumour -> at(ith_clone) -> Clone_Size ;
			
								
		CE -> Population_Size = Population_Size;
		Size_Dependent_Penalty(CE);
	}

	void Tumour_Growth(unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID)
	{

		compute_Tumour_Evolution( CE,  BasePath, myID );
        //MPI_Barrier(MPI_COMM_WORLD);
        Purge_Death_Clones( CE );
        Store_ALL_Population_Stats(CE, BasePath +"/Stats/Alive_Clones_Prior.txt");

	}

	void load_Tumour_Population( unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID )
	{
		string temp = ReplaceAll( BasePath, ("ID_"+to_string(myID)), "ID_0") + "/Stats/Alive_Clones_Prior.txt";
        open_Tumour_Population(temp, CE);
        set_Expansion_Struct_Parameters(CE);

	}

	void Drug_Trial( unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID )
	{
		string drug_file = "./DRUG/CTX_Scheme_ID_" + to_string(myID) + ".drug";
		treatment( CE, BasePath, myID, drug_file );
        cout << "Porcess " << myID << " DONE WITH TREATMENT " <<endl;

        //cout << "Porcess " << myID << " PRINT CLONE " <<endl;
        //sleep((unsigned int)myID);
        //for( unsigned int ith_clone = 0; ith_clone < CE-> Tumour -> size() ; ith_clone ++)
        //{
        //	usleep((unsigned int)myID+200);
        //	cout << "Porcess " << myID << " PRINT CLONE " <<endl;
        //	print_Clone_DS(CE -> Tumour -> at(ith_clone));   
        //}   
        //MPI_Barrier(MPI_COMM_WORLD); 	
        
	}

	void Simulate_Tumour_Evolution( unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID )
	{

		if(MULTIPLE_TUMOURS)
        	Multiple_Tumour_Scheme(CE, BasePath, myID);
        else
        {
        	if(myID == 0)
        		Tumour_Growth( CE, BasePath, myID );

        	MPI_Barrier(MPI_COMM_WORLD);
        	
        	if(myID != 0)
        		load_Tumour_Population( CE, BasePath, myID );

        	MPI_Barrier(MPI_COMM_WORLD);
        	Drug_Trial( CE, BasePath, myID );
        
        	MPI_Finalize();
        }//else

	}//Function

	void Simulate_StandAlone_Tumour_Evolution( unique_ptr<Clonal_Expansion> const & CE, string BasePath, int myID )
	{
		cout << " JUST ONE PROCES " << endl;
		// This always should be zero
		if(myID == 0)
        {
        	compute_Tumour_Evolution( CE, BasePath, myID );
        	Purge_Death_Clones( CE );
        	Drug_Trial( CE, BasePath, myID );
        }
        MPI_Finalize();
	}//Function



}  // namespace blah




int main( int argc, char** argv )
{
	using namespace std;
	using namespace core;

	long seed;

	int			myID;
	int			N_Procs; 


 	string BasePath;

	r_global = gsl_rng_alloc (gsl_rng_rand48);     // pick random number generator
  	seed = time (NULL) * getpid();    
  	gsl_rng_set (r_global, seed);  

  	unique_ptr<Clonal_Expansion> const CE = get_Clonnal_Expansion_DS();

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,	&myID); 
	MPI_Comm_size(MPI_COMM_WORLD,	&N_Procs);



  	if(myID == 0)
  		print_Clonal_Expansion_DS(CE);

    try
    {
        if(myID == 0)
        	BasePath = Path_Bcast_From_Master( MPI_COMM_WORLD );
        else
        	BasePath = Path_Bcast_From_Salves( myID, MPI_COMM_WORLD);
        	
        /*
			Wait all processes to enter the follwing section of code
        */
        MPI_Barrier(MPI_COMM_WORLD);
        BasePath = generate_subfolders( BasePath,  myID );
        MPI_Barrier(MPI_COMM_WORLD);// This Barrier may be omitted

        if(N_Procs > 1)
        	Simulate_Tumour_Evolution( CE, BasePath, myID );
        else
        	Simulate_StandAlone_Tumour_Evolution( CE, BasePath,  myID );

        return EXIT_SUCCESS;

    }//try code

    catch( exception const& error )
    {
        cerr << "!" << error.what() << endl;
    }
    gsl_rng_free (r_global);

    
    MPI_Finalize();
	return EXIT_SUCCESS;

}

