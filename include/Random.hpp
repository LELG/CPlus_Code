#ifndef MY_RANDOM
#define MY_RANDOM

#include <random>

namespace core
{
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
}

#endif
