#include "Random.hpp"
#include <random>

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

unsigned int Random::binomial_sample(unsigned int n, double p)
{
    return std::binomial_distribution<unsigned int>{n, p}(eng);
}

double Random::Uniform_Mutation_Rate(double mu_rate)
{
    return uniform_real_distribution<double>{mu_rate/2.0, mu_rate*2.0}(eng);
}
double Random::Update_Proliferation_Rate(double Parent_Proliferation_Rate)
{
    double U = normal_distribution<double>{Parent_Proliferation_Rate, 0.001 }(eng);
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

    double U = normal_distribution<double>{Parent_mu_rate, mut_rate }(eng);
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
