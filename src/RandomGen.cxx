/**
 * @file RandomGen.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-06-27
 */
#include "RandomGen.h"

namespace larnest
{
    RandomGen* RandomGen::sInstance = nullptr;

    RandomGen::RandomGen()
    {
        sSeed = sRandomDevice();
        sGenerator.seed(sSeed);
    }

    RandomGen* RandomGen::GetRandomGen() 
    {
        if (!sInstance) {
            sInstance = new RandomGen;
        }
        return sInstance;
    }

    void RandomGen::SetSeed(uint64_t seed) 
    {
        sSeed = seed;
        sGenerator.seed(seed);
    }

    int64_t RandomGen::rand_int(int64_t low, int64_t high)
    {
        std::uniform_int_distribution<int> uniform_int(low, high);
        return uniform_int(sGenerator);
    }

    double RandomGen::rand_uniform() 
    {
        return rand_uniform(0.0, 1.0);
    }

    double RandomGen::rand_uniform(double low, double high)
    {
        std::uniform_real_distribution<double> uniform_real(low, high);
        return uniform_real(sGenerator);
    }

    double RandomGen::rand_gaussian(double mean, double sigma) 
    {
        std::normal_distribution<double> normal(mean, sigma);
        return normal(sGenerator);
    }

    double RandomGen::rand_trunc_gaussian(double mean, double sigma, double low, double high) 
    {
        std::normal_distribution<double> normal(mean, sigma);
        double value = normal(sGenerator);
        while (value < low || value > high) {
            value = normal(sGenerator);
        }
        return value;
    }

    double RandomGen::rand_zero_trunc_gaussian(double mean, double sigma)
    {
        std::normal_distribution<double> normal(mean, sigma);
        double value = normal(sGenerator);
        while (value < 0.0) {
            value = normal(sGenerator);
        }
        return value;
    }

    int64_t RandomGen::rand_poisson(double mean)
    {
        std::poisson_distribution<int> poisson(mean);
        return poisson(sGenerator);
    }

    int64_t RandomGen::rand_binomial(int64_t N0, double prob)
    {
        std::binomial_distribution<int> binomial(N0, prob);
        return binomial(sGenerator);
    }

    double RandomGen::rand_exponential(double alpha)
    {
        std::exponential_distribution<double> exponential(alpha);
        return exponential(sGenerator);
    }
}
