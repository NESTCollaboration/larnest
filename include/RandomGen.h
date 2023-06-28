/**
 * @file    RandomGen.h
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   Random generator for LArNEST. Based on the NEST RandomGen 
 *          created by Jacob Cutter and edited by others.
 * @version 0.1
 * @date 2023-06-27
 */
#pragma once
#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

namespace larnest
{
    class RandomGen 
    {
    public:
        static RandomGen* GetRandomGen();

        void SetSeed(uint64_t seed);

        int64_t rand_int(int64_t low, int64_t high);
        double rand_uniform();
        double rand_uniform(double low, double high);
        double rand_gaussian(double mean, double sigma);
        double rand_trunc_gaussian(double mean, double sigma, double low, double high);
        double rand_zero_trunc_gaussian(double mean, double sigma);
        int64_t rand_poisson(double mean);
        int64_t rand_binomial(int64_t N0, double prob);
        double rand_exponential(double alpha);

    private:
        RandomGen();
        static RandomGen* sInstance;

        uint64_t sSeed;
        std::random_device sRandomDevice;
        std::mt19937_64 sGenerator;
    };
}