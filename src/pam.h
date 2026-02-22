// Author: Xun Li <lixun910@gmail.com>
// May 27, 2020
//
// Code ported from elki project: http://elki-project.github.io/
// Copyright follows elki project: http://elki-project.github.io/
// AGPLv3: https://github.com/elki-project/elki/blob/master/LICENSE.md
//
// Minimal version for isopam: FastPAM only
// Warnings fixed by isopam maintainers

#ifndef __XL_PAM_H
#define __XL_PAM_H

#include <vector>
#include <unordered_map>
#include <cstdlib>

using namespace std;

/**
 * Xoroshiro128+ random number generator
 * http://xoroshiro.di.unimi.it/
 */
class Xoroshiro128Random
{
    signed long long s0;
    signed long long s1;
public:
    Xoroshiro128Random() : s0(0), s1(0) {}
    virtual ~Xoroshiro128Random() {}
    
    void SetSeed(signed long long xor64) {
        xor64 ^= (unsigned long long)xor64 >> 12;
        xor64 ^= xor64 << 25;
        xor64 ^= (unsigned long long)xor64 >> 27;
        s0 = xor64 * 2685821657736338717L;
        xor64 ^= (unsigned long long)xor64 >> 12;
        xor64 ^= xor64 << 25;
        xor64 ^= (unsigned long long)xor64 >> 27;
        s1 = xor64 * 2685821657736338717L;
    }
    
    int nextInt(int n) {
        if (n <= 0) return 0;
        // Fixed: added parentheses around (n - 1)
        int r = (int)((n & -n) == n ? nextLong() & (n - 1)
            : (unsigned long long)(((unsigned long long)nextLong() >> 32) * n) >> 32);
        return r;
    }
    
    signed long long nextLong() {
        signed long long t0 = s0, t1 = s1;
        signed long long result = t0 + t1;
        t1 ^= t0;
        s0 = (t0 << 55) | ((unsigned long long)t0 >> (64 - 55));
        s0 = s0 ^ t1 ^ (t1 << 14);
        s1 = (t1 << 36) | ((unsigned long long)t1 >> (64 - 36));
        return result;
    }
    
    double nextDouble() {
        char tempStr[] = "0x1.0p-53";
        double nd = std::strtod(tempStr, NULL);
        return ((unsigned long long)nextLong() >> 11) * nd;
    }
    
    std::vector<int> randomSample(int samplesize, int n)
    {
        std::vector<int> samples(samplesize);
        int i = 0;
        unordered_map<int, bool> sample_dict;
        // Fixed: cast samplesize to size_t for comparison
        while (sample_dict.size() < static_cast<size_t>(samplesize)) {
            int rnd = nextInt(n);
            if (sample_dict.find(rnd) == sample_dict.end()) {
                samples[i++] = rnd;
            }
            sample_dict[rnd] = true;
        }
        return samples;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Distance Matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DistMatrix
{
protected:
    std::vector<int> ids;
    bool has_ids;
public:
    DistMatrix(const std::vector<int>& _ids = std::vector<int>())
        : ids(_ids), has_ids(!_ids.empty()) {}
    virtual ~DistMatrix() {}
    virtual double getDistance(int i, int j) = 0;
    virtual void setIds(const std::vector<int>& _ids) {
        ids = _ids;
        has_ids = !ids.empty();
    }
};

class RDistMatrix : public DistMatrix
{
    int num_obs;
    int n;
    const std::vector<double>& dist;
public:
    RDistMatrix(int num_obs, const std::vector<double>& dist, 
                const std::vector<int>& _ids = std::vector<int>())
        : DistMatrix(_ids), num_obs(num_obs), dist(dist) {
        n = (num_obs - 1) * num_obs / 2;
    }
    virtual ~RDistMatrix() {}
    
    virtual double getDistance(int i, int j) {
        if (i == j) return 0;
        if (has_ids) {
            i = ids[i];
            j = ids[j];
        }
        int r = i > j ? i : j;
        int c = i < j ? i : j;
        int idx = n - (num_obs - c - 1) * (num_obs - c) / 2 + (r - c) - 1;
        return dist[idx];
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initializer
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PAMInitializer
{
protected:
    DistMatrix* dist;
public:
    PAMInitializer(DistMatrix* dist) : dist(dist) {}
    virtual ~PAMInitializer() {}
    virtual std::vector<int> run(const std::vector<int>& ids, int k) = 0;
};

class BUILD : public PAMInitializer
{
public:
    BUILD(DistMatrix* dist);
    virtual ~BUILD() {}
    virtual std::vector<int> run(const std::vector<int>& ids, int k);
};

class LAB : public PAMInitializer
{
public:
    LAB(DistMatrix* dist, int seed);
    virtual ~LAB() {}
    virtual std::vector<int> run(const std::vector<int>& ids, int k);
    void shuffle(std::vector<int>& samples, int ssize, int end);
    double getMinDist(int j, std::vector<int>& medids, std::vector<double>& mindist);
    Xoroshiro128Random random;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// PAM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PAM
{
public:
    PAM(int num_obs, DistMatrix* dist_matrix, PAMInitializer* init,
        int k, int maxiter, const std::vector<int>& ids = std::vector<int>());
    virtual ~PAM();
    
    virtual double run();
    virtual std::vector<int> getResults();
    virtual std::vector<int> getMedoids();
    virtual std::vector<int> getAssignement() { return assignment; }

protected:
    double getDistance(int i, int j);
    virtual double run(std::vector<int>& medoids, int maxiter);
    virtual double assignToNearestCluster(std::vector<int>& means);
    virtual double computeReassignmentCost(int h, int mnum);

protected:
    int num_obs;
    DistMatrix* dist_matrix;
    PAMInitializer* initializer;
    int k;
    int maxiter;
    std::vector<int> ids;
    std::vector<int> assignment;
    std::vector<double> nearest;
    std::vector<double> second;
    std::vector<int> medoids;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// FastPAM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class FastPAM : public PAM
{
public:
    FastPAM(int num_obs, DistMatrix* dist_matrix, PAMInitializer* init,
            int k, int maxiter, double fasttol, 
            const std::vector<int>& ids = std::vector<int>());
    virtual ~FastPAM();
    
    virtual double run() { return PAM::run(); }

protected:
    virtual double run(std::vector<int>& medoids, int maxiter);
    virtual void computeReassignmentCost(int h, std::vector<double>& cost);
    virtual double assignToNearestCluster(std::vector<int>& means);
    virtual double computeReassignmentCost(int h, int mnum);
    
    void findBestSwaps(std::vector<int>& medoids,
                       std::vector<int>& bestids,
                       std::vector<double>& best,
                       std::vector<double>& cost);
    
    bool isMedoid(int id);
    int argmin(const std::vector<double>& best);
    void updateAssignment(std::vector<int>& medoids, int h, int m);
    int updateSecondNearest(int j, std::vector<int>& medoids,
                            int h, double dist_h, int n);

protected:
    double fastswap;
    double fasttol;
};

#endif
