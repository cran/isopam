// Author: Xun Li <lixun910@gmail.com>
// May 27, 2020
//
// Code ported from elki project: http://elki-project.github.io/
// Copyright follows elki project: http://elki-project.github.io/
// AGPLv3: https://github.com/elki-project/elki/blob/master/LICENSE.md
//
// Minimal version for isopam: FastPAM only
// Warnings fixed by isopam maintainers

#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <cfloat>
#include <algorithm>

#include "pam.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// BUILD
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BUILD::BUILD(DistMatrix* dist) : PAMInitializer(dist) {
}

std::vector<int> BUILD::run(const std::vector<int>& ids, int k)
{
    int nn = static_cast<int>(ids.size());
    std::vector<int> medids;
    std::set<int> medids_set;
    
    int ssize = 10 + static_cast<int>(ceil(sqrt(static_cast<double>(nn))));
    if (ssize < nn) ssize = nn;
    
    std::vector<double> mindist(nn), bestd(nn), tempd(nn), temp;
    
    int bestid = -1;  // Fixed: initialize bestid
    {
        double best = DBL_MAX;
        for (int i = 0; i < nn; ++i) {
            double sum = 0, d;
            for (int j = 0; j < nn; ++j) {
                d = dist->getDistance(ids[i], ids[j]);
                sum += d;
                tempd[ids[j]] = d;
            }
            if (sum < best) {
                best = sum;
                bestid = ids[i];
                temp = mindist;
                mindist = tempd;
                tempd = temp;
            }
        }
        medids.push_back(bestid);
        medids_set.insert(bestid);
    }
    
    for (int i = 1; i < k; i++) {
        double best = DBL_MAX;
        bestid = -1;
        for (int j = 0; j < nn; ++j) {
            if (medids_set.find(ids[j]) != medids_set.end()) {
                continue;
            }
            double sum = 0., v;
            for (int kk = 0; kk < nn; ++kk) {
                v = std::min(dist->getDistance(ids[j], ids[kk]), mindist[ids[kk]]);
                sum += v;
                tempd[ids[kk]] = v;
            }
            if (sum < best) {
                best = sum;
                bestid = ids[j];
                temp = bestd;
                bestd = tempd;
                tempd = temp;
            }
        }
        if (bestid == -1) {
            return std::vector<int>();
        }
        medids.push_back(bestid);
        medids_set.insert(bestid);
        temp = bestd;
        bestd = mindist;
        mindist = temp;
    }
    return medids;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// LAB
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LAB::LAB(DistMatrix* dist, int seed) : PAMInitializer(dist)
{
    random.SetSeed(seed);
}

void LAB::shuffle(std::vector<int>& samples, int ssize, int end) {
    ssize = ssize < end ? ssize : end;
    
    for (int i = 1; i < ssize; i++) {
        int a = i - 1, b = i + random.nextInt(end - i);
        int tmp = samples[b];
        samples[b] = samples[a];
        samples[a] = tmp;
    }
}

double LAB::getMinDist(int j, std::vector<int>& medids, std::vector<double>& mindist)
{
    double prev = mindist[j];
    if (prev == DBL_MIN) {
        prev = DBL_MAX;
        for (size_t i = 0; i < medids.size(); ++i) {
            double d = dist->getDistance(j, medids[i]);
            prev = d < prev ? d : prev;
        }
        mindist[j] = prev;
    }
    return prev;
}

std::vector<int> LAB::run(const std::vector<int>& ids, int k)
{
    int nn = static_cast<int>(ids.size());
    std::vector<int> medids;
    std::set<int> medids_set;
    
    int ssize = 10 + static_cast<int>(ceil(sqrt(static_cast<double>(nn))));
    if (ssize > nn) ssize = nn;
    
    std::vector<double> mindist(nn, DBL_MIN), bestd(nn), tempd(nn, DBL_MIN), tmp;
    std::vector<int> sample(nn);
    for (int i = 0; i < nn; ++i) sample[i] = ids[i];
    int range = static_cast<int>(sample.size());
    
    shuffle(sample, ssize, range);
    
    {
        double best = DBL_MAX;
        int bestoff = -1;
        
        for (int i = 0; i < ssize; ++i) {
            double sum = 0, d;
            for (size_t j = 0; j < tempd.size(); ++j) tempd[j] = DBL_MIN;
            for (int j = 0; j < ssize; ++j) {
                d = dist->getDistance(sample[i], sample[j]);
                sum += d;
                tempd[sample[j]] = d;
            }
            if (sum < best) {
                best = sum;
                bestoff = i;
                tmp = mindist;
                mindist = tempd;
                tempd = tmp;
            }
        }
        medids.push_back(sample[bestoff]);
        medids_set.insert(sample[bestoff]);
        int tmp_val = sample[--range];
        sample[range] = sample[bestoff];
        sample[bestoff] = tmp_val;
    }
    
    while (medids.size() < static_cast<size_t>(k)) {
        ssize = range < ssize ? range : ssize;
        shuffle(sample, ssize, range);
        double best = DBL_MAX;
        int bestoff = -1;
        for (int i = 0; i < ssize; ++i) {
            if (medids_set.find(sample[i]) != medids_set.end()) {
                continue;
            }
            double sum = 0., v;
            for (size_t j = 0; j < tempd.size(); ++j) tempd[j] = DBL_MIN;
            for (int j = 0; j < ssize; ++j) {
                double prev = getMinDist(sample[j], medids, mindist);
                if (prev == 0) {
                    continue;
                }
                v = std::min(dist->getDistance(sample[i], sample[j]), prev);
                sum += v;
                tempd[sample[j]] = v;
            }
            if (sum < best) {
                best = sum;
                bestoff = i;
                tmp = bestd;
                bestd = tempd;
                tempd = tmp;
            }
        }
        if (bestoff < 0) {
            return medids;
        }
        medids.push_back(sample[bestoff]);
        medids_set.insert(sample[bestoff]);
        int tmp_val = sample[--range];
        sample[range] = sample[bestoff];
        sample[bestoff] = tmp_val;
        tmp = bestd;
        bestd = mindist;
        mindist = tmp;
    }
    
    return medids;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// PAM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PAM::PAM(int num_obs, DistMatrix* dist_matrix, PAMInitializer* init, 
         int k, int maxiter, const std::vector<int>& _ids)
    : num_obs(num_obs), dist_matrix(dist_matrix), initializer(init),
      k(k), maxiter(maxiter), ids(_ids)
{
    if (initializer == NULL) {
        initializer = new BUILD(dist_matrix);
    }
    
    if (ids.empty()) {
        ids.resize(num_obs);
        for (int i = 0; i < num_obs; ++i) {
            ids[i] = i;
        }
    }
}

PAM::~PAM() {
}

double PAM::run() {
    std::vector<int> seq_ids(num_obs);
    for (int i = 0; i < num_obs; ++i) {
        seq_ids[i] = i;
    }
    
    medoids = initializer->run(seq_ids, k);
    
    assignment.resize(num_obs, -1);
    nearest.resize(num_obs, -1);
    second.resize(num_obs, -1);
    
    double cost = run(medoids, maxiter);
    
    return cost;
}

std::vector<int> PAM::getMedoids()
{
    return medoids;
}

std::vector<int> PAM::getResults()
{
    std::vector<int> cluster_result(num_obs, 0);
    for (int i = 0; i < num_obs; ++i) {
        cluster_result[i] = assignment[ids[i]] + 1;
    }
    return cluster_result;
}

double PAM::run(std::vector<int>& medoids, int maxiter) {
    int k = static_cast<int>(medoids.size());
    double tc = assignToNearestCluster(medoids);
    int bestid = -1;  // Fixed: initialize bestid
    int iteration = 0;
    while (iteration < maxiter || maxiter <= 0) {
        ++iteration;
        double best = DBL_MAX;
        int bestcluster = -1;
        for (int h = 0; h < num_obs; ++h) {
            if (medoids[assignment[h]] == h) {
                continue;
            }
            double hdist = nearest[h];
            if (hdist <= 0.) {
                continue;
            }
            for (int pi = 0; pi < k; pi++) {
                double cpi = computeReassignmentCost(h, pi) - hdist;
                if (cpi < best) {
                    best = cpi;
                    bestid = h;
                    bestcluster = pi;
                }
            }
        }
        if (!(best < -1e-12 * tc)) {
            break;
        }
        medoids[bestcluster] = bestid;
        double nc = assignToNearestCluster(medoids);
        if (nc > tc) {
            if (nc - tc < 1e-7 * tc) {
                break;
            }
            break;
        }
        tc = nc;
    }
    return tc;
}

double PAM::getDistance(int i, int j) {
    return dist_matrix->getDistance(i, j);
}

double PAM::assignToNearestCluster(std::vector<int>& means) {
    double cost = 0.;
    for (int i = 0; i < num_obs; ++i) {
        double mindist = DBL_MAX, mindist2 = DBL_MAX;
        int minindx = -1;
        for (size_t j = 0; j < means.size(); ++j) {
            double dist = getDistance(i, means[j]);
            if (dist < mindist) {
                mindist2 = mindist;
                minindx = static_cast<int>(j);
                mindist = dist;
            } else if (dist < mindist2) {
                mindist2 = dist;
            }
        }
        if (minindx < 0) {
            return 0;
        }
        assignment[i] = minindx;
        nearest[i] = mindist;
        second[i] = mindist2;
        cost += mindist;
    }
    return cost;
}

double PAM::computeReassignmentCost(int h, int mnum) {
    double cost = 0.;
    for (int j = 0; j < num_obs; ++j) {
        if (h == j) {
            continue;
        }
        double distcur = nearest[j];
        double dist_h = getDistance(h, j);
        if (assignment[j] == mnum) {
            cost += std::min(dist_h, second[j]) - distcur;
        } else if (dist_h < distcur) {
            cost += dist_h - distcur;
        }
    }
    return cost;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// FastPAM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FastPAM::FastPAM(int num_obs, DistMatrix* dist_matrix, PAMInitializer* init,
                 int k, int maxiter, double fasttol, const std::vector<int>& _ids)
    : PAM(num_obs, dist_matrix, init, k, maxiter, _ids), fasttol(fasttol)
{
    fastswap = 1 - fasttol;
}

FastPAM::~FastPAM() {
}

double FastPAM::run(std::vector<int>& medoids, int maxiter) {
    int k = static_cast<int>(medoids.size());
    double tc = assignToNearestCluster(medoids);
    int fastswaps = 0;
    std::vector<int> bestids(k);
    int bestid = -1;  // Fixed: initialize bestid
    
    std::vector<double> best(k), cost(k);
    
    int iteration = 0;
    while (iteration < maxiter || maxiter <= 0) {
        ++iteration;
        findBestSwaps(medoids, bestids, best, cost);
        int min = argmin(best);
        if (!(best[min] < -1e-12 * tc)) {
            break;
        }
        while (min >= 0 && best[min] < -1e-12 * tc) {
            bestid = bestids[min];
            updateAssignment(medoids, bestid, min);
            
            tc += best[min];
            best[min] = DBL_MAX;
            
            while ((min = argmin(best)) >= 0 && best[min] < -1e-12 * tc) {
                bestid = bestids[min];
                if (medoids[assignment[bestid] & 0x7FFF] == bestid) {
                    best[min] = DBL_MAX;
                    continue;
                }
                double hdist = nearest[bestid];
                double c = computeReassignmentCost(bestid, min) - hdist;
                if (c <= best[min] * fastswap) {
                    best[min] = c;
                    ++fastswaps;
                    break;
                }
                best[min] = DBL_MAX;
            }
        }
    }
    
    for (int i = 0; i < num_obs; ++i) {
        assignment[i] = assignment[i] & 0x7FFF;
    }
    return tc;
}

void FastPAM::findBestSwaps(std::vector<int>& medoids, std::vector<int>& bestids,
                            std::vector<double>& best, std::vector<double>& cost)
{
    size_t n_medoids = medoids.size();
    
    best.resize(n_medoids, DBL_MAX);
    cost.resize(n_medoids, 0);
    
    for (int h = 0; h < num_obs; ++h) {
        if (medoids[assignment[h] & 0x7FFF] == h) {
            continue;
        }
        
        for (size_t j = 0; j < n_medoids; ++j) cost[j] = -nearest[h];
        
        computeReassignmentCost(h, cost);
        
        for (size_t i = 0; i < cost.size(); i++) {
            double costi = cost[i];
            if (costi < best[i]) {
                best[i] = costi;
                bestids[i] = h;
            }
        }
    }
}

bool FastPAM::isMedoid(int id) {
    return false;
}

void FastPAM::computeReassignmentCost(int h, std::vector<double>& cost) {
    for (int j = 0; j < num_obs; ++j) {
        if (h == j) {
            continue;
        }
        double distcur = nearest[j];
        double distsec = second[j];
        double dist_h = getDistance(h, j);
        int pj = assignment[j] & 0x7FFF;
        
        cost[pj] += std::min(dist_h, distsec) - distcur;
        if (dist_h < distcur) {
            double delta = dist_h - distcur;
            for (int pi = 0; pi < pj; pi++) {
                cost[pi] += delta;
            }
            for (size_t pi = static_cast<size_t>(pj) + 1; pi < cost.size(); pi++) {
                cost[pi] += delta;
            }
        }
    }
}

double FastPAM::computeReassignmentCost(int h, int mnum)
{
    double cost = 0.;
    for (int j = 0; j < num_obs; ++j) {
        if (h == j) {
            continue;
        }
        double distcur = nearest[j];
        double dist_h = getDistance(h, j);
        if ((assignment[j] & 0x7FFF) == mnum) {
            cost += std::min(dist_h, second[j]) - distcur;
        } else if (dist_h < distcur) {
            cost += dist_h - distcur;
        }
    }
    return cost;
}

double FastPAM::assignToNearestCluster(std::vector<int>& means) {
    double cost = 0.;
    for (int j = 0; j < num_obs; ++j) {
        int iditer = j;
        double mindist = DBL_MAX, mindist2 = DBL_MAX;
        int minindx = -1, minindx2 = -1;
        for (size_t h = 0; h < means.size(); ++h) {
            double dist = getDistance(iditer, means[h]);
            if (dist < mindist) {
                minindx2 = minindx;
                mindist2 = mindist;
                minindx = static_cast<int>(h);
                mindist = dist;
            } else if (dist < mindist2) {
                minindx2 = static_cast<int>(h);
                mindist2 = dist;
            }
        }
        if (minindx < 0) {
            return 0;
        }
        assignment[iditer] = minindx | (minindx2 << 16);
        nearest[iditer] = mindist;
        second[iditer] = mindist2;
        
        cost += mindist;
    }
    return cost;
}

int FastPAM::argmin(const std::vector<double>& best) {
    double min = DBL_MAX;
    int ret = -1;
    for (size_t i = 0; i < best.size(); i++) {
        double v = best[i];
        if (v < min) {
            min = v;
            ret = static_cast<int>(i);
        }
    }
    return ret;
}

void FastPAM::updateAssignment(std::vector<int>& medoids, int h, int m) {
    medoids[m] = h;
    
    double hdist = nearest[h];
    nearest[h] = 0;
    
    int olda = assignment[h];
    
    if ((olda & 0x7FFF) != m) {
        assignment[h] = m | ((olda & 0x7FFF) << 16);
        second[h] = hdist;
    } else {
        assignment[h] = m | (olda & 0x7FFF0000);
    }
    
    for (int i = 0; i < num_obs; ++i) {
        int j = i;
        if (h == j) {
            continue;
        }
        double distcur = nearest[j];
        double distsec = second[j];
        double dist_h = getDistance(h, j);
        int pj = assignment[j];
        int po = static_cast<unsigned int>(pj) >> 16;
        pj &= 0x7FFF;
        if (pj == m) {
            if (dist_h < distsec) {
                nearest[j] = dist_h;
                assignment[j] = m | (po << 16);
            } else {
                nearest[j] = distsec;
                assignment[j] = po | (updateSecondNearest(j, medoids, m, dist_h, po) << 16);
            }
        } else {
            if (dist_h < distcur) {
                nearest[j] = dist_h;
                second[j] = distcur;
                assignment[j] = m | (pj << 16);
            } else if (po == m) {
                assignment[j] = pj | (updateSecondNearest(j, medoids, m, dist_h, pj) << 16);
            } else if (dist_h < distsec) {
                second[j] = dist_h;
                assignment[j] = pj | (m << 16);
            }
        }
    }
}

int FastPAM::updateSecondNearest(int j, std::vector<int>& medoids, int h, double dist_h, int n)
{
    double sdist = dist_h;
    int sbest = h;
    for (size_t i = 0; i < medoids.size(); ++i) {
        if (static_cast<int>(i) != h && static_cast<int>(i) != n) {
            double d = getDistance(j, medoids[i]);
            if (d < sdist) {
                sdist = d;
                sbest = static_cast<int>(i);
            }
        }
    }
    second[j] = sdist;
    return sbest;
}
