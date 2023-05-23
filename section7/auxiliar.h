#ifndef auxiliar_h
#define auxiliar_h

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <random> 
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "algorithms.h"
#include "network.h"

class GameData {
public:
    GameData(int N, int T);
    int Simulate_Game(int run, std::vector<Player*>& Players, int T, const NetworkData& network, std::vector<std::vector<std::vector<int>>>& Strategy_vectors, std::vector<double>& sigmas, std::vector<std::vector<double>>& Capacities, std::vector<std::vector<double>>& Total_occupancies, std::vector<std::vector<double>>& addit_Congestions, const std::vector<int>& Contexts, Player* compPlayer, int id, std::vector<std::pair<int, int>>& chosen_arms, std::vector<double>& lossesdef, std::vector<double>& lossescomp);
    const std::vector<std::vector<int>>& getPlayedActions() const { return Played_actions; }
    const std::vector<std::vector<double>>& getIncurredLosses() const { return Incurred_losses; }
    const std::vector<std::vector<double>>& getCumLosses() const { return Cum_losses; }

private:
    std::vector<std::vector<int>> Played_actions;
    std::vector<double> Mixed_strategies;
    std::vector<std::vector<double>> Incurred_losses;
    std::vector<double> Regrets; 
    std::vector<std::vector<double>> Cum_losses;
};

class SiouxNetwork_data_original {};
class Strategy_vectors {};

void Initialize_Players(int N, const std::vector<std::pair<int, int>>& od_Pairs, std::vector<std::vector<std::vector<int>>>& Strategy_vectors, std::vector<double> min_traveltimes, std::vector<double> max_traveltimes, std::vector<int> idxs_controlled, double T, std::string Algo, int version, std::vector<double> Sigma, std::vector<Eigen::MatrixXd>& Kernels, std::vector<double> sigmas, int numberofcontexts, std::vector<std::vector<double>> Capacities, std::vector<Player*>& Players, int& id, Player* &compPlayer);

std::vector<Eigen::MatrixXd> Optimize_Kernels(bool reoptimize, std::string Algo, const std::vector<int>& idxs_controlled, const std::vector<std::vector<std::vector<int>>>& Strategy_vectors, const std::vector<double>& sigmas, int poly_degree, const std::vector<std::vector<double>>& Outcomes, const std::vector<std::vector<double>>& Capacities, const std::vector<std::vector<double>>& Payoffs, std::vector<std::vector<double>>& list_of_param_arrays);

// cargar parametros
std::vector<std::vector<double>> loadParamsFromFile(std::string fileName);

#endif
