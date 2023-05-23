#pragma once
#include <vector>
#include <iostream>
#include <numeric>
#include <random> 
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "network.h"
#include "auxiliar.h"

const int num_routes = 5;
const int multFactor = 1;


class Simulation {

public:
	Simulation(const NetworkData& network);
	void init();

private:
	// jugadores
	int numplayers, controlledplayers;
	std::vector<int> idcontrolledplayers;
	std::string Algo;

	// network
	NetworkData network;
	std::vector<std::vector<OD_Demand>> od_Demands;
	std::vector<std::vector<std::vector<int>>> Strategy_vectors;
	std::vector<std::pair<int, int>> od_Pairs;

	// parametros simulacion
	int rondas, numcontextos;

	// kernel
	int polykernel;
	bool reoptimize;


	void selectParameters();

};