#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

//includes para el algoritmo k_shortest_paths
#include <limits>
#include <set>
#include <map>
#include <queue>

const int NUM_NODOS = 24;
const int K = 5; // max num routes

std::vector<std::vector<int>> read_csv_toInteger(std::string filename);
std::vector<std::vector<std::string>> read_csv(std::string filename);

typedef int CoordenadaIzq;
typedef int CoordenadaDer;
typedef int OD_Demand;


typedef struct 
{
    int numNodo;
    std::pair<CoordenadaIzq, CoordenadaDer> nodo;
} Nodo;

typedef struct
{
    int init_node;
    int term_node;
    double capacity;
    double lenght; 
    double freeFlowTime; 
    int power;
} Carretera;

class NetworkData {
public: 
    NetworkData(std::vector<Nodo> n, std::vector<Carretera> c);
    std::vector<Nodo> getNodos() const;
    std::vector<Carretera> getCarreteras() const;
    int getNumCarreteras() const;
    int getNumNodos() const;

private:
    std::vector<Nodo> nodos;
    std::vector<Carretera> carreteras;
};

// Create Network
NetworkData createNetwork();
std::vector<Nodo> crearListaNodos(const std::vector<std::vector<std::string>> & datosNodos);
std::vector<Carretera> crearListaCarreteras(const std::vector<std::vector<std::string>> & datosCarreteras);

// Compute Strategy Vectors
std::vector<std::vector<OD_Demand>> createOD_Demands();
std::vector<std::vector<std::vector<int>>> computeStrategyVectors(const NetworkData& network, std::vector<std::vector<OD_Demand>>& od_Demands, std::vector<std::pair<int, int>>& od_Pairs, int numRoutes, int multFactor);
std::vector<std::vector<int>> k_shortest_paths(const NetworkData& network, const int& init_node, const int& term_node, const int& k_paths);

// Compute Travel Times
std::vector<double> Compute_traveltimes(const NetworkData& networkData, const std::vector<std::vector<std::vector<int>>>& Strategy_vectors, const std::vector<int>& played_actions, int player_id, std::vector<double> Capacities);

// Auxiliares
int getidx(const NetworkData& network, int nodo1, int nodo2);
int dot_product(std::vector<int> vec1, std::vector<int> vec2);
double dot_product2(std::vector<int> vec1, std::vector<double> vec2);

// Getters
std::vector<double> getCapacities(const NetworkData &n);
std::vector<int> getFreeFlowTimes(const NetworkData& n);
int get_edge_idx(std::vector<Carretera>);

void printNodos(std::vector<Nodo> n);
void printCarreteras(std::vector<Carretera> c);
void printOD_Demands(std::vector<std::vector<OD_Demand>> d);