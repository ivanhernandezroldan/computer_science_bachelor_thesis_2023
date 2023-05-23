#include "simulation.h"


Simulation::Simulation(const NetworkData& network) : network(network) {
	this->network = network;
	this->od_Demands = createOD_Demands(); 
	this->Strategy_vectors = computeStrategyVectors(network, od_Demands, this->od_Pairs, num_routes, multFactor);
	this->numplayers = this->Strategy_vectors.size();
	selectParameters();
}

void Simulation::selectParameters(){
	// parámetros que se establecen para controlar el juego->por defecto toman estos valores NO MODIFICAR 
	this->controlledplayers = 10;
	this->rondas = 80;
	this->numcontextos = 7;
	this->polykernel = 4;
	this->reoptimize = false;
	this->Algo = "cGPMW";
}

void Simulation::init(){
	// índices de jugadores controlados, se obtienen aleatoriamente
	std::vector<int> aux(this->controlledplayers, 0);
	this->idcontrolledplayers = aux;
	std::vector<int> idxs_all(this->numplayers);
	std::iota(idxs_all.begin(), idxs_all.end(), 0);

	std::random_device rd;
	std::mt19937 g(rd()); 
	std::shuffle(idxs_all.begin(), idxs_all.end(), g);

	std::copy(idxs_all.begin(), idxs_all.begin() + this->controlledplayers, this->idcontrolledplayers.begin()); // elegir los controlledPlayers entre todos los players
	
	// capacidades y contextos aleatorios
	std::vector<double> initCapacities = getCapacities(network);
	std::vector<std::vector<double>> Capacities;
	for (int c = 0; c < this->numcontextos; c++) {
		std::vector<double> perturbed_capacities;
		for (int i = 0; i < initCapacities.size(); i++) {
			double perturbation_factor = 0.7 + 0.5 * (double)std::rand() / RAND_MAX;
			double perturbed_capacity = perturbation_factor * initCapacities[i];
			perturbed_capacities.push_back(perturbed_capacity);
		}
		Capacities.push_back(perturbed_capacities); // unas capacidades para cada contexto
		// 76 capacidades (1 por carretera) distintas para cada contexto (10)
	}
	int runs = 1;
	std::srand(runs); // inicializa el generador de numeros aleatorios --> rand()
	std::vector<int> Contexts(this->rondas); // un contexto de los 10 posibles para cada ronda
	for (int i = 0; i < this->rondas; i++) {
		Contexts[i] = std::rand() % this->numcontextos; // el contexto de la ronda i es el resto de dividir un numero aleatorio entre el numero de contextos
	}


	// escalar recompensas entre 0 1, para ello se calculan maxtraveltimes y mintraveltimes
	std::vector<double> max_traveltimes(this->numplayers, 0.0); 
	std::vector<double> min_traveltimes(this->numplayers, 1e8);
	std::vector<std::vector<double>> Capacities_rand;
	std::vector<std::vector<double>> Outcomes;
	std::vector<std::vector<double>> Payoffs;
	for (int i = 0; i < this->rondas; i++) {
		std::vector<double> outcome(this->numplayers, 0.0); // todos los jugadores aplican la estrategia 0 por defecto
		for (int p : this->idcontrolledplayers) {
			int num_actions = Strategy_vectors[p].size();
			outcome[p] = std::rand() % num_actions; // el outcome va a ser la accion que juege el jugador: un numero entre 1 y 76 carreteras
		}
		std::vector<double> capacities = Capacities[std::rand() % this->numcontextos]; // se pasan las capacidades de un determinado contexto aleatorio
		std::vector<int> outcomeint(outcome.size());
		for (int i = 0; i < outcome.size(); i++) {
			outcomeint[i] = static_cast<int>(outcome[i]);
		}
	
		std::vector<double> traveltimes = Compute_traveltimes(this->network, Strategy_vectors, outcomeint, -1, capacities);
		for (int n = 0; n < this->numplayers; n++) {
			max_traveltimes[n] = std::max(max_traveltimes[n], traveltimes[n] + 0.01);
			min_traveltimes[n] = std::min(min_traveltimes[n], traveltimes[n] - 0.01);
		}
		Capacities_rand.push_back(capacities);
		Outcomes.push_back(outcome);
		std::vector<double> payoffs(this->numplayers);
		for (int n = 0; n < this->numplayers; n++) {
			payoffs[n] = -traveltimes[n];
		}
		Payoffs.push_back(payoffs);
	}

	/*
	Entendemos que los maxtraveltimes y los mintravel times se van a usar en cada algoritmo para escalar los premios entre 0 y 1. Payoff es la recompensa, y, 
	usando la recompensa se escalaran los premios así: pero esto se entiende que se hará dentro del propio algoritmo
	scaled_payoffs = (max_traveltimes - Payoffs) / (max_traveltimes - min_traveltimes)
	*/

	// Kernel
	std::vector<double> sigmas(max_traveltimes.size()); // tamaño = numJugadores --> una desviación típica por cada jugador

	for (int i = 0; i < max_traveltimes.size(); ++i) {
		sigmas[i] = 0.001 * (max_traveltimes[i] - min_traveltimes[i]); // cambiado para ajustar menor valor de sigma
	}

	std::vector<std::vector<double>> list_of_param_arrays;
	std::vector<Eigen::MatrixXd> Kernels (this->numplayers);
	if(this->Algo == "cGPMW" || this->Algo == "GPMW")
		Kernels = Optimize_Kernels(this->reoptimize, this->Algo, this->idcontrolledplayers, this->Strategy_vectors, sigmas, this->polykernel, Outcomes, Capacities, Payoffs, list_of_param_arrays);
	

	// inicializar jugadores
	std::vector<Player*> players(this->numplayers);
	int id = -1;
	Player* compPlayer = nullptr;
	Initialize_Players(this->numplayers, this->od_Pairs, this->Strategy_vectors, min_traveltimes, max_traveltimes, this->idcontrolledplayers, this->rondas, this->Algo, 0, sigmas, Kernels, sigmas, this->numcontextos, Capacities, players, id, compPlayer);
	
	
	// Simulate Game
	int run = 1;
	std::vector<std::vector<double>> addit_Congestions (this->rondas);
	std::vector<std::vector<double>> Total_occupancies;
	GameData game = GameData(this->numplayers, this->rondas);
	std::vector<std::pair<int, int>> chosen_arms(this->rondas);
	std::vector<double> lossesdef(this->rondas);
	std::vector<double> lossescomp(this->rondas);


	int cumLossescmp = game.Simulate_Game(run, players, this->rondas, this->network, this->Strategy_vectors, sigmas, Capacities, Total_occupancies, addit_Congestions, Contexts, compPlayer, id, chosen_arms,lossesdef, lossescomp); // 
	
	for (int p = 0; p < players.size(); p++) {
		if (players[p]->getType() == PlayerType::cGPMW) {
			std::cout << "Jugador cGPMW " << p << ": " << game.getCumLosses()[this->rondas - 1][p] << std::endl;
		}
		else if (players[p]->getType() == PlayerType::GPMW) {
			std::cout << "Jugador GPMW " << p << ": " << game.getCumLosses()[this->rondas - 1][p] << std::endl;
		}
		else {
			std::cout << "Jugador Hedge " << p << ": " << game.getCumLosses()[this->rondas - 1][p] << std::endl;
		}
	}
	
	if (id != -1) {
		double ind = game.getCumLosses()[this->rondas - 1][id] / cumLossescmp;
		ind = 1 - ind;
		ind *= 100;
		std::cout << "-----------Resultados-----------" << std::endl << "En este caso tenemos para el " << id << ":" << std::endl;
		std::cout << "Unas perdidas acumuladas GPMW: " << cumLossescmp << " vs cGPMW: " << game.getCumLosses()[this->rondas - 1][id] << std::endl;
		std::cout << "Estamos reduciendo las perdidas con el bandido contextual un " << ind << "%" << std::endl;

		std::ofstream file;
		file.open("output.txt");

		for (int i = 0; i < this->rondas; i++) {
			file << chosen_arms[i].first << "," << chosen_arms[i].second << ",";
			file << lossesdef[i] << "," << lossescomp[i] << "\n";
		}

		file.close();
	}

}

