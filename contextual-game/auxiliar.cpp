#include "auxiliar.h"

GameData:: GameData(int N, int T){
    this->Played_actions = std::vector<std::vector<int>>(T);
    this->Mixed_strategies = std::vector<double>(N); 
    this->Incurred_losses = std::vector<std::vector<double>>(T);
    this->Regrets = std::vector<double>(N);
    this->Cum_losses = std::vector<std::vector<double>>(T);
}

int GameData::Simulate_Game(int run, std::vector<Player*>& Players, int T, const NetworkData& network, std::vector<std::vector<std::vector<int>>>& Strategy_vectors, std::vector<double>& sigmas, std::vector<std::vector<double>>& Capacities, std::vector<std::vector<double>>& Total_occupancies, std::vector<std::vector<double>>& addit_Congestions, const std::vector<int>& Contexts, Player* compPlayer, int id, std::vector<std::pair<int, int>> &chosen_arms, std::vector<double> &lossesdef, std::vector<double> &lossescomp)
{
    int incurredlosses = 0;
    bool haycomp = true;
    std::vector<double> cumlosscomp(T);
    std::vector<double> cumlossdef(T);
    cumlosscomp[0] = 0;
    cumlossdef[0] = 0;

    
    if (id == -1) haycomp = false;
    int N = Players.size();
    std::vector<int> playedactions(T);

    std::vector<std::vector<int>> Played_actionsecond(T);


    std::vector<double> original_capacities = getCapacities(network);

    for (int t = 0; t < T; ++t) {
        std::cout << "Ronda " << t << std::endl;
        std::vector<double> Capacities_t = Capacities[Contexts[t]]; // se coge la capacidad en la ronda t en funcion del contexto en esa ronda
        std::vector<int> played_actions_t(N);

        // 1 - Cada jugador juega una acción
        for (int i = 0; i < N; ++i) {
            if (Players[i]->getType() == PlayerType::cGPMW && t > 0) {
                Players[i]->computeStrategys(Capacities_t); 
            }
            played_actions_t[i] = Players[i]->sample_action();  // los jugadores no controlados van a usar siempre su único brazo que es el 0
        }
        
        int asfd = 243;
        int actioncomp = 0;
        if (haycomp) actioncomp = compPlayer->sample_action();
        
        this->Played_actions[t] = played_actions_t; // guarda las acciones de todos los jugadores de la ronda t
       
        if (haycomp) playedactions[t] = actioncomp;

        // 2 - Asignar Recompensas/losses
                int identificador = -1; // se asigna el identificador -1 para que compute travel times calcule los tiempos para todos los jugadores
        std::vector<double> losses_t = Compute_traveltimes(network, Strategy_vectors, this->Played_actions[t], identificador, Capacities_t);
        double lossesrondat = 0;
        if (haycomp) {
            std::vector<int> tmp = Played_actions[t];
            tmp[id] = actioncomp;
             std::vector<double> lossestmp = Compute_traveltimes(network, Strategy_vectors, tmp, id, Capacities_t);
            lossesrondat = lossestmp[id];
            chosen_arms[t] = {actioncomp, Played_actions[t][id]};
            lossesdef[t] = losses_t[id];
            lossescomp[t] = lossesrondat;
            if (t > 0) {
                cumlosscomp[t] = lossescomp[t] + cumlosscomp[t - 1];
                cumlossdef[t] = lossesdef[t] + cumlossdef[t - 1];
            }
            else {
                cumlosscomp[t] = lossescomp[t];
                cumlossdef[t] = lossesdef[t];
            }
        }
        this->Incurred_losses[t] = losses_t ; // Incurred_losses[t][player_id] --> para la ronda t devuelve el ARREPENTIMIENTO del jugador player_id

        // calculamos perdidas acumuladas 
        if (t > 0) {
            std::vector<double> cum_losses_t(N);
            for (int i = 0; i < N; i++) {
                cum_losses_t[i] = this->Cum_losses[t - 1][i] + this->Incurred_losses[t][i];
            }
            if(haycomp) incurredlosses += lossesrondat;
            this->Cum_losses[t] = cum_losses_t;
        }
        else this->Cum_losses[t] = losses_t;

        int E = Strategy_vectors[0][0].size(); // numero de carreteras
        Total_occupancies.push_back(std::vector<double>(E, 0.0)); // --> Total_occupancies[t][e] para cada ronda (t) muestra la ocupación de la carretera (e) que es la suma de la ocupación que realiza cada jugador (i) de esa carretera en el momento t 
        
        /* Cada occupancies (cada ronda) guarda lo que se ocupa cada carretera sumando 
        * lo que usan los jugadores estas carreteras mediante las estrategias o brazos elegidos
        */
        std::vector<double> congestions(E, 0.0); 

        for (int i = 0; i < N; ++i) {
            int aux = Strategy_vectors[i][this->Played_actions[t][i]].size(); // aux = nº de carreteras --> hubiese valido con sustituir this->Played_actions[t][i] por 0 pq size() siempre va a devolver el nº de carreteras (76 para todos) y también se contempla el caso en el que solo juegue una acción
            for (int j = 0; j < aux; ++j) { // 76 carreteras: sumar lo que ocupa en total la estrategia de cada jugador en las 76 carreteras
                Total_occupancies[t][j] += Strategy_vectors[i][this->Played_actions[t][i]][j]; // en cada total occupancies se guarda todo lo que ocupa un jugador en todas sus carreteras del brazo

            }

            for (int e = 0; e < Capacities_t.size(); ++e) {
                congestions[e] = 0.15 * std::pow(Total_occupancies[t][e] / Capacities_t[e], 4);
            }

        }
        addit_Congestions[t] = congestions; // addit_Congestions[t][e] muestra la congestion en el momento t de la carretera e
            
        // 3 - Actualizar estrategias
        for (int i = 0; i < N; ++i) {
            int as = 2;
            if (Players[i]->getType() == PlayerType::Hedge) {
                if(Players[i]->getK() > 1)
                    Players[i]->Update(this->Played_actions[t], i, network, Capacities_t, Strategy_vectors);
                    /*Se pasan las elecciones de todos los jugadores (retroalimentación alta), datos de la red, las capacidades de cada carretera en 
                    el momento t (en cada momento t hay un contexto distinto y por tanto capacidades distintas), y las estrategias de todos los jugadores */
            }
            if (Players[i]->getType() == PlayerType::GPMW) {
                double mean = 0.0;  // media
                double std_dev = sigmas[i];  // desviación estándar
                std::mt19937 gen(1234);  // semilla del generador de números aleatorios
                std::normal_distribution<double> dist(mean, std_dev);  // distribución normal
                double noise = dist(gen);  // generar una muestra aleatoria. dist es un objeto de la clase std::normal_distribution que representa la distribución normal con los parámetros especificados. 
                double noisy_loss = Incurred_losses[t][i] + noise;
                Players[i]->Update(t, this->Played_actions[t][i], Total_occupancies.back(), -noisy_loss, Capacities_t); 
                /*Se pasan las elecciones de todos los jugadores (retroalimentación baja), la ocupación de las carretereas en la última ronda, errepentimiento ruidoso, y las capacidades de cada carretera en 
                    el momento t (en cada momento t hay un contexto distinto y por tanto capacidades distintas)*/
            }

            if (Players[i]->getType() == PlayerType::cGPMW && t > 0) {
                double mean = 0.0;  // media
                double std_dev = sigmas[i];  // desviación estándar
                std::mt19937 gen(1234);  // semilla del generador de números aleatorios
                std::normal_distribution<double> dist(mean, std_dev);  // distribución normal
                double noise = dist(gen);  // generar una muestra aleatoria. dist es un objeto de la clase std::normal_distribution que representa la distribución normal con los parámetros especificados. 
                double noisy_loss = Incurred_losses[t][i] + noise;
                
                Players[i]->UpdateHistory(t, this->Played_actions[t][i], Total_occupancies.back(), -noisy_loss,Capacities_t);
            }

        }
        if (haycomp) {
            double mean = 0.0;  // media
            double std_dev = sigmas[id];  // desviación estándar
            std::mt19937 gen(1234);  // semilla del generador de números aleatorios
            std::normal_distribution<double> dist(mean, std_dev);  // distribución normal
            double noise = dist(gen);  // generar una muestra aleatoria. dist es un objeto de la clase std::normal_distribution que representa la distribución normal con los parámetros especificados. 
            double noisy_loss = lossesrondat + noise;
            compPlayer->Update(t, actioncomp, Total_occupancies.back(), -noisy_loss, Capacities_t);
        }
        


        double avg_cong = 0;
        for (int i = 0; i < addit_Congestions.size(); i++) {
            double sum = 0;
            for (int j = 0; j < addit_Congestions[i].size(); j++) {
                sum += addit_Congestions[i][j];
            }
            avg_cong += sum / addit_Congestions[i].size();
        }

        avg_cong /= addit_Congestions.size();
        }
        return incurredlosses;
}

void Initialize_Players(int N, const std::vector<std::pair<int, int>>& od_Pairs, std::vector<std::vector<std::vector<int>>> &Strategy_vectors, std::vector<double> min_traveltimes, std::vector<double> max_traveltimes, std::vector<int> idxs_controlled, double T, std::string Algo, int version, std::vector<double> Sigma, std::vector<Eigen::MatrixXd>& Kernels, std::vector<double> sigmas, int numberofcontexts, std::vector<std::vector<double>> Capacities, std::vector<Player*>& Players, int &id, Player*& compPlayer) {
    bool one = false;
    for (int i = 0; i < N; i++) {
        int K_i = Strategy_vectors[i].size();
        double min_payoff = -max_traveltimes[i]; // min recompensa = - max tiempo viaje
        double max_payoff = -min_traveltimes[i];
        // cambiar la forma de tratar players
        // idxs son los ids de los agentes que son controlados
        if (find(idxs_controlled.begin(), idxs_controlled.end(), i) != idxs_controlled.end() && K_i > 1) { // si el agente está controlado por el agente y tiene más de un brazo
            if (Algo == "Hedge") {
                Players[i] = new Player_Hedge(K_i, T, min_payoff, max_payoff);
            }
            else if (Algo == "GPMW") {
                Players[i] = new Player_GPMW(K_i, T, min_payoff, max_payoff, Strategy_vectors[i], Kernels[i], sigmas[i]);
            }
            else if (Algo == "cGPMW") {
                Players[i] = new Player_cGPMW(K_i, T, min_payoff, max_payoff, Capacities, Strategy_vectors[i], Kernels[i], sigmas[i]);
                if (!one) {
                    compPlayer = new Player_GPMW(K_i, T, min_payoff, max_payoff, Strategy_vectors[i], Kernels[i], sigmas[i]);
                    id = i;
                    one = true;
                }
            }
        }
        else {
            K_i = 1;
            Players[i] = new Player_Hedge(K_i, T, min_payoff, max_payoff);
            for(int brazoborrado = 0;brazoborrado < 4; brazoborrado++)
                Strategy_vectors[i].pop_back();
        }
    }

}


std::vector<Eigen::MatrixXd> Optimize_Kernels(bool reoptimize, std::string Algo,  const std::vector<int>& idxs_controlled, const std::vector<std::vector<std::vector<int>>>& Strategy_vectors, const std::vector<double>& sigmas, int poly_degree, const std::vector<std::vector<double>>& Outcomes, const std::vector<std::vector<double>>& Capacities, const std::vector<std::vector<double>>& Payoffs, std::vector<std::vector<double>>& list_of_param_arrays)
{
    std::vector<Eigen::MatrixXd> Kernels(Strategy_vectors.size());

    // cargar parametros para el algoritmo Algo
    std::string filename = "list_of_param_arrays_" + Algo + ".txt";
    list_of_param_arrays = loadParamsFromFile(filename);

    // Kernel tiene N jugadores pero solo se usan los indices de los jugadores controlados, el resto de los 500 jugadores estan vacios.
    // Por lo tanto, para acceder a cada kernel hay que usar el indice de los jugadores controlados como la posición del kernel
    for (int jugador = 0; jugador < idxs_controlled.size(); jugador++) { // para cada jugador de los 20
        int ind = idxs_controlled[jugador];
        if (Kernels[ind].isZero()) { 

            std::vector<int> idx_nonzeros; // se usa en reoptimizacion -> entiendo que guarda las estrategias o brazos (de los 5 que hay) si no tienen valores de 0
            for (int carr = 0; carr < Strategy_vectors[ind][0].size(); carr++) { // para cada carretera se suman los valores de los caminos para ver si el jugador pasa por la carretera en algun camino
                int suma = 0;
                for (int camino = 0; camino < Strategy_vectors[ind].size(); ++camino) { 
                    suma += Strategy_vectors[ind][camino][carr];
                }
                if (suma != 0) idx_nonzeros.push_back(suma);
            }
            const int dim = idx_nonzeros.size();

            if (reoptimize == false) { // se hace en el init
                std::vector<double> loaded_params = list_of_param_arrays[ind]; // cargando parametros desde lista para Algo
                Eigen::VectorXd variances(2);
                variances << loaded_params[0], loaded_params[3];
                Eigen::VectorXd scales(2);
                scales << loaded_params[1], loaded_params[4];
                Eigen::VectorXd biases(2);
                biases << loaded_params[2], loaded_params[5];
                Eigen::MatrixXd active_dims(2, dim);
                for (int i = 0; i < dim; i++) {
                    active_dims(0, i) = i;
                    active_dims(1, i) = i + dim;
                }

                Eigen::MatrixXd X = Eigen::MatrixXd::Identity(dim, dim); // se usa la matriz identidad

                Eigen::MatrixXd kernel_1 = Eigen::MatrixXd::Zero(dim, dim);
                // Kernels polinomicos
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j) {
                        double dot_product = X.row(i).dot(X.row(j));
                        kernel_1(i, j) = std::pow(loaded_params[1] * dot_product + loaded_params[2], 1.0);
                    }
                }
                Eigen::MatrixXd kernel_2 = Eigen::MatrixXd::Zero(dim, dim);
                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j) {
                        double dot_product = X.row(i).dot(X.row(j));
                        kernel_2(i, j) = std::pow(variances(1) * dot_product + biases(1), static_cast<double>(poly_degree));
                    }
                }

                Kernels[ind] = kernel_1.cwiseProduct(kernel_2);
                


            }
        }

    }
    return Kernels;
}


std::vector<std::vector<double>> loadParamsFromFile(std::string fileName)
{
    std::ifstream file(fileName);
    std::vector<std::vector<double>> params;

    if (file.is_open()) {
        std::string line;
        while (getline(file, line)) {
            std::istringstream iss(line);
            std::vector<double> paramLine;
            double val;
            while (iss >> val) {
                paramLine.push_back(val);
            }
            params.push_back(paramLine);
        }
        file.close();
    }
    else {
        std::cerr << "No se pudo abrir el archivo: " << fileName << std::endl;
    }

    return params;
}
