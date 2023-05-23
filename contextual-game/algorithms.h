#ifndef alogrithms_h
#define algorithms_h


#include <vector>
#include <cmath>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <string>
#include <numeric>
#include <iostream>
#include <random>

#include <dlib/matrix.h>
#include <dlib/statistics.h>
#include <dlib/svm_threaded.h>
#include <dlib/svm.h>
#include <dlib/optimization.h>

#include "network.h"

typedef dlib::matrix<double, 0, 1> sample_type;
typedef dlib::radial_basis_kernel<sample_type> kernel_type;
void print_dlib_X_train(const std::vector<sample_type>& dlib_X_train, int brazo,  std::vector<double> payoffs);
std::vector<sample_type> history_to_dlib_samples(const std::vector<std::vector<double>>& history);
std::vector<double> history_payoffs_to_dlib_labels(const std::vector<double>& history_payoffs);
double calculate_residual_variance(dlib::decision_function<kernel_type>& model, const std::vector<sample_type>& dlib_X_train, const std::vector<double>& dlib_y_train);
std::vector<std::vector<double>> normalize_vector_of_vectors(const std::vector<std::vector<double>>& vec_of_vecs, double& min_value, double& max_value);
enum class PlayerType {
    cGPMW,
    Hedge,
    GPMW
};

class Player {
protected:
    PlayerType type_;
    int K_;
    double min_payoff_;
    double max_payoff_;
    std::vector<double> weights_;
    double T_;
    double gamma_t_;

public:

    virtual int sample_action();
    //Hedge
    virtual void Update(std::vector<int> played_actions, int player_idx, const NetworkData& network, std::vector<double> Capacities_t, std::vector<std::vector<std::vector<int>>> Strategy_vectors);
    //GPMW
    virtual void Update(int ronda, int played_action, std::vector<double> total_occupancies, double payoff, std::vector<double> Capacities_t);
    //cGPMW
    virtual void UpdateHistory(int ronda, int played_action, std::vector<double> total_occupancies, double payoff, std::vector<double> capacitites);
    virtual void computeStrategys(const std::vector<double>& capacities_t);
    int getK();
    PlayerType getType();

};


class Player_Hedge : public Player {
public:
    Player_Hedge(int K, double T, double min_payoff, double max_payoff) {
        this->K_ = K;
        this->T_ = T;
        this->min_payoff_ = min_payoff;
        this->max_payoff_ = max_payoff;
        this->gamma_t_ = (sqrt(8 * log(K) / T));// tasa aprendizaje
        this->type_ = PlayerType::Hedge;
        this->weights_ = std::vector<double>(K, 1); // para cada brazo el valor inicial en la distribución es 1
    }

    std::vector<double> mixed_strategy();
    int sample_action() override;
    void Update(std::vector<int> played_actions, int player_idx, const NetworkData& network, std::vector<double> Capacities_t, std::vector<std::vector<std::vector<int>>> Strategy_vectors) override;
    
    // funciones auxiliares 
    int get_K() const { return K_; }
    double get_T() const { return T_; }
    double get_min_payoff() const { return min_payoff_; }
    double get_max_payoff() const { return max_payoff_; }
};



class Player_GPMW : public Player {
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXr;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;


public:
    Player_GPMW(int K, int T, double min_payoff, double max_payoff, std::vector<std::vector<int>>  my_strategy_vector, Eigen::MatrixXd kernel, double sigma_e) {
        this->type_ = PlayerType::GPMW;
        this->K_ = K;
        this->T_ = T;
        this->min_payoff_ = min_payoff;
        this->max_payoff_ = max_payoff;
        this->weights_ = std::vector<double>(K, 1);
        // para cada carretera se suman los valores de los caminos para ver si el jugador pasa por esa carretera en algun camino
        this->idx_nonzeros;
        for (int carr = 0; carr < my_strategy_vector[0].size(); carr++) { //Strategy_vectors[ind][0].size() da igual si cogiese Strategy_vectors[ind][1].size() pq todos tienen size = 76 (carreteras)
            int suma = 0;
            for (int camino = 0; camino < my_strategy_vector.size(); ++camino) {
                suma += my_strategy_vector[camino][carr];
            }
            if (suma != 0) this->idx_nonzeros.push_back(carr); 
        }

        this->cum_losses = std::vector<double>(K, 0.0);
        this->ucb_rewards_est = std::vector<double>(K, 0.0);
        this->gamma_t_ = std::sqrt(8 * log(K) / T); // tasa aprendizaje
        this->kernel = kernel;
        this->sigma_e = sigma_e;
        this->strategy_vecs = my_strategy_vector;

        demand = *std::max_element(my_strategy_vector[0].begin(), my_strategy_vector[0].end());
    }

    std::vector<double> mixed_strategy();
    int sample_action() override;
    void Update(int ronda, int played_action, std::vector<double> total_occupancies, double payoff, std::vector<double> Capacities_t) override;


private:
    std::vector<int> idx_nonzeros;
    std::vector<int> played_actions;

    std::vector<double> cum_losses;
    std::vector<double> mean_rewards_est;
    std::vector<double> std_rewards_est;
    std::vector<double> ucb_rewards_est;
    Eigen::MatrixXd kernel;
    double sigma_e;
    std::vector<std::vector<int>> strategy_vecs;

    std::vector<double> history_payoffs;
    std::vector<std::vector<double>> history;
    double demand;
};




class Player_cGPMW : public Player {
public:
    Player_cGPMW(int K, int T, double min_payoff, double max_payoff, std::vector<std::vector<double>> Capacities, std::vector<std::vector<int>>  my_strategy_vector, Eigen::MatrixXd kernel, double sigma_e) {
        this->type_ = PlayerType::cGPMW;
        this->K_ = K;
        this->T_ = T;
        this->min_payoff_ = min_payoff;
        this->max_payoff_ = max_payoff;
        this->weights_ = std::vector<double>(K, 1);
        // para cada carretera se suman los valores de los caminos para ver si el jugador pasa por esa carretera en algun camino
        this->idx_nonzeros;
        for (int carr = 0; carr < my_strategy_vector[0].size(); carr++) { // Strategy_vectors[ind][0].size() da igual si cogiese Strategy_vectors[ind][1].size() pq todos tienen size = 76 (carreteras)
            int suma = 0;
            for (int camino = 0; camino < my_strategy_vector.size(); ++camino) {
                suma += my_strategy_vector[camino][carr];
            }
            if (suma != 0) this->idx_nonzeros.push_back(carr); 
        }

        this->gamma_t_ = std::sqrt(8 * log(K) / T); // tasa aprendizaje
        this->kernel = kernel;
        this->sigma_e = sigma_e;
        this->strategy_vecs = my_strategy_vector;
        this->Capacities = Capacities;
    }

    std::vector<double> mixed_strategy();
    int sample_action() override;
    void UpdateHistory(int ronda, int played_action, std::vector<double> total_occupancies, double payoff, std::vector<double> capacitites) override;
    void computeStrategys(const std::vector<double>& capacities_t) override;

private:
    std::vector<int> idx_nonzeros;

    Eigen::MatrixXd kernel;
    double sigma_e;
    std::vector<std::vector<int>> strategy_vecs;

    std::vector<double> history_payoffs;
    std::vector<std::vector<double>> history;
    std::vector<int> played_actions;


    std::vector<std::vector<double>> Capacities;
    std::vector<std::vector<double>> history_occupancies;
    std::vector<double> contexts;
    std::vector<double> idx_balls;

};

#endif
