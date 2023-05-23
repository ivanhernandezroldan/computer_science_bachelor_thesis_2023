#pragma once
#include "simulation.h"
#include "network.h"

int main(...)
{
	NetworkData network = createNetwork();
	Simulation simulation(network);
	simulation.init();
	return 0;
}