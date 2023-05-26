# computer-science-final-degree-project-2023
Repository for the source code of my Computer Science Final Degree Project. The project investigates the contextual multi-armed bandits.
This work involves the development of several implementations (each of these is referenced in the order they appear in the main body of the thesis):

Sections 6.1 and 6.2 pertain to the initial implementations of contextual multi-armed bandits.
Section 6.1 discusses a system with few contexts using instances of the Upper Confidence Bounds (UCB) algorithm, a stochastic bandit method, as an introductory contextual bandit.
Section 6.2 delves into Lipschitzian contextual bandits which employ a discretization process and utilize instances of UCB.

Section 7 involves a contextual game, encompassing the development of the Hedge, GPMW, and cGPMW algorithms. In a network setting, the algorithms are contrasted, studying how cGPMW employs contextual information. A contextual game is applied, with graph-related algorithms for short path computation and GPMW and cGPMW algorithms utilizing an RVM regression model.

Section 8 presents a trading bot system utilizing EXP4. Multiple financial experts are tested using advanced concepts from the capital market to predict market trends. An analysis of exploration and learning rates is conducted.

Section 9 details a movie recommendation system that leverages EXP4. It employs various policies such as Epsilon-Greedy, a neural network adapted as a contextual bandit, and collaborative filtering using a neighborhood algorithm to construct a high-quality recommendation system.

Authors: Iván Hernández Roldán and Alejandro Magarzo Gonzalo.
