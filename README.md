# computer-science-final-degree-project-2023
Source code of my Computer Science Bachelor's Thesis (Contexutal_Bandits_Fundamentals_and_Applications.pdf). This project develops a study of the contextual multi-armed bandits, a type of problem in the scope of reinforcement learning. 

This repository includes the implementation of some practical applications of the contextual multi-armed bandits' algorithms (each of them are referenced in the same order as they appear in the main body of the thesis):

- **Section 6.1** develops the first recommendation system in a context with few types of users using instances of the Upper Confidence Bounds (UCB) algorithm, a stochastic bandit method, as an introductory contextual bandit.
- **Section 6.2** delves into Lipschitzian contextual bandits in order to employ a discretization process to create a recommendation system in a context with several types of user. After the discretization process, the situation is similiar to the one in Section 6.1 and it is possible to use instances of UCB.

- **Section 7** involves a contextual routing game based on a network of roads and nodes, encompassing the development of the Hedge, GPMW, and cGPMW algorithms. In this network, the algorithms are contrasted, studying how the cGPMW algorithm, which uses an RVM regression model, takes advantage of the available contextual information to make better decisions about the shortest path for a specific journey.

- **Section 8** presents a trading bot system using EXP4 algorithm. Multiple financial experts, who use advanced concepts from the capital market to predict market trends, are tested in order to select depending on the contest the best one and follow its investment advice. Also, an analysis of exploration and learning rates is conducted.

- **Section 9** details a movie recommendation system that also leverages EXP4 algorithm. It employs various policies such as Epsilon-Greedy, a neural network adapted as a contextual bandit, and collaborative filtering using a neighborhood algorithm to construct a high-quality recommendation system.

Authors: Iván Hernández Roldán and Alejandro Magarzo Gonzalo.
