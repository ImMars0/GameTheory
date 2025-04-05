# GameTheory

This C++ project implements a solution for a two-player matrix game using the Simplex method for determining optimal mixed strategies and the game value. The code reads a payoff matrix, determines whether a pure strategy equilibrium exists, and if not, transforms the problem into a linear programming problem for both players. Then, it solves the problem using the Simplex algorithm.

Features:

    Payoff Matrix Input: Allows the user to input the payoff matrix for the two players.

    Pure Strategy Check: Checks whether the game has a pure strategy equilibrium by calculating the lower and upper bounds of the game value.

    Linear Programming Transformation: Transforms the matrix game into a linear programming problem for the players, with the problem stated in standard form for the Simplex algorithm.

    Simplex Method: Implements the Simplex method to solve the linear programming problem step-by-step, with detailed tableau outputs at each iteration.

    Optimal Strategy Calculation: Once the solution is found, it calculates the optimal mixed strategies for both players and outputs the game value.

How it works:

    Input Payoff Matrix: The user is prompted to enter the payoff matrix Q for the game.

    Equilibrium Analysis: The program calculates whether a pure strategy equilibrium exists based on the payoff matrix.

    Linear Programming Formulation: If no equilibrium is found, it formulates the problem into linear programming (maximization for Player B, minimization for Player A).

    Simplex Algorithm: The Simplex method is applied to find the optimal mixed strategies and the game value.

    Results: The optimal strategies for both players and the value of the game are displayed.

Example Output:

    Optimal Strategy for Player A: The mixed strategy for Player A, showing the probability of choosing each pure strategy.

    Optimal Strategy for Player B: The mixed strategy for Player B, showing the probability of choosing each pure strategy.

    Game Value: The optimal value of the game.

This project demonstrates how to apply the Simplex method to solve matrix games, and can be a valuable tool for those studying game theory, linear programming, and optimization.
