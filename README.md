# GeneticAlgorithmInC

This is my implementation of a genetic algorithm in C.
I've created this out of curiosity, as I've implemented several evolutionary algorithms in Python before, and just wanted to compare the performance of C and Python on a evolutionary algorithm. I also wanted to program in C again, as it is one of the coolest programming languages out there.

## The optimization problem

This code currently optimizes the Rosenbrock function with n = 3.
The Rosenbrock function is a tough function to optimize, read more about it on [Wikipedia](https://en.wikipedia.org/wiki/Rosenbrock_function). You could also change the `double calculate_function(Chromosome chromosome)` function to any other [Test function for optimization](https://en.wikipedia.org/wiki/Test_functions_for_optimization).

## Output

The program creates a .txt with average and best fitness values for each generation. See `docs/fitness_values.txt` for an example. To visualize the progress better, a short python script has been used.

![Optimization Plot](docs/optimization_plot.png)