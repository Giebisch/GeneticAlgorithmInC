#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define GENE_COUNT          3
#define GENERATIONS         750
#define POPULATION_SIZE     100
#define TOURNAMENT_SIZE     10
#define MUTATION_CHANCE     0.3
#define CROSSOVER_CHANCE    0.7
#define LOW                 -5.
#define HIGH                5.

typedef struct chromosome {
    double genes[GENE_COUNT];
    double fitness;
} Chromosome;

// helper functions
int random_int();
double normal_distribution();
double random_double(double low, double high);
double calculate_fitness(Chromosome chromosome);
double calculate_function(Chromosome chromosome);
void printf_chromosome(Chromosome chromosome);
double get_average_fitness(Chromosome *population);
Chromosome get_best_chromosome(Chromosome *population, int size);

// genetic algorithm functions
Chromosome *generate_starting_population();
Chromosome *tournament_selection(Chromosome *population);
Chromosome *uniform_crossover(Chromosome *population);
Chromosome *gaussian_mutation(Chromosome *population);

int main(int argc, char *argv[]) {
    // set random seed
    srand((unsigned int)time(NULL));
    FILE *out_file = fopen("fitness_values.txt", "w");
    fprintf(out_file, "AVERAGE\tBEST\n");

    Chromosome* population = generate_starting_population();

    for(int i=0; i < GENERATIONS; i++){
        Chromosome* selected_chromosomes = tournament_selection(population);
        Chromosome* crossover_chromosomes = uniform_crossover(selected_chromosomes);
        Chromosome* mutated_chromosomes = gaussian_mutation(crossover_chromosomes);
        
        Chromosome best = get_best_chromosome(mutated_chromosomes, POPULATION_SIZE);
        // printf_chromosome(best);
        
        for(int t = 0; t < POPULATION_SIZE; t++){
            population[t] = mutated_chromosomes[t];
        }
        fprintf(out_file, "%lf\t%lf\n", get_average_fitness(population), best.fitness);
    }

    Chromosome best = get_best_chromosome(population, POPULATION_SIZE);
    printf_chromosome(best);
    printf("%lf", calculate_function(best));
    free(population);

    return 0;
}

double random_double(double low, double high){
    // returns a random double in low <= x < high
    return ((double)rand()/(double)(RAND_MAX)) * (high - low) + low;
}

int random_int(int low, int high){
    // returns a random int in low <= x < high
    int diff = high - low;
    return (rand() % diff) + low;
}

double calculate_function(Chromosome chromosome){
    // Rosenbrock function with n = 3
    double value = 100. * pow((chromosome.genes[1] - pow(chromosome.genes[0], 2.)),2.) + pow((1. - chromosome.genes[0]), 2.)
                 + 100. * pow((chromosome.genes[2] - pow(chromosome.genes[1], 2.)),2.) + pow((1. - chromosome.genes[1]), 2.);
    return value;
}

double calculate_fitness(Chromosome chromosome){
    double value = calculate_function(chromosome);
    return 1./(1. + fabs(value));
}

double get_average_fitness(Chromosome *population){
    double average = 0;
    for(int i=0; i < POPULATION_SIZE; i++){
        average = average + population[i].fitness;
    }
    return average/POPULATION_SIZE;
}

void printf_chromosome(Chromosome chromosome){
    printf("%lf : %lf, %lf, %lf\n", chromosome.fitness, chromosome.genes[0], chromosome.genes[1], chromosome.genes[2]);
}

double normal_distribution(){
    // Box-Muller transform
    double y1 = random_double(0., 1.);
    double y2 = random_double(0., 1.);
    return cos(2.*3.14*y2)*sqrt(-2.*log(y1));
}

Chromosome *generate_starting_population(){
    Chromosome *chromosomes = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE; i++){
        for(int g=0; g < GENE_COUNT; g++){
            chromosomes[i].genes[g] = random_double(LOW, HIGH);
        }
        chromosomes[i].fitness = calculate_fitness(chromosomes[i]);
    }
    return chromosomes;
}

Chromosome *tournament_selection(Chromosome *population){
    Chromosome *new_population = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE; i++){
        // select TOURNAMENT_SIZE-many random chromosomes
        Chromosome tournament_chromosomes[TOURNAMENT_SIZE];
        for(int t=0; t < TOURNAMENT_SIZE; t++){
            int index = random_int(0, POPULATION_SIZE);
            tournament_chromosomes[t] = population[index];
        }
        // find chromosome with the highest fitness value in tournament and add to new population
        new_population[i] = get_best_chromosome(tournament_chromosomes, TOURNAMENT_SIZE);
    }
    return new_population;
}

Chromosome get_best_chromosome(Chromosome *population, int size){
    double max_fitness = 0.0;
    int winner_index = 0;
    for(int c=0; c < size; c++){
        if(population[c].fitness > max_fitness){
            max_fitness = population[c].fitness;
            winner_index = c;
        }
    }
    return population[winner_index];
}

Chromosome *gaussian_mutation(Chromosome *population){
    Chromosome *new_population = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE; i++){
        double new_genes[GENE_COUNT];
        for(int g=0; g < GENE_COUNT; g++){
            double old_gene = population[i].genes[g];
            if(random_double(0., 1.) < MUTATION_CHANCE){
                new_genes[g] = old_gene + normal_distribution() * 0.1 * (HIGH - LOW);
            } else {
                new_genes[g] = old_gene;
            }
        }
        Chromosome new_chromosome;
        for(int gene = 0; gene < GENE_COUNT; gene++){
            new_chromosome.genes[gene] = new_genes[gene];
        }
        new_chromosome.fitness = calculate_fitness(new_chromosome);
        new_population[i] = new_chromosome;
    }
    return new_population;
}

Chromosome *uniform_crossover(Chromosome *population){
    Chromosome *new_population = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE/2; i++){
        Chromosome rand_chr1 = population[random_int(0, POPULATION_SIZE)];
        Chromosome rand_chr2 = population[random_int(0, POPULATION_SIZE)];
        if(random_double(0., 1.) < CROSSOVER_CHANCE){
            double genes1[GENE_COUNT];
            double genes2[GENE_COUNT];
            for(int g=0; g < GENE_COUNT; g++){
                if(random_double(0.,1.) < 0.5){
                    genes1[g] = rand_chr1.genes[g];
                } else {
                    genes1[g] = rand_chr2.genes[g];
                }
                if(random_double(0.,1.) < 0.5){
                    genes2[g] = rand_chr1.genes[g];
                } else {
                    genes2[g] = rand_chr2.genes[g];
                }
            }
            Chromosome new_chr1, new_chr2;
            for(int gene = 0; gene < GENE_COUNT; gene++){
                new_chr1.genes[gene] = genes1[gene];
                new_chr2.genes[gene] = genes2[gene];
            }
            new_chr1.fitness = calculate_fitness(new_chr1);
            new_chr2.fitness = calculate_fitness(new_chr2);
            new_population[i] = new_chr1;
            new_population[POPULATION_SIZE - 1 - i] = new_chr2;
        } else {
            new_population[i] = rand_chr1;
            new_population[POPULATION_SIZE - 1 - i] = rand_chr2;
        }
    }
    return new_population;
}