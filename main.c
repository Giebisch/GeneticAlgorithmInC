#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define GENE_COUNT          3
#define GENERATIONS         1000
#define POPULATION_SIZE     100
#define TOURNAMENT_SIZE     10
#define MUTATION_CHANCE     0.3
#define LOW                 -5.12
#define HIGH                5.12

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

Chromosome get_best_chromosome(Chromosome *population, int size);

// genetic algorithm functions
Chromosome *generate_starting_population();
Chromosome *tournament_selection(Chromosome *population);
Chromosome *gaussian_mutation(Chromosome *population);

int main(int argc, char *argv[]) {
    // set random seed
    srand((unsigned int)time(NULL));

    Chromosome* population = generate_starting_population();

    for(int i=0; i < GENERATIONS; i++){
        population = gaussian_mutation(population);
        population = tournament_selection(population);
    }

    Chromosome best = get_best_chromosome(population, POPULATION_SIZE);
    printf_chromosome(best);
    printf("%lf", calculate_function(best));
    free(population);
    return 0;
}

double random_double(double low, double high){
    // returns a random double in 0 <= x < 1
    return ((double)rand()/(double)(RAND_MAX)) * (high - low) + low;
}

int random_int(int low, int high){
    // returns a random int in low <= x < high
    int diff = high - low;
    return (rand() % diff) + low;
}

double calculate_function(Chromosome chromosome){
    double value = pow(chromosome.genes[0], 2) + pow(chromosome.genes[1], 2) + pow(chromosome.genes[2], 2);
    return value;
}

double calculate_fitness(Chromosome chromosome){
    double value = calculate_function(chromosome);
    return 1/(1 + abs(value));
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

void printf_chromosome(Chromosome chromosome){
    printf("%lf : %lf, %lf, %lf\n", chromosome.fitness, chromosome.genes[0], chromosome.genes[1], chromosome.genes[2]);
}

Chromosome *tournament_selection(Chromosome *population){
    Chromosome *new_population = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE; i++){
        // select TOURNAMENT_SIZE-many random chromosomes
        Chromosome *tournament_chromosomes = (Chromosome*)calloc(TOURNAMENT_SIZE, sizeof(Chromosome));
        for(int t=0; t < TOURNAMENT_SIZE; t++){
            int index = random_int(0, POPULATION_SIZE);
            tournament_chromosomes[t] = population[index];
        }

        // find chromosome with the highest fitness value in tournament and add to new population
        new_population[i] = get_best_chromosome(tournament_chromosomes, TOURNAMENT_SIZE);

        free(tournament_chromosomes);
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
            if(random_double(0, 1) < MUTATION_CHANCE){
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

double normal_distribution(){
    // Box-Muller transform
    double y1 = random_double(0, 1);
    double y2 = random_double(0, 1);
    return cos(2*3.14*y2)*sqrt(-2.*log(y1));
}