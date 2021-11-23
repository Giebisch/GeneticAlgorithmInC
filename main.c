#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define GENE_COUNT          3
#define GENERATIONS         500
#define POPULATION_SIZE     500
#define TOURNAMENT_SIZE     2

typedef struct chromosome {
    float genes[GENE_COUNT];
    float fitness;
} Chromosome;

// helper functions
int random_int();
float random_float(float low, float high);
float calculate_fitness(Chromosome chromosome);
void printf_chromosome(Chromosome chromosome);
Chromosome get_best_chromosome(Chromosome *population, int size);

// genetic algorithm functions
Chromosome *generate_starting_population();
Chromosome *tournament_selection(Chromosome *population);

int main(int argc, char *argv[]) {
    // set random seed
    srand((unsigned int)time(NULL));

    Chromosome* startingPopulation = generate_starting_population();

    for(int i=0; i < GENERATIONS; i++){
        Chromosome best = get_best_chromosome(startingPopulation, POPULATION_SIZE);
        printf_chromosome(best);
        startingPopulation = tournament_selection(startingPopulation);
    }

    free(startingPopulation);
    return 0;
}

float random_float(float low, float high){
    // returns a random float in 0 <= x < 1
    return ((float)rand()/(float)(RAND_MAX)) * (high - low) + low;
}

int random_int(int low, int high){
    // returns a random int in low <= x < high
    int diff = high - low;
    return (rand() % diff) + low;
}

float calculate_fitness(Chromosome chromosome){
    // dummy function for now
    float value = chromosome.genes[0] * chromosome.genes[0] + chromosome.genes[1] * chromosome.genes[1] + chromosome.genes[2] * chromosome.genes[2];
    return 1/(1 + abs(value));
}

Chromosome *generate_starting_population(){
    Chromosome *chromosomes = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE; i++){
        for(int g=0; g < GENE_COUNT; g++){
            chromosomes[i].genes[g] = random_float(-10, 10);
        }
        chromosomes[i].fitness = calculate_fitness(chromosomes[i]);
    }
    return chromosomes;
}

void printf_chromosome(Chromosome chromosome){
    printf("%f : %f, %f, %f\n", chromosome.fitness, chromosome.genes[0], chromosome.genes[1], chromosome.genes[2]);
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
    int max_fitness = 0;
    int winner_index = 0;
    for(int c=0; c < size; c++){
        if(population[c].fitness > max_fitness){
            max_fitness = population[c].fitness;
            winner_index = c;
        }
    }
    return population[winner_index];
}
