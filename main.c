#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define POPULATION_SIZE 10
#define GENE_COUNT 3

typedef struct chromosome {
    float genes[GENE_COUNT];
    float fitness;
} Chromosome;

float random_float();
int random_int();
float calculate_fitness(Chromosome chromosome);
Chromosome *generate_starting_population();
void printf_chromosome(Chromosome chromosome);

int main(int argc, char *argv[]) {
    // set random seed
    srand((unsigned int)time(NULL));

    Chromosome* startingPopulation = generate_starting_population();

    for(int i=0; i < POPULATION_SIZE; i++){
        printf_chromosome(startingPopulation[i]);
    }

    free(startingPopulation);
    return 0;
}

float random_float(){
    // returns a random float in 0 <= x < 1
    return (float)rand()/(float)(RAND_MAX);
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
            chromosomes[i].genes[g] = random_float();
        }
        chromosomes[i].fitness = calculate_fitness(chromosomes[i]);
    }
    return chromosomes;
}

void printf_chromosome(Chromosome chromosome){
    printf("%f : %f, %f, %f\n", chromosome.fitness, chromosome.genes[0], chromosome.genes[1], chromosome.genes[2]);
}