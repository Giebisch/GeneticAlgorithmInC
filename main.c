#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define GENE_COUNT          3
#define GENERATIONS         100
#define POPULATION_SIZE     10
#define TOURNAMENT_SIZE     2
#define MUTATION_CHANCE     0.3
#define CROSSOVER_CHANCE    0.7
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
Chromosome *uniform_crossover(Chromosome *population);
Chromosome *gaussian_mutation(Chromosome *population);

int main(int argc, char *argv[]) {
    // set random seed
    srand((unsigned int)time(NULL));

    Chromosome* population = generate_starting_population();

    for(int i=0; i < GENERATIONS; i++){
        population = tournament_selection(population);
        population = uniform_crossover(population);
        population = gaussian_mutation(population);
        Chromosome best = get_best_chromosome(population, POPULATION_SIZE);
        printf_chromosome(best);
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
    // g1 ** 2 + g2 ** 2 + g3 ** 2
    double value = pow(chromosome.genes[0], 2) + pow(chromosome.genes[1], 2) + pow(chromosome.genes[2], 2);
    return value;
}

double calculate_fitness(Chromosome chromosome){
    double value = calculate_function(chromosome);
    return 1./(1. + fabs(value));
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
    printf("%.20lf : %lf, %lf, %lf\n", chromosome.fitness, chromosome.genes[0], chromosome.genes[1], chromosome.genes[2]);
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

Chromosome *uniform_crossover(Chromosome *population){
    Chromosome *new_population = (Chromosome*)calloc(POPULATION_SIZE, sizeof(Chromosome));
    for(int i=0; i < POPULATION_SIZE/2; i++){
        Chromosome rand_chr1 = population[random_int(0, POPULATION_SIZE)];
        Chromosome rand_chr2 = population[random_int(0, POPULATION_SIZE)];
        if(random_double(0, 1) < CROSSOVER_CHANCE){
            int genes1[GENE_COUNT];
            int genes2[GENE_COUNT];
            for(int g=0; g < GENE_COUNT; g++){
                if(random_double(0,1) < 0.5){
                    genes1[g] = rand_chr1.genes[g];
                } else {
                    genes1[g] = rand_chr2.genes[g];
                }
                if(random_double(0,1) < 0.5){
                    genes2[g] = rand_chr1.genes[g];
                } else {
                    genes2[g] = rand_chr2.genes[g];
                }
            }
            Chromosome new_chr1, new_chr2;
            for(int gene = 0; gene < GENE_COUNT; gene++){
                new_chr1.genes[gene] = genes1[gene];
            }
            for(int gene = 0; gene < GENE_COUNT; gene++){
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