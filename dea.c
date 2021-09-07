// Fork of https://gitlab.com/simoesusp/disciplinas/-/tree/master/SSC0713-Sistemas-Evolutivos-Aplicados-a-Robotica/AG-Grafico
//
// Author: Alice
// Author: Eduardo do Valle Simoes @simoesusp
// Author: Breno Cunha Queiroz @Brenocq
// Author: Cho young Lim @alperceire

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <errno.h>
#include <limits.h>
#include <float.h>

#define maxx 1000 // inner genes range [0, 1000]
#define pop 10 // inner EA population
#define pop2 5 // outer EA population
#define inner_loop 1000 // inner EA max generation
#define heritage_max 10 // max hereditary store array
// mutation rate of mutation rate (outer)
#define mutation_rate_squared 10 // per 1000  == 1% of range
#define predation_mutation_rate 10 // per 1000 == 1% of range

// ------------- variable for inner ---------------- //
int gen = 1;
// matrix of inner chromosome with 5 genes
double inner_chromosome[pop][5];
double temp[pop][5];
// array of inner fitness result
double fit[pop];
double bestInd[5];
double maxfit = 0.0;
double mean_fit = 0.0;
int maxi = 0;
unsigned char flag_fine_tune = 0;
double ultimate_inner_fit = 0.0;
double ultimate_inner_chromosome[5];
double fit_steps[inner_loop / 5] = {0};

float mutation_rate = 2;
bool bool_tournament = false;
bool delta_cross_over = true;
int period_predation = 10;
int period_synthesis = 10;

// ------------- variable for outer ---------------- //
int gen2 = 1;
float gene_mutation_rate[pop2] = {0}; // gene0
bool gene_selection[pop2] = {true}; // gene1
bool gene_crossover[pop2] = {true}; // gene2
int gene_predation[pop2] = {0}; // gene3
int gene_syntesis[pop2] = {0}; // gene4
double fit2[pop2];
double maxfit2 = 0.0;
double mean_fit2 = 0.0;
int maxi2 = 0;
// matrix of hereditary scores
double heritage[pop2][heritage_max] = {1};
double RMS_heritage[pop2];

// ------------- functions ---------------- //

double fitness(double x[5]);
void print_result(void);
void evaluate_fitness(void);
void populate(void);
void predation(void);
double cross_with_delta(double chromox, double chromoy);
void elitism(void);
void elitism2(void);
void tournament(void);
void pseudo_synthesize(void);
void repopulate(void);
void evalutate_mutation(void);
void iterate_evolutionary_algorithm(void);
void populate2(void);
void evaluate_fitness2(void);
void loop_outer_algorithm(int outer_loop);
void print_result(void);

// ------------- main ---------------- //

int main(int argc, char *argv[]) {
    srand(time(NULL));

    char *p;
    int num;
    errno = 0;
    long conv = strtol(argv[1], &p, 10);
    if (errno != 0 || *p != '\0' || conv > INT_MAX || conv < INT_MIN) {
        printf("need arg: integer: number of generation\n");
        return(1);
    } else {
        num = conv;
    }
    loop_outer_algorithm(num);
    print_result();
    return(0);
}

// int outer_loop: maximum generation of outer EA
void loop_outer_algorithm(int outer_loop) {
    gen2 = 1;
    populate2();
    evaluate_fitness2();
    // print the outer EA generation status
    // printf("generation, best_RMS_fitness, mean_fitness\n");
    // printf("%d, %f, %f\n", gen2, maxfit2, mean_fit2);
    for (int i = 1; i < outer_loop; i++) {
        elitism2();
        evaluate_fitness2();
        gen2++;
        // printf("%d, %f, %f\n", gen2, maxfit2, mean_fit2);
    }
}

// initiate outer chromosomes
void populate2(void) {
    for (int i = 0; i < pop2; i++) {
        gene_mutation_rate[i] = (double) (rand() % 100);
        gene_selection[i] = (bool) (rand() % 2);
        gene_crossover[i] = (bool) (rand() % 2);
        gene_predation[i] = (int) (rand() % (inner_loop - 1)) + 1;
        gene_syntesis[i] = (int) (rand() % (inner_loop - 1)) + 1;
    }
}

// evalutate fitness of outer chromosomes by running inner EA
void evaluate_fitness2(void) {
    double sum_fit = 0.0;
    double sum_heritage = 0.0;
    for (int i = 0; i < pop2; i++) {
        gen = 1;
        mutation_rate = (double) gene_mutation_rate[i];
        bool_tournament = (bool) gene_selection[i];
        delta_cross_over = (bool) gene_crossover[i];
        period_predation = (int) gene_predation[i];
        period_synthesis = (int) gene_syntesis[i];
        populate();
        evaluate_fitness();
        // print the inner EA generation status
        //printf("%d, %f, %f, %f\n", gen, maxfit, mean_fit, mutation_rate);
        for (int j = 1; j < inner_loop; j++) {
            iterate_evolutionary_algorithm();
            if ((j % (inner_loop / 5)) == 0) {
                fit_steps[j/(inner_loop / 5)] = maxfit;
            }
        }
        for (int j = 1; j < (inner_loop / 5); j++) {
            fit2[i] += (fit_steps[j] * 0.1);
        }
        fit2[i] += maxfit;
        if (ultimate_inner_fit < maxfit) {
            ultimate_inner_fit = maxfit;
            for (int inj = 0; inj < 5; inj++) {
                ultimate_inner_chromosome[inj] = inner_chromosome[maxi][inj];
            }
        }
        sum_heritage = fit2[i] * fit2[i];
        for (int j = 1; j < heritage_max; j++) {
            heritage[i][j - 1] = heritage[i][j];
            sum_heritage += (heritage[i][j - 1] * heritage[i][j - 1]);
        }
        heritage[i][heritage_max - 1] = fit2[i];
        sum_fit += fit2[i];
        RMS_heritage[i] = sqrt(sum_heritage / heritage_max);
    }
    mean_fit2 = sum_fit / pop2;

    maxfit2 = RMS_heritage[0];
    maxi2 = 0;
    for (int k = 1; k < pop2; k++) {
        if (RMS_heritage[k] > maxfit2) {
            maxfit2 = RMS_heritage[k];
            maxi2 = k;
        }
    }
}

// populate inner chromosomes
void populate(void) {
    for (int i = 0; i < pop; i++) {
        for (int j = 0; j < 5; j++) {
            inner_chromosome[i][j] = (double) (rand() % maxx);
        }
    }
}

// inner fitness function
double fitness_alt(double x[5]) {
    double y = 0.1*x[0]*x[1] + 0.1*x[2]*x[3] - 0.1*x[4]*x[0] -0.1*x[2]*x[1];
    if (y > FLT_MAX) {
        y == FLT_MAX;
    }
    return(y);
}
// evalute inner chromosomes fitness
void evaluate_fitness(void) {
    double sum_fit = 0.0;

    for (int i = 0; i < pop; i++) {
        fit[i] = fitness(inner_chromosome[i]);
        sum_fit += fit[i];
    }
    mean_fit = sum_fit / pop;

    maxfit = fit[0];
    maxi = 0;
    for (int i = 1; i < pop; i++) {
        if (fit[i] > maxfit) {
            maxfit = fit[i];
            maxi = i;
        }
    }
    for (int i = 1; i < 5; i++) {
        bestInd[i - 1] = bestInd[i];
    }
    bestInd[4] = maxfit;
}

// inner fitness function
double fitness(double x[5]) {
    double y = ((3*cos(0.49*x[0]) + 10*sin(0.25*x[0]) + 0.5*cos(0.01*x[0]) +
               1*sin(0.07*x[1]) + 5*sin(0.25*x[1]) + 4*cos(0.2*x[1]) +
               1*sin(0.12*x[2]) + 3*sin(0.4*x[2]) + 3*cos(0.4*x[2]) +
               5*sin(0.4*x[3]) + 5*sin(0.1*x[3]) + 2*cos(0.5*x[3]) +
               3*sin(0.5*x[4]) + 5*sin(0.8*x[4]) + 1*cos(0.1*x[4])) /
               (2*cos(0.04*x[0]) + 5*sin(0.03*x[3]) + 10*cos(0.05*x[4])));
    if (y > FLT_MAX) {
        y == FLT_MAX;
    }
    return(y);
}


// select crossover and mutate
void iterate_evolutionary_algorithm(void) {
    if (bool_tournament) {
        tournament();
    } else {
        elitism();
    }
    evaluate_fitness();
    evalutate_mutation();
    predation();
    pseudo_synthesize();
    // printf("%d, %f, %f, %f\n", gen, maxfit, mean_fit, mutation_rate);
    gen++;
}

void tournament(void) {
    int a, b, c, parent1, parent2;

    for (int i = 0; i < pop; i++) {
        for (int j = 0; j < 5; j++) {
            temp[i][j] = inner_chromosome[i][j];  // Backup
        }
    }
    // tournament
    for (int i = 0; i < pop; i++) {
        if (i == maxi)    // protect the best chromosome
            continue;

        // select first parent
        a = (rand() % pop);
        b = (rand() % pop);
        if (fit[a] > fit[b])
            parent1 = a;
        else
            parent1 = b;

        // select second parent
        a = (rand() % pop);
        b = (rand() % pop);
        if (fit[a] > fit[b])
            parent2 = a;
        else
            parent2 = b;

        // Cross_over
        for (int j = 0; j < 5; j++) {
            if (delta_cross_over) {
                inner_chromosome[i][j] = cross_with_delta(temp[parent1][j], temp[parent2][j]);
            } else {
                inner_chromosome[i][j] = (inner_chromosome[parent1][j] + inner_chromosome[parent2][j]) / 2.0;
            }
        }
        // mutation
        // 50% to mutate 1 gene
        // 30% to mutate 2 gene
        // 10% to mutate 3 gene
        int percent_mutation = (rand() % 100);
        if (percent_mutation > 50) {
            a = (rand() % 5);
            inner_chromosome[i][a] = inner_chromosome[i][a] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][a] > maxx) inner_chromosome[i][a] = maxx;
            if (inner_chromosome[i][a] < 0)    inner_chromosome[i][a] = 0.1;
        } else if (percent_mutation > 70) {
            b = (a + (rand() % 4) + 1 % 5);
            inner_chromosome[i][b] = inner_chromosome[i][b] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][b] > maxx) inner_chromosome[i][b] = maxx;
            if (inner_chromosome[i][b] < 0)    inner_chromosome[i][b] = 0.1;
        } else if (percent_mutation > 90) {
            c = (b + (rand() % 3) + 1 % 5);
            inner_chromosome[i][c] = inner_chromosome[i][c] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][c] > maxx) inner_chromosome[i][c] = maxx;
            if (inner_chromosome[i][c] < 0)    inner_chromosome[i][c] = 0.1;
        }
    }
}

// inner selection
void elitism(void) {
    int a, b, c;
    for (int i = 0; i < pop; i++) {
        if (i == maxi) {
            continue;
        }

        // Cross_over
        for (int j = 0; j < 5; j++) {
            if (delta_cross_over) {
                inner_chromosome[i][j] = cross_with_delta(inner_chromosome[i][j], inner_chromosome[maxi][j]);
            } else {
            inner_chromosome[i][j] = (inner_chromosome[i][j] + inner_chromosome[maxi][j])/ 2.0;
            }
        }

        // mutation
        // 50% to mutate 1 gene
        // 30% to mutate 2 gene
        // 10% to mutate 3 gene
        int percent_mutation = (rand() % 100);
        if (percent_mutation > 50) {
            a = (rand() % 5);
            inner_chromosome[i][a] = inner_chromosome[i][a] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][a] > maxx) inner_chromosome[i][a] = maxx;
            if (inner_chromosome[i][a] < 0)    inner_chromosome[i][a] = 0.1;
        } else if (percent_mutation > 70) {
            b = (a + (rand() % 4) + 1 % 5);
            inner_chromosome[i][b] = inner_chromosome[i][b] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][b] > maxx) inner_chromosome[i][b] = maxx;
            if (inner_chromosome[i][b] < 0)    inner_chromosome[i][b] = 0.1;
        } else if (percent_mutation > 90) {
            c = (b + (rand() % 3) + 1 % 5);
            inner_chromosome[i][c] = inner_chromosome[i][c] + ((double) (rand() % maxx)-maxx/2)*mutation_rate*2/1000.0f;
            if (inner_chromosome[i][c] > maxx) inner_chromosome[i][c] = maxx;
            if (inner_chromosome[i][c] < 0)    inner_chromosome[i][c] = 0.1;
        }
    }
}

// cross over with delta
double cross_with_delta(double chromox, double chromoy) {
    int switch_delta;
    float delta;
    double offspring = 0.0;
    delta = chromox + chromoy / 2.0;
    switch_delta = rand() % 7;
    switch (switch_delta) {
        case 0:
            offspring = delta;
            break;
        case 1:
            offspring = chromox;
            break;
        case 2:
            offspring = chromoy;
            break;
        case 3:
            offspring = chromox - delta;
            break;
        case 4:
            offspring = chromoy - delta;
            break;
        case 5:
            offspring = chromox + delta;
            break;
        case 6:
            offspring = chromoy + delta;
            break;
    }
    if (offspring < 0) offspring = 0.01;
    if (offspring > maxx) offspring = maxx;

    return(offspring);
}

void predation(void) {
    if(gen % period_predation == 0) {
        int worst_index = 0;
        int worst_fitness = fit[0];
        for(int i = 1; i < pop; i++) {
            if(fit[i] < worst_fitness) {
                worst_index = i;
                worst_fitness = fit[i];
            }
        }
        for (int j = 0; j < 5; j++) {
            inner_chromosome[worst_index][j] = (double) (rand() % maxx);
        }
    }
}

// Pseudo-synthesize:
// The concept of synthesizing predation is to take the most common gene
// from the entire population.
// However, in this project, all five genes from the inner EA are
// floating-point numbers, and it is almost impossible to have a "common gene"
// since every gene most certainly will have a unique value.
// The pseudo-synthesize will address this issue by taking the average
// from the entire population instead of the most common.
void pseudo_synthesize(void) {
    if(gen % period_synthesis == 0) {
        int worst_index = 0;
        int worst_fitness = fit[0];
        for(int i = 1; i < pop; i++) {
            if(fit[i] < worst_fitness) {
                worst_index = i;
                worst_fitness = fit[i];
            }
        }
        for (int j = 0; j < 5; j++) {
            double gene_mean[5] = {0};
            for (int k = 0; k < pop; k++) {
                gene_mean[j] += inner_chromosome[k][j];
            }
            inner_chromosome[worst_index][j] = (double) gene_mean[j] / pop;
        }
    }
}

// repolulate or "genocide":
// with a proper trigger of stagnation, eradicate the population
// except for the most successful chromosome.
void repopulate(void) {
    for (int i = 0; i < pop; i++) {
        if(i!=maxi) {
            for (int j = 0; j < 5; j++)
                inner_chromosome[i][j] = (double) (rand() % maxx);
        }
    }
}

// evaluate-mutation:
// evaluate the progress of the past five generations and verify the stagnation.
// With the stagnation, make two generations of "fine-tuning"
// by narrowing the mutation rate. After two generations,
// the next generations will expand the mutation rate to
// explore ranges in the gene farther from the most successful chromosome.
void evalutate_mutation(void) {
    if((bestInd[4] - bestInd[0] < 0.00001) && (flag_fine_tune < 2)) {
        mutation_rate = mutation_rate / 10;
        flag_fine_tune++;
    } else if((bestInd[4] - bestInd[0] < 0.00001) && (flag_fine_tune == 2)) {
        mutation_rate *= 120;
        flag_fine_tune++;
    } else if((bestInd[4] - bestInd[0] < 0.00001) && (flag_fine_tune >= 3)) {
        mutation_rate *= 1.2;
        flag_fine_tune++;
    }
    if(mutation_rate > 100) {
        repopulate();
        flag_fine_tune = 0;
        mutation_rate = 2;
    }
}

void elitism2(void) {
    int a, b, c;
    for (int i = 0; i < pop2; i++) {
        if (i == maxi2) {
            continue;
        }
        // case last fitness is better than maxfit2
        if (heritage[i][heritage_max - 1] > maxfit2) {
            continue;
        }

        // mutation
        int mutation_gene = (rand() % 5);
        switch (mutation_gene) {
            case 0:
                gene_mutation_rate[i] += ((double) ((rand() % 1000) - 1000/2) * mutation_rate_squared*2/1000.0f);
                if (gene_mutation_rate[i] > 500) gene_mutation_rate[i] -= 500;
                if (gene_mutation_rate[i] < 0) gene_mutation_rate[i] += 500;
                break;
            case 1:
                gene_selection[i] = !(gene_selection[i]);
                break;
            case 2:
                gene_crossover[i] = !(gene_crossover[i]);
                break;
            case 3:
                gene_predation[i] += ((double) ((rand() % inner_loop) - inner_loop/2) * predation_mutation_rate*2/1000.0f);
                if (gene_predation[i] > inner_loop) gene_predation[i] -= inner_loop;
                if (gene_predation[i] < 0) gene_predation[i] += inner_loop;
                if (gene_predation[i] == 0) gene_predation[i] = 1;
                break;
            case 4:
                gene_syntesis[i] += ((double) ((rand() % inner_loop) - inner_loop/2) * predation_mutation_rate*2/1000.0f);
                if (gene_syntesis[i] > inner_loop) gene_syntesis[i] -= inner_loop;
                if (gene_syntesis[i] < 0) gene_syntesis[i] += inner_loop;
                if (gene_syntesis[i] == 0) gene_syntesis[i] = 1;
                break;
        }
    }
}

void print_result(void) {
    printf("best fitness: %f\n", ultimate_inner_fit);
    for (int i = 0; i < 5; i++) {
        printf("best x[%d]: %f\n", i, ultimate_inner_chromosome[i]);
    }
    printf("best mutation rate: %f%%\n", (gene_mutation_rate[maxi2]/10));
    if (gene_selection[maxi2] == 0) {
        printf("best selection method: elitism\n");
    } else {
        printf("best selection method: tournament\n");
    }
    if (gene_crossover[maxi2] == 0) {
        printf("best cross-over method: mean\n");
    } else {
        printf("best cross-over method: delta\n");
    }
    printf("best random predation interval: %f\n", ((float) gene_predation[maxi2]/ inner_loop ));
    printf("best sythensis predation interval: %f\n", ((float) gene_syntesis[maxi2] / inner_loop ));
}
