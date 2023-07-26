#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define ROWS 7
#define COLS 7
#define GEN_COUNT 40000
#define POP_SIZE 100
#define CHROM_SIZE 49
#define SEL_RATE 0.9
#define CROSS_RATE 0.8
#define MUT_RATE 0.8

typedef struct {
    char row;
    char col;
    char data;
} Coordinate;

typedef int fixed_point_t;
#define SCALE_FACTOR 1000
typedef struct {
    fixed_point_t lower_limit;
    fixed_point_t upper_limit;
    fixed_point_t percentage;
    char p;
} Roulette_range;

Coordinate previous_population[POP_SIZE][CHROM_SIZE] = {0, 0, 0}; //Matriz de cromossomo
Coordinate new_population[POP_SIZE][CHROM_SIZE] = {0, 0, 0}; //Matriz de cromossomo

unsigned int partial_evaluation_population[POP_SIZE];  // Matriz de avalia��es parcial
unsigned char scores_index[POP_SIZE];
Roulette_range roulette_fx[POP_SIZE];

Coordinate start = {0, 0, 1};
Coordinate finish = {5, 3, 39};
int current_generation = 0;
unsigned int weight_sum = 0;
unsigned int highest_weight_sum_value = 0;

// Prototipos Funcoes
void shuffle_alleles();
void init_population();
void print_population();
void evaluate_population();
unsigned int reassess(char);
unsigned int calculate_fitness(Coordinate, Coordinate);
void print_partial_evaluation_population();
void sort_chromo(unsigned int v[POP_SIZE]);
void reproduce_new_generation();
void get_parents(char* , char*);
void roulette();
void print_roulette();
int check_stop();
void print_file_route();
void best_chromo();
bool cross_parents(char, char);
void mutate_parents(int j);
void check_repetitions(unsigned int, int, int, int);

fixed_point_t float_to_fixed_point(float x) {
    return (fixed_point_t)(x * SCALE_FACTOR);
}

float fixed_point_to_float(fixed_point_t x) {
    return ((float)x / SCALE_FACTOR);
}

int cont_rept = 0;

int main()
{
    srand((unsigned long long)time(NULL) );
    init_population();
    print_population();

    while(true) {
        evaluate_population();
        reproduce_new_generation();

        if (check_stop()){
            evaluate_population();
            break;
		}
    }
    print_population();
    best_chromo();
    print_file_route();
    printf("\nMaior valor de soma pesos = %llu\n", highest_weight_sum_value);
    printf("Cont rept = %d\n", cont_rept);
    return 0;
}

int check_stop() {

	return (current_generation == GEN_COUNT - 1);
}

unsigned int reassess(char j)
{
    return (previous_population[j][0].data != start.data) ? (partial_evaluation_population[j]) * 100 : 0;
}

unsigned int calculate_fitness(Coordinate start, Coordinate atual)
{
    unsigned char lin = abs(start.row - atual.row);
    unsigned char col = abs(start.col - atual.col);
    return pow((lin + col) * 20, 2);
}

void debug(fixed_point_t n) {
    // Essa row nunca deve ser alcançada.
    printf("Erro: n = %.3f não está em nenhuma faixa\n", n);
    exit(10);
}

void roulette()
{
    char k;
    fixed_point_t acumulado = 0;
    for(k = 0; k < POP_SIZE; k++) {
        roulette_fx[k].percentage = float_to_fixed_point((float)partial_evaluation_population[k] / weight_sum);
        roulette_fx[k].lower_limit = acumulado;
        acumulado += roulette_fx[k].percentage;
        roulette_fx[k].upper_limit = acumulado;
        roulette_fx[k].p = scores_index[(POP_SIZE - 1) - k];
    }
}

void transfer_pop()
{
    for (int j = 0; j < POP_SIZE; j++) {
        for (int k = 0; k < CHROM_SIZE; k++) {
            previous_population[j][k] = new_population[j][k];

        }
    }
}

void mutate_parents(int j) {
    if(((double)rand() / RAND_MAX) <= MUT_RATE) {
	    int pt_mt_1 = rand() % CHROM_SIZE;
	    int pt_mt_2 = rand() % CHROM_SIZE;
	    Coordinate aux = new_population[j][pt_mt_1];
        new_population[j][pt_mt_1] = new_population[j][pt_mt_2];
        new_population[j][pt_mt_2] = aux;
	}
}

void get_parents(char* parent_1, char* parent_2)
{
    char j;
    fixed_point_t n;
    bool debug_t = false;
    while (*parent_1 == *parent_2) {

        n = float_to_fixed_point((float)rand() / RAND_MAX);
        for (j = 0; j < POP_SIZE; j++) {
            if (n >= roulette_fx[j].lower_limit && (n < roulette_fx[j].upper_limit || j == POP_SIZE - 1)) {
                *parent_1 = roulette_fx[j].p;
                debug_t = true;
            }
        }
        if (!debug_t) debug(n);
        debug_t = false;
        n = float_to_fixed_point((float)rand() / RAND_MAX);
        for (j = 0; j < POP_SIZE; j++) {
            if (n >= roulette_fx[j].lower_limit && (n < roulette_fx[j].upper_limit || j == POP_SIZE - 1)) {
                *parent_2 = roulette_fx[j].p;
                debug_t = true;
            }
        }
        if (!debug_t) debug(n);
        debug_t = false;
    }
}

void reproduce_new_generation() {

	char new_pop = 0;
    current_generation += 1;
    char parent_1, parent_2;
    // elitismo();

	while(new_pop < POP_SIZE) {
        do {
            parent_1 = 0; parent_2 = 0;
            get_parents(&parent_1, &parent_2);

        }while(!cross_parents(parent_1, parent_2));

		mutate_parents(new_pop);
		mutate_parents(new_pop + 1);
		new_pop += 2;
	}
	transfer_pop();
}

bool cross_parents(char parent_1, char parent_2) {
	char static crossing_id = 0;

	if (((double)rand() / RAND_MAX) <= CROSS_RATE) {
        char cut_point_1, cut_point_2;
        cut_point_1 = rand() % (CHROM_SIZE);
        cut_point_2 = cut_point_1 + rand() % (CHROM_SIZE  - cut_point_1);
        for (char k = 0; k < CHROM_SIZE; k++) {
            if (k >= cut_point_1 && k < cut_point_2) {
                new_population[crossing_id + 1][k] = previous_population[parent_1][k]; //filho mae
                new_population[crossing_id][k] = previous_population[parent_2][k]; //filho pai
            }
            else {
                new_population[crossing_id][k] = previous_population[parent_1][k];
                new_population[crossing_id + 1][k] = previous_population[parent_2][k];
            }
        }
        check_repetitions(parent_1, crossing_id, cut_point_1, cut_point_2);
        check_repetitions(parent_2, crossing_id + 1, cut_point_1, cut_point_2);
        crossing_id  += 2;
        crossing_id = (crossing_id == POP_SIZE) ? 0 : crossing_id;

        return true;
    }
    return false;
}

void check_repetitions(unsigned int parent_j, int child_j, int cut_point_1, int cut_point_2)
{
    int k, k1, count = 0;
    bool checked_1 = false, checked_2 = false;

    while (!checked_1 || !checked_2) {
        checked_1 = true;
        checked_2 = true;

        for (k = CHROM_SIZE - 1; k >= cut_point_2; k--) {
            for (k1 = cut_point_1; k1 < cut_point_2; k1++) {
                if (new_population[child_j][k].data == new_population[child_j][k1].data) {
                    new_population[child_j][k] = previous_population[parent_j][k1];
                    checked_1 = false;
                    break;
                }
            }
        }

        for (k = 0; k < cut_point_1; k++) {
            for (k1 = cut_point_1; k1 < cut_point_2; k1++) {
                if (new_population[child_j][k].data == new_population[child_j][k1].data) {
                    new_population[child_j][k] = previous_population[parent_j][k1];
                    checked_2 = false;
                    break;
                }
            }
        }
        count++;
    }

    if (count >= cont_rept) cont_rept = count;
}

void evaluate_population() {
    char j, k;
    weight_sum = 0;
    for (j = 0; j < POP_SIZE; j++) {
        unsigned int peso = 0;
         for (k = 0;  k < CHROM_SIZE; k++) {
            Coordinate pos = previous_population[j][k]; //cromo-init
            if (k == 0) {
                partial_evaluation_population[j] = 0;
            }
            else {
                peso = pow(calculate_fitness(previous_population[j][k - 1], previous_population[j][k]), 1);
                partial_evaluation_population[j] += peso;
            }
            if (k != 0 && pos.data == finish.data) {
                partial_evaluation_population[j] += k;
                partial_evaluation_population[j] += reassess(j);
                break;
            }
        }
        scores_index[j] = j; //guarda o endereco do cromosso
        weight_sum += partial_evaluation_population[j]; //conteudo da nota
    }

    // printf_indice_notas();
    // print_partial_evaluation_population();
    sort_chromo(partial_evaluation_population);
    // print_partial_evaluation_population();
    // printf_indice_notas();
    roulette();
    // print_roulette();
    // printf("soma = %llu\n", weight_sum);
    if (weight_sum >= highest_weight_sum_value) highest_weight_sum_value = weight_sum;
}

void sort_chromo(unsigned int v[POP_SIZE])
{
    char i, j;
    int x, p;
	for (j = 1; j < POP_SIZE; ++j) {
		x = v[j];
		p = scores_index[j];
		for (i = j - 1; i >= 0 && v[i] < x; --i) {
            v[i + 1] = v[i];
            scores_index[i + 1] = scores_index[i];
		}
		v[i + 1] = x;
		scores_index[i + 1] = p;
	}
}

void init_population() {
    char p, i, j, k, count = 0;
    for (p = 0; p < POP_SIZE; p++) {
        k = 0; count = 0;
        for (i = 0; i < ROWS; i++) {
            for (j = 0; j < COLS; j ++) {
                previous_population[p][k].col = j;
                previous_population[p][k].row = i;
                previous_population[p][k++].data = count++;
            }
        }
    }
}

/// PRINTS

void print_file_route()
{
    unsigned long long int k;
    unsigned long long int j;
    FILE *ptr_arq;
    ptr_arq = fopen("result_rota_otimizado.json","w");
    fprintf(ptr_arq , "%s", "[");

    // pos_cromo(scores_index[POP_SIZE - 1], i_geraativa, &j);
    j = scores_index[POP_SIZE - 1];
    for (k = 0; previous_population[j][k].data != finish.data; k++) {
        if (previous_population[j][k].data == start.data)
            fprintf(ptr_arq , "%d,", previous_population[j][k].data);
        else fprintf(ptr_arq , "%d,", previous_population[j][k].data);
    }
    fprintf(ptr_arq , "%d", previous_population[j][k].data);
    fprintf(ptr_arq , "%s", "]\n");
    fclose(ptr_arq);
}

void print_cromo(int j)
{
    printf("\n[%s [%d] G = %d\n", __func__, j, current_generation);
    int k;
    for(k = 0; k < CHROM_SIZE; k++)
        printf("C.%d.%d data= %d\n", j, k + 1, previous_population[j][k].data);
}

void print_cromo_nova(int j)
{
    printf("\n[%s [%d] G = %d\n", __func__, j, current_generation);
    int k;
    for(k = 0; k < CHROM_SIZE; k++)
        printf("C.%d.%d data= %d [%d][%d]\n", j, k + 1, new_population[j][k].data,  new_population[j][k].row, new_population[j][k].col );
}

void best_chromo()
{
    int j;
    j = scores_index[POP_SIZE - 1];
    printf("\ni = %d, j = %d\n", current_generation, j);
    print_cromo(j);
}

void print_roulette()
{
    printf("\n[%s]\n", __func__);
    int k;
    for(k = 0; k < POP_SIZE; k++) {
        printf("[.%d - percentage = %f ; lower_limit = %f ; upper_limit = %f ; p = %d\n", k + 1, fixed_point_to_float(roulette_fx[k].percentage),
        fixed_point_to_float(roulette_fx[k].lower_limit), fixed_point_to_float(roulette_fx[k].upper_limit), roulette_fx[k].p);
    }
}

void print_population() {
    printf("\n\n**PRINT POPULACAO***\n\n");
    int j, k;
    for(j = 0; j < POP_SIZE; j++) {
        for(k = 0; k < CHROM_SIZE; k++) {
            printf("G.%d, %d.%d data= %d, [%d][%d]\n",
                    current_generation + 1, j + 1, k + 1, previous_population[j][k].data,
                     previous_population[j][k].row, previous_population[j][k].col);
        }
        printf("\n");
    }
	return;
}

void print_partial_evaluation_population()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < POP_SIZE; i++) {
        printf("[Individuo %d. = %llu\n", i, partial_evaluation_population[i]);
    }
    printf("Soma dos pesos = %llu\n\n", weight_sum);
}

void printf_indice_notas()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < POP_SIZE; i++)
    {
       print_cromo(scores_index[i]);
    }
}
