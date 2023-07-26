#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define LIN 7
#define COL 7
#define QTD_GERACAO 40000
#define TAM_POPULACAO 100
#define TAM_CROMOSSOMO 49
#define TAXA_SEL 0.9
#define TAXA_CRUZAMENTO 0.8
#define TAXA_MUTACAO 0.8

typedef struct {
    char linha;
    char col;
    char dado;
} Coordenada;

typedef int fixed_point_t;
#define SCALE_FACTOR 1000
typedef struct {
    fixed_point_t inf;
    fixed_point_t sup;
    fixed_point_t porc;
    char p;
} Faixas_roleta;

Coordenada populacao_anterior[TAM_POPULACAO][TAM_CROMOSSOMO] = {0, 0, 0}; //Matriz de cromossomo
Coordenada populacao_nova[TAM_POPULACAO][TAM_CROMOSSOMO] = {0, 0, 0}; //Matriz de cromossomo

unsigned int avaliacao_parcial_populacao[TAM_POPULACAO];  // Matriz de avalia��es parcial
unsigned char indice_notas[TAM_POPULACAO];
Faixas_roleta fx_roleta[TAM_POPULACAO];

Coordenada inicio = {0, 0, 1};
Coordenada final = {5, 3, 39};
int geracao_atual = 0;
unsigned int soma_pesos = 0;
unsigned int maior_valor_soma_pesos = 0;

// Prototipos Funcoes
void embaralha_alelos();
void init_populacao();
void free_mapa();
void print_populacao();
void avaliar_populacao();
unsigned int reavalia(char);
unsigned int distancia(Coordenada, Coordenada);
void print_avaliacao_parcial_populacao();
void ordenar_cromo(unsigned int v[TAM_POPULACAO]);
void reproduzir_nova_geracao();
void get_pais(char* , char*);
void roleta();
void print_roleta();
int check_stop();
void print_arq_rota();
void bests_cromo();
bool cruzapais(char pai_1, char pai_2);





fixed_point_t float_to_fixed_point(float x) {
    return (fixed_point_t)(x * SCALE_FACTOR);
}

float fixed_point_to_float(fixed_point_t x) {
    return ((float)x / SCALE_FACTOR);
}

int main()
{
    srand((unsigned long long)time(NULL) );
    init_populacao();
    print_populacao();

    while(true) {
        avaliar_populacao();
        reproduzir_nova_geracao();

        if (check_stop()){
            avaliar_populacao();
            break;
		}
    }
    print_populacao();
    bests_cromo();
    print_arq_rota();
    printf("\nMaior valor de soma pesos = %llu\n", maior_valor_soma_pesos);
    return 0;
}


int check_stop() {

	return (geracao_atual == QTD_GERACAO - 1);
}

unsigned int reavalia(char j)
{
    return (populacao_anterior[j][0].dado != inicio.dado) ? (avaliacao_parcial_populacao[j]) * 100 : 0;
}

unsigned int distancia(Coordenada inicio, Coordenada atual)
{
    unsigned char lin = abs(inicio.linha - atual.linha);
    unsigned char col = abs(inicio.col - atual.col);
    return pow((lin + col) * 20, 2);
}

void debug(fixed_point_t n) {
    // Essa linha nunca deve ser alcançada.
    printf("Erro: n = %.3f não está em nenhuma faixa\n", n);
    exit(10);
}

void roleta()
{
    char k;
    fixed_point_t acumulado = 0;
    for(k = 0; k < TAM_POPULACAO; k++) {
        fx_roleta[k].porc = float_to_fixed_point((float)avaliacao_parcial_populacao[k] / soma_pesos);
        fx_roleta[k].inf = acumulado;
        acumulado += fx_roleta[k].porc;
        fx_roleta[k].sup = acumulado;
        fx_roleta[k].p = indice_notas[(TAM_POPULACAO - 1) - k];
    }
}

void verifica_repeticoes(unsigned int j_pai, int j_filho, int pt_corte_1, int pt_corte_2)
{
    int k, k1;
    bool cked_1, cked_2;
    do {
        cked_1 = true;
        cked_2 = true;
        for (k = TAM_CROMOSSOMO - 1; k >= pt_corte_2; k--) {
            for (k1 = pt_corte_1; k1 < pt_corte_2; k1++) {
                if (populacao_nova[j_filho][k].dado == populacao_nova[j_filho][k1].dado){//, comparar end ou dado?
                     populacao_nova[j_filho][k] = populacao_anterior[j_pai][k1];
                    cked_1 = false;
                    break;
                }
            }
        }
        for (k = 0; k < pt_corte_1; k++){
            for (k1 = pt_corte_1; k1 < pt_corte_2; k1++){
                if (populacao_nova[j_filho][k].dado == populacao_nova[j_filho][k1].dado){//, comparar end ou dado?
                     populacao_nova[j_filho][k] = populacao_anterior[j_pai][k1];
                    cked_2 = false;
                    break;
                }
            }
        }
    } while(!(cked_1 && cked_2));
}

void transfer_pop()
{
    for (int j = 0; j < TAM_POPULACAO; j++) {
        for (int k = 0; k < TAM_CROMOSSOMO; k++) {
            populacao_anterior[j][k] = populacao_nova[j][k];

        }
    }
}

void mutapais(int j) {
    if(((double)rand() / RAND_MAX) <= TAXA_MUTACAO) {
	    int pt_mt_1 = rand() % TAM_CROMOSSOMO;
	    int pt_mt_2 = rand() % TAM_CROMOSSOMO;
	    Coordenada aux = populacao_nova[j][pt_mt_1];
        populacao_nova[j][pt_mt_1] = populacao_nova[j][pt_mt_2];
        populacao_nova[j][pt_mt_2] = aux;
	}
}

void get_pais(char* pai_1, char* pai_2)
{
    char j;
    fixed_point_t n;
    bool debug_t = false;
    while (*pai_1 == *pai_2) {

        n = float_to_fixed_point((float)rand() / RAND_MAX);
        for (j = 0; j < TAM_POPULACAO; j++) {
            if (n >= fx_roleta[j].inf && (n < fx_roleta[j].sup || j == TAM_POPULACAO - 1)) {
                *pai_1 = fx_roleta[j].p;
                debug_t = true;
            }
        }
        if (!debug_t) debug(n);
        debug_t = false;
        n = float_to_fixed_point((float)rand() / RAND_MAX);
        for (j = 0; j < TAM_POPULACAO; j++) {
            if (n >= fx_roleta[j].inf && (n < fx_roleta[j].sup || j == TAM_POPULACAO - 1)) {
                *pai_2 = fx_roleta[j].p;
                debug_t = true;
            }
        }
        if (!debug_t) debug(n);
        debug_t = false;
    }
}

void reproduzir_nova_geracao() {

	char _i_novapop = 0;
    geracao_atual += 1;
    char pai_1, pai_2;
    // elitismo();

	while(_i_novapop < TAM_POPULACAO) {
        do {
            pai_1 = 0; pai_2 = 0;
            get_pais(&pai_1, &pai_2);

        }while(!cruzapais(pai_1, pai_2));

		mutapais(_i_novapop);
		mutapais(_i_novapop + 1);
		_i_novapop += 2;
	}
	transfer_pop();
}

bool cruzapais(char pai_1, char pai_2) {
	char static id_cruz = 0;

	if (((double)rand() / RAND_MAX) <= TAXA_CRUZAMENTO) {
        char pt_corte_1, pt_corte_2;
        pt_corte_1 = rand() % (TAM_CROMOSSOMO);
        pt_corte_2 = pt_corte_1 + rand() % (TAM_CROMOSSOMO  - pt_corte_1);
        for (char k = 0; k < TAM_CROMOSSOMO; k++) {
            if (k >= pt_corte_1 && k < pt_corte_2) {
                populacao_nova[id_cruz + 1][k] = populacao_anterior[pai_1][k]; //filho mae
                populacao_nova[id_cruz][k] = populacao_anterior[pai_2][k]; //filho pai
            }
            else {
                populacao_nova[id_cruz][k] = populacao_anterior[pai_1][k];
                populacao_nova[id_cruz + 1][k] = populacao_anterior[pai_2][k];
            }
        }
        verifica_repeticoes(pai_1, id_cruz, pt_corte_1, pt_corte_2);
        verifica_repeticoes(pai_2, id_cruz + 1, pt_corte_1, pt_corte_2);
        id_cruz  += 2;
        id_cruz = (id_cruz == TAM_POPULACAO) ? 0 : id_cruz;

        return true;
    }
    return false;
}

void avaliar_populacao() {
    char j, k;
    soma_pesos = 0;
    for (j = 0; j < TAM_POPULACAO; j++) {
        unsigned int peso = 0;
         for (k = 0;  k < TAM_CROMOSSOMO; k++) {
            Coordenada pos = populacao_anterior[j][k]; //cromo-init
            if (k == 0) {
                avaliacao_parcial_populacao[j] = 0;
            }
            else {
                peso = pow(distancia(populacao_anterior[j][k - 1], populacao_anterior[j][k]), 1);
                avaliacao_parcial_populacao[j] += peso;
            }
            if (k != 0 && pos.dado == final.dado) {
                avaliacao_parcial_populacao[j] += k;
                avaliacao_parcial_populacao[j] += reavalia(j);
                break;
            }
        }
        indice_notas[j] = j; //guarda o endereco do cromosso
        soma_pesos += avaliacao_parcial_populacao[j]; //conteudo da nota
    }

    // printf_indice_notas();
    // print_avaliacao_parcial_populacao();
    ordenar_cromo(avaliacao_parcial_populacao);
    // print_avaliacao_parcial_populacao();
    // printf_indice_notas();
    roleta();
    // print_roleta();
    // printf("soma = %llu\n", soma_pesos);
    if (soma_pesos >= maior_valor_soma_pesos) maior_valor_soma_pesos = soma_pesos;
}

void ordenar_cromo(unsigned int v[TAM_POPULACAO])
{
    char i, j;
    int x, p;
	for (j = 1; j < TAM_POPULACAO; ++j) {
		x = v[j];
		p = indice_notas[j];
		for (i = j - 1; i >= 0 && v[i] < x; --i) {
            v[i + 1] = v[i];
            indice_notas[i + 1] = indice_notas[i];
		}
		v[i + 1] = x;
		indice_notas[i + 1] = p;
	}
}

void init_populacao() {
    char p, i, j, k, cont = 0;
    for (p = 0; p < TAM_POPULACAO; p++) {
        k = 0; cont = 0;
        for (i = 0; i < LIN; i++) {
            for (j = 0; j < COL; j ++) {
                populacao_anterior[p][k].col = j;
                populacao_anterior[p][k].linha = i;
                populacao_anterior[p][k++].dado = cont++;
            }
        }
    }
}

/// PRINTS

void print_arq_rota()
{
    unsigned long long int k;
    unsigned long long int j;
    FILE *ptr_arq;
    ptr_arq = fopen("result_rota_otimizado.json","w");
    fprintf(ptr_arq , "%s", "[");

    // pos_cromo(indice_notas[TAM_POPULACAO - 1], i_geraativa, &j);
    j = indice_notas[TAM_POPULACAO - 1];
    for (k = 0; populacao_anterior[j][k].dado != final.dado; k++) {
        if (populacao_anterior[j][k].dado == inicio.dado)
            fprintf(ptr_arq , "%d,", populacao_anterior[j][k].dado);
        else fprintf(ptr_arq , "%d,", populacao_anterior[j][k].dado);
    }
    fprintf(ptr_arq , "%d", populacao_anterior[j][k].dado);
    fprintf(ptr_arq , "%s", "]\n");
    fclose(ptr_arq);
}

void print_cromo(int j)
{
    printf("\n[%s [%d] G = %d\n", __func__, j, geracao_atual);
    int k;
    for(k = 0; k < TAM_CROMOSSOMO; k++)
        printf("C.%d.%d dado= %d\n", j, k + 1, populacao_anterior[j][k].dado);
}

void print_cromo_nova(int j)
{
    printf("\n[%s [%d] G = %d\n", __func__, j, geracao_atual);
    int k;
    for(k = 0; k < TAM_CROMOSSOMO; k++)
        printf("C.%d.%d dado= %d [%d][%d]\n", j, k + 1, populacao_nova[j][k].dado,  populacao_nova[j][k].linha, populacao_nova[j][k].col );
}

void bests_cromo()
{
    int j;
    j = indice_notas[TAM_POPULACAO - 1];
    printf("\ni = %d, j = %d\n", geracao_atual, j);
    print_cromo(j);
}

void print_roleta()
{
    printf("\n[%s]\n", __func__);
    int k;
    for(k = 0; k < TAM_POPULACAO; k++) {
        printf("[.%d - porc = %f ; inf = %f ; sup = %f ; p = %d\n", k + 1, fixed_point_to_float(fx_roleta[k].porc),
        fixed_point_to_float(fx_roleta[k].inf), fixed_point_to_float(fx_roleta[k].sup), fx_roleta[k].p);
    }
}

void print_populacao() {
    printf("\n\n**PRINT POPULACAO***\n\n");
    int j, k;
    for(j = 0; j < TAM_POPULACAO; j++) {
        for(k = 0; k < TAM_CROMOSSOMO; k++) {
            printf("G.%d, %d.%d dado= %d, [%d][%d]\n",
                    geracao_atual + 1, j + 1, k + 1, populacao_anterior[j][k].dado,
                     populacao_anterior[j][k].linha, populacao_anterior[j][k].col);
        }
        printf("\n");
    }
	return;
}

void print_avaliacao_parcial_populacao()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < TAM_POPULACAO; i++) {
        printf("[Individuo %d. = %llu\n", i, avaliacao_parcial_populacao[i]);
    }
    printf("Soma dos pesos = %llu\n\n", soma_pesos);
}

void printf_indice_notas()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < TAM_POPULACAO; i++)
    {
       print_cromo(indice_notas[i]);
    }
}
