#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define LIN 7
#define COL 7
#define QTD_GERACAO 40000
#define TAM_POPULACAO 60
#define TAM_CROMOSSOMO LIN * COL
#define TAXA_SEL 0.9
#define TAXA_CRUZAMENTO 1
#define TAXA_MUTACAO 1

typedef struct {
    int linha;
    int col;
    int dado;
} Coordenada;

typedef struct {
    double inf;
    double sup;
    double porc;
    unsigned int p;
} Faixas_roleta;

Coordenada populacao_anterior [TAM_POPULACAO][TAM_CROMOSSOMO + 1] = {0, 0, 0}; //Matriz de cromossomo
Coordenada populacao_nova [TAM_POPULACAO][TAM_CROMOSSOMO + 1] = {0, 0, 0}; //Matriz de cromossomo

unsigned long long int avaliacao_parcial_populacao[TAM_POPULACAO];  // Matriz de avaliações parcial
unsigned int indice_notas[TAM_POPULACAO];
Faixas_roleta fx_roleta[TAM_POPULACAO];


// Prototipos Funcoes
void print_mapa();
void init_mapa();
void embaralha_alelos();
void init_populacao();
void free_mapa();
void print_populacao();
void avaliar_populacao();
unsigned long long int reavalia(int , int);
unsigned long long int distancia(Coordenada, Coordenada);
void print_avaliacao_parcial_populacao();
void ordenar_cromo(unsigned long long int v[][TAM_POPULACAO]);
void reproduzir_nova_geracao();
unsigned int get_pais();
void roleta();
void print_roleta();
int check_stop();


Coordenada inicio = {0, 0, 1};
Coordenada final = {5, 5, 41};
Coordenada** mapa;
int geracao_atual = 0;
unsigned long long int soma_pesos = 0;




int main()
{
    init_mapa();
    print_mapa();
    init_populacao();
    print_populacao();

    while(true) {
        avaliar_populacao();
        reproduzir_nova_geracao();
        // print_populacao();
        // system("pause");
        if (check_stop()){
            avaliar_populacao();
            break;
		}
    }
    return 0;
}

int check_stop() {

	return (geracao_atual == QTD_GERACAO - 1);
}

unsigned long long int reavalia(int _final, int j)
{
    unsigned long long int peso = 0;
    return peso = (populacao_anterior[j][0].dado != inicio.dado) ? (avaliacao_parcial_populacao[j]) * 100 : peso;
}

unsigned long long int distancia(Coordenada inicio, Coordenada atual)
{
    unsigned int lin = sqrt(pow((inicio.linha - atual.linha), 2));
    unsigned int col = sqrt(pow((inicio.col - atual.col), 2));
    return pow((lin + col) * 20, 2);
}

void print_avaliacao_parcial_populacao()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < TAM_POPULACAO; i++) {
        printf("[Individuo %d. = %llu\n", i, avaliacao_parcial_populacao[geracao_atual % 2][i]);
    }
    printf("Soma dos pesos = %llu\n\n", soma_pesos);
}

void printf_indice_notas()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < TAM_POPULACAO; i++)
    {
        printf(".%d - Nota = %d\n", i, indice_notas[i]);
    }
}

unsigned int get_pais()
{
    double n = rand() / (double)RAND_MAX;
    for(int i = 0; i < TAM_POPULACAO; i++) {
        if(n > fx_roleta[i].inf && n <= fx_roleta[i].sup) {
            return fx_roleta[i].p;
        }
    }
    return NULL;
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
                if (populacao[geracao_atual % 2][j_filho][k].dado == populacao[geracao_atual % 2][j_filho][k1].dado){//go do, comparar end ou dado?
                     populacao[geracao_atual % 2][j_filho][k] = populacao[(geracao_atual % 2) - 1][j_pai][k1];
                    cked_1 = false;
                    break;
                }
            }
        }
        for(k = 0; k < pt_corte_1; k++){
            for(k1 = pt_corte_1; k1 < pt_corte_2; k1++){
                if(populacao[geracao_atual % 2][j_filho][k].dado == populacao[geracao_atual % 2][j_filho][k1].dado){//go do, comparar end ou dado?
                     populacao[geracao_atual % 2][j_filho][k] = populacao[(geracao_atual % 2) - 1][j_pai][k1];
                    cked_2 = false;
                    break;
                }
            }
        }
    } while(!(cked_1 && cked_2));
}

bool cruzapais(unsigned int pai_1, unsigned int pai_2) {
	int static id_cruz = 0;
	if (((double)rand() / RAND_MAX)<=TAXA_CRUZAMENTO) {
        int j_pai = pai_1, k, l_mae = pai_2;
        int pt_corte_1, pt_corte_2;
        pt_corte_1 = rand() % (TAM_CROMOSSOMO);
        pt_corte_2 = pt_corte_1 + rand() % (TAM_CROMOSSOMO  - pt_corte_1);
        for (k = 0; k < TAM_CROMOSSOMO; k++) {
            if (k >= pt_corte_1 && k < pt_corte_2) {
                populacao_nova[id_cruz + 1][k] = populacao_anterior[j_pai][k]; //filho mae
                populacao_nova[id_cruz][k] = populacao_anterior[l_mae][k]; //filho pai
            }
            else {
                populacao_nova[id_cruz][k] = populacao_anterior[j_pai][k];
                populacao_nova[id_cruz + 1][k] = populacao_anterior[l_mae][k];
            }
        }
        verifica_repeticoes(j_pai, id_cruz, pt_corte_1, pt_corte_2);
        verifica_repeticoes(l_mae, id_cruz + 1, pt_corte_1, pt_corte_2);


        id_cruz  += 2;
        id_cruz = (id_cruz == TAM_POPULACAO) ? 0 : id_cruz;
        return true;
    }
    return false;
}

void mutapais(int j) {
    if(((double)rand() / RAND_MAX) <= TAXA_MUTACAO) {
	    int pt_mt_1 = rand() % TAM_CROMOSSOMO;
	    int pt_mt_2 = rand() % TAM_CROMOSSOMO;
	    Coordenada aux = populacao[geracao_atual % 2][j][pt_mt_1];
        populacao[geracao_atual % 2][j][pt_mt_1] = populacao[geracao_atual % 2][j][pt_mt_2];
        populacao[geracao_atual % 2][j][pt_mt_2] = aux;
	}
}

void reproduzir_nova_geracao() {

	int _i_novapop = 0;
    geracao_atual += 1;
    unsigned int i_pai1_, i_pai2_;
    // elitismo();

	while(_i_novapop < TAM_POPULACAO) {
        do {
            while(!((i_pai1_ = get_pais()) != NULL));
            while(!((i_pai2_ = get_pais()) != NULL));
            do {
                if (i_pai2_ == i_pai1_)
                    while(!((i_pai2_ = get_pais()) != NULL));
                else
                    break;
            }while(!(i_pai2_ != i_pai1_));

        } while(!cruzapais(i_pai1_, i_pai2_));
		mutapais(_i_novapop);
		mutapais(_i_novapop + 1);
		_i_novapop+=2;
	}
}


void avaliar_populacao() {
    int j, k;
    unsigned long long int peso = 0;
    soma_pesos = 0;
    for (j = 0; j < TAM_POPULACAO; j++) {
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
                avaliacao_parcial_populacao[j] += reavalia(k, j);
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
    roleta(soma_pesos);
    //print_roleta();
}

void print_roleta()
{
    printf("\n[%s]\n", __func__);
    int k;
    for(k = 0; k < TAM_POPULACAO; k++) {
        printf("[.%d - porc = %.2f ; inf = %.2f ; sup = %.2f ; p = %d\n", k + 1, fx_roleta[k].porc, fx_roleta[k].inf, fx_roleta[k].sup, fx_roleta[k].p);
    }
}

void roleta()
{
    int k;
    for(k = 0; k < TAM_POPULACAO; k++) {
        fx_roleta[k].porc = (double)avaliacao_parcial_populacao[k] / soma_pesos;
        fx_roleta[k].inf = k == 0 ? 0 : fx_roleta[k - 1].sup ;
        fx_roleta[k].sup = fx_roleta[k].inf + (fx_roleta[k].porc);
        fx_roleta[k].p = indice_notas[(TAM_POPULACAO - 1) - k];
    }
}

void ordenar_cromo(unsigned long long int v[TAM_POPULACAO])
{
    int i, j, x, p;
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

void print_mapa()
{
    printf("\n\n*** PRINT MAPA ***\n\n");
    int i, j;
    for(i = 0; i < LIN; i++) {
        for(j = 0; j < COL; j++) {
            printf("[%d][%d] - %d, ", mapa[i][j].linha, mapa[i][j].col, mapa[i][j].dado);
        }
        printf("\n");
    }
}

void free_mapa()
{
    int i;

    for(i = 0; i < LIN; i++) {
        free(mapa[i]);
    }

    free(mapa);
}

void init_mapa()
{
    int i, j, cont = 0;

    mapa = malloc(LIN * sizeof(Coordenada*));

    for (i = 0; i < LIN; i++) {
        mapa[i] = malloc(COL * sizeof(Coordenada));
        for (j = 0; j < COL; j++) {
            mapa[i][j].dado = ++cont;
            mapa[i][j].linha = i;
            mapa[i][j].col = j;
        }
    }
}

void init_populacao() {
    int k = 0;
    for (int p = 0; p < TAM_POPULACAO; p++) {
        for (int i = 0; i < LIN; i++) {
            for (int j = 0; j < COL; j ++) {
                populacao_anterior[p][k].col = mapa[i][j].col;
                populacao_anterior[p][k].linha = mapa[i][j].linha;
                populacao_anterior[p][k++].dado = mapa[i][j].dado;

                populacao_nova[p][k].col = mapa[i][j].col;
                populacao_nova[p][k].linha = mapa[i][j].linha;
                populacao_nova[p][k++].dado = mapa[i][j].dado;
            }
        }
        k = 0;
    }
    free_mapa();
}

void print_populacao() {
    printf("\n\n**PRINT POPULACAO***\n\n");
    int j, k;
    for(j = 0; j < TAM_POPULACAO; j++) {
        for(k = 0; k < TAM_CROMOSSOMO; k++) {
            printf("G.%d, %d.%d dado= %d, [%d][%d]\n",
                    (geracao_atual % 2) + 1, j + 1, k + 1, populacao[geracao_atual % 2][j][k].dado,
                     populacao[geracao_atual % 2][j][k].linha, populacao[geracao_atual % 2][j][k].col);
        }
        printf("\n");
    }
	return;
}

