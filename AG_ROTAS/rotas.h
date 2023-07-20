#ifndef _ROTA_H
#define _ROTA_H
#include <stdbool.h>
//#include <gmp.h>

#define LIN 9
#define COL 9

#define TRUE 1
#define FALSE 0
#define QTGERA 30000
#define TAMPOP 80
#define TAMCROMO 81 //BUG
#define TAXASEL 0.9
#define TAXACRUZ 1
#define TAXAMUTA 1


//#define INSTALL_DEBUG

typedef struct {
    int linha;
    int col;
    int dado;
}posicao;

typedef struct {
    double inf;
    double sup;
    double porc;
    posicao** p;
}faixas_roleta;

void init_mapa();
void criapop(void);
void embaralha_alelos(int i);
void avaliapop(void);
void ordenar_cromo(unsigned long long int v[][TAMPOP]);
void pos_cromo(posicao** cromo, int g, int * j);
void roleta();
unsigned long long int dis(posicao * inicio, posicao * atual);
posicao** selecionapais();
void reproduzpop(void);
void elitismo();
bool cruzapais(posicao** pai_1, posicao** pai_2);
void mutapais(int j);
void verifica_repeticoes(int j_pai, int j_filho, int pt_corte_1, int pt_corte_2);
int checaparada(void);
void mostrapop(void);
void motrar_nota_geracao();
unsigned long long int reavalia(int _final, int j);
void print_avaliacao_parcial_populacao();

#endif
