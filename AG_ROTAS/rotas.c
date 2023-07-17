#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "rotas.h"

faixas_roleta fx_roleta[TAMPOP];
posicao inicio = {0, 0, 1};
posicao final = {5, 5, 41};
int i_geraativa = 0; //geracao atual
posicao* m_i_pop_anterior[1][TAMPOP][TAMCROMO] = {0, 0, 0}; //matriz de cromossomo
posicao* m_i_pop_nova[1][TAMPOP][TAMCROMO] = {0, 0, 0}; //matriz de cromossomo
unsigned long long int m_f_popaval[1][TAMPOP]; //matriz de avalia��es
posicao** indice_notas[TAMPOP];
float m_f_estataval[QTGERA][3]; //matriz de estatisticas: fitness m�nimo, m�ximo, m�dio
posicao** i_pai1; //primeiro pai selecionado
posicao** i_pai2; //segundo pai selecionado
unsigned long long int soma_pesos = 0;
posicao matriz[LIN][COL];

/* Função responsavel por inicializar o mapa com as respectivas cordenadas de cada posição */
void init_mapa()
{
    int i, j, cont = 0;
    for(i = 0; i < LIN; i++) {
        for(j = 0; j < COL; j++) {
            matriz[i][j].dado = ++cont;
            matriz[i][j].linha = i;
            matriz[i][j].col = j;
        }
    }
}

void print_mapa()
{
    printf("\n\n*** PRINT MAPA ***\n\n");
    int i, j;
    for(i = 0; i < LIN; i++) {
        for(j = 0; j < COL; j++) {
            printf("[%d][%d] - %d, ", matriz[i][j].linha, matriz[i][j].col, matriz[i][j].dado);
        }
        printf("\n");
    }
}

/*Inicia população com 2 individuos
 Um individuou contem todo o mapa completo
 Os dois individuos sao criados iguais
*/
void criapop(void) {
    int i, j, k = 0, i_, j_;
    for (i = 0; i < 1; i++) { //for gerecao
        for (j = 0; j < TAMPOP; j++) { //for cromosso
            for (i_ = 0; i_ < LIN; i_++) {
                for (j_ = 0; j_ < COL; j_++) {
                    m_i_pop_anterior[i][j][k++] = &matriz[i_][j_];
                }
            }
            k = 0;
        }
        // embaralha_alelos(i);
    }
}

void embaralha_alelos(int i)
{
    int j, k;
    for(j = 0; j < TAMPOP; j++) {
        for(k = 0; k < TAMCROMO; k++) {
            posicao * troca_aux = m_i_pop_anterior[i][j][k];
            int troca = rand() % TAMCROMO;
            m_i_pop_anterior[i][j][k] = m_i_pop_anterior[i][j][troca];
            m_i_pop_anterior[i][j][troca] = troca_aux;
        }
    }
}

unsigned long long int dis(posicao * inicio, posicao * atual)
{
    unsigned int lin = sqrt(pow((inicio->linha - atual->linha), 2));
    unsigned int col = sqrt(pow((inicio->col - atual->col), 2));
    return pow((lin + col) * 20, 2);
}

void print_avaliacao_parcial_populacao()
{
    printf("\n[%s]\n", __func__);
    for (int i = 0; i < TAMPOP; i++) {
        printf("[Individuo %d. = %llu\n", i, m_f_popaval[0][i]);
    }
    printf("Soma dos pesos = %llu\n\n", soma_pesos);
}

void print_roleta()
{
    printf("\n[%s]\n", __func__);
    int k, j;
    for(k = 0; k < TAMPOP; k++) {
        pos_cromo(fx_roleta[k].p, 0, &j);
        printf("[.%d - porc = %.2f ; inf = %.2f ; sup = %.2f ; p = %d\n", k + 1, fx_roleta[k].porc, fx_roleta[k].inf, fx_roleta[k].sup, j);
    }
}

void avaliapop(void) {
    int j, k;
    unsigned long long int peso = 0;
    soma_pesos = 0;
    for (j = 0; j < TAMPOP; j++) {
         for (k = 0;  k < TAMCROMO; k++) {
            posicao* pos = m_i_pop_anterior[0][j][k]; //cromo-init
            if (k == 0) {
                m_f_popaval[0][j] = 0;
            }
            else {
                peso = pow(dis(m_i_pop_anterior[0][j][k - 1], m_i_pop_anterior[0][j][k]), 1);
                m_f_popaval[0][j] += peso;
            }
            if (k != 0 && pos->dado == final.dado) {
                m_f_popaval[0][j] += k;
                m_f_popaval[0][j] += reavalia(k, j);
                break;
            }
        }
        indice_notas[j] = &(m_i_pop_anterior[0][j][0]); //guarda o endereco do cromosso
        soma_pesos += m_f_popaval[0][j]; //conteudo da nota
    }
        int j2;
        if (i_geraativa >= 30000) {
             //printf("\n_DEBUG_\n");
            // int j;
             //pos_cromo(&m_i_pop[i_geraativa][j][0], i_geraativa, &j2);
             //print_cromo(j2, i_geraativa);
             //system("pause");
        }

    // print_avaliacao_parcial_populacao();
    ordenar_cromo(m_f_popaval);
    // print_avaliacao_parcial_populacao();
    // printf_indice_notas();
    roleta(soma_pesos);
   // print_roleta();
}

unsigned long long int reavalia(int _final, int j)
{
    unsigned long long int peso = 0;
    return peso = (m_i_pop_anterior[0][j][0]->dado != inicio.dado) ? (m_f_popaval[0][j]) * 100 : peso;
}

void ordenar_cromo(unsigned long long int v[][TAMPOP])
{
    int i, j, x;
    posicao * p;
	for (j = 1; j < TAMPOP; ++j) {
		x = v[0][j];
		p = indice_notas[j];
		for (i = j - 1; i >= 0 && v[0][i] < x; --i) {
            v[0][i + 1] = v[0][i];
            indice_notas[i + 1] = indice_notas[i];
		}
		v[0][i + 1] = x;
		indice_notas[i + 1] = p;
	}
}

void bests_cromo()
{
    int j;
    pos_cromo(indice_notas[TAMPOP - 1], 0, &j);
    printf("\ni = %d, j = %d\n", i_geraativa, j);
    print_cromo(j, 0);
    //exit(111);
}

void pos_cromo(posicao** cromo, int g, int * j)
{
    int i;
    for(i = 0; i < TAMPOP; i++) {
        posicao** pos  = &m_i_pop_anterior[0][i][0];
        if(pos == cromo) {
            *j = i;
             return;
        }
    }
}

void printf_indice_notas()
{
    printf("\n[%s]\n", __func__);
    int j;
    for (int i = 0; i < TAMPOP; i++)
    {
        pos_cromo(indice_notas[i], 0, &j);
        print_cromo(j, i_geraativa);
    }
}
void roleta()
{
    int k;
    for(k = 0; k < TAMPOP; k++) {
        fx_roleta[k].porc = (double)m_f_popaval[0][k] / soma_pesos;
        fx_roleta[k].inf = k == 0 ? 0 : fx_roleta[k - 1].sup ;
        fx_roleta[k].sup = fx_roleta[k].inf + (fx_roleta[k].porc);
        fx_roleta[k].p = indice_notas[(TAMPOP - 1) - k];
    }
}

posicao** selecionapais()
{
    double n = rand() / (double)RAND_MAX;
    // printf("\n[%s] n = %.3f\n", __func__, n);
    int j;
    for (int i = 0; i < TAMPOP; i++) {
        if (n > fx_roleta[i].inf && n <= fx_roleta[i].sup) {
            pos_cromo(fx_roleta[i].p, 0, &j);
            // printf("\n[%s] fx_roleta = %d\n", __func__, j);
            return fx_roleta[i].p;
        }
    }
    return NULL;
}

void print_cromo(int j, int i)
{
    printf("\n[%s [%d] G = %d\n", __func__, j, i);
    int k;
    for(k = 0; k < TAMCROMO; k++)
        printf("C.%d.%d dado= %d\n", j, k + 1, m_i_pop_anterior[0][j][k]->dado);
}

void transfer_pop()
{
    for (int j = 0; j < TAMPOP; j++) {
        for (int k = 0; k < TAMCROMO; k++) {
            m_i_pop_anterior[0][j][k] = m_i_pop_nova[0][j][k];

        }
    }
}
void reproduzpop(void) {

	int _i_novapop = 0, j1, j2;
    i_geraativa += 1;
    posicao** i_pai1_;
    posicao** i_pai2_;
    // elitismo();

	while(_i_novapop < TAMPOP) {
        do {
            while(!((i_pai1_ = selecionapais()) != NULL));
            while(!((i_pai2_ = selecionapais()) != NULL));
            do {
                if(i_pai2_ == i_pai1_)
                    while(!((i_pai2_ = selecionapais()) != NULL));
                else
                    break;
            } while(!(i_pai2_ != i_pai1_));
            // printf("\n_DEBUG_\n");
            // pos_cromo(i_pai1_, i_geraativa - 1, &j1);
            // pos_cromo(i_pai2_, i_geraativa - 1, &j2);
            // print_cromo(j1, i_geraativa - 1);
            // print_cromo(j2, i_geraativa - 1);
            // system("pause");

        } while(!cruzapais(i_pai1_, i_pai2_));
		mutapais(_i_novapop);
		mutapais(_i_novapop + 1);
		_i_novapop+=2;
	}

	transfer_pop();
}

void elitismo()
{
    int j1, j2, k;
    pos_cromo(indice_notas[TAMPOP - 1], 0 - 1, &j1);
    pos_cromo(indice_notas[TAMPOP - 2], 0 - 1, &j2);
    for(k = 0; k < TAMCROMO; k++) {
        posicao * pos1 = m_i_pop_anterior[0][j1][k];
        posicao * pos2 = m_i_pop_anterior[0][j2][k];
        m_i_pop_anterior[0][0][k] = pos1;
        m_i_pop_anterior[0][1][k] = pos2;
    }
}

bool cruzapais(posicao** pai_1, posicao** pai_2) {
	int static id_cruz = 0;
	if(((double)rand() / RAND_MAX)<=TAXACRUZ) {
        int j_pai, k, l_mae;
        int pt_corte_1, pt_corte_2;
        pt_corte_1 = rand() % (TAMCROMO);
        pt_corte_2 = pt_corte_1 + rand() % (TAMCROMO  - pt_corte_1);
        pos_cromo(pai_1, 0, &j_pai);
        pos_cromo(pai_2, 0, &l_mae);
        for(k = 0; k < TAMCROMO; k++) {
            if(k >= pt_corte_1 && k < pt_corte_2) {
                m_i_pop_nova[0][id_cruz + 1][k] = m_i_pop_anterior[0][j_pai][k]; //filho mae
                m_i_pop_nova[0][id_cruz][k] = m_i_pop_anterior[0][l_mae][k]; //filho pai
            }
            else {
                m_i_pop_nova[0][id_cruz][k] = m_i_pop_anterior[0][j_pai][k];
                m_i_pop_nova[0][id_cruz + 1][k] = m_i_pop_anterior[0][l_mae][k];
            }
        }
        verifica_repeticoes(j_pai, id_cruz, pt_corte_1, pt_corte_2);
        verifica_repeticoes(l_mae, id_cruz + 1, pt_corte_1, pt_corte_2);
        // *j1 = id_cruz;
        // *j2 = id_cruz + 1;
        // if (i_geraativa == 30000) {
        //     printf("\n_DEBUG_\n");
        //     int j;
        //     pos_cromo(&m_i_pop[i_geraativa][id_cruz][0], i_geraativa, &j);
        //     print_cromo(j, i_geraativa);
        //     system("pause");
        // }

        id_cruz  += 2;
        id_cruz = (id_cruz == TAMPOP) ? 0 : id_cruz;

        return true;
    }
    return false;
}

void verifica_repeticoes(int j_pai, int j_filho, int pt_corte_1, int pt_corte_2)
{
    int k, k1;
    bool cked_1, cked_2;
    do {
        cked_1 = true;
        cked_2 = true;
        for(k = TAMCROMO - 1; k >= pt_corte_2; k--) {
            for(k1 = pt_corte_1; k1 < pt_corte_2; k1++) {
                if(m_i_pop_nova[0][j_filho][k]->dado == m_i_pop_nova[0][j_filho][k1]->dado){//go do, comparar end ou dado?
                     m_i_pop_nova[0][j_filho][k] = m_i_pop_anterior[0][j_pai][k1];
                    cked_1 = false;
                    break;
                }
            }
        }
        for(k = 0; k < pt_corte_1; k++){
            for(k1 = pt_corte_1; k1 < pt_corte_2; k1++){
                if(m_i_pop_nova[0][j_filho][k]->dado == m_i_pop_nova[0][j_filho][k1]->dado){//go do, comparar end ou dado?
                     m_i_pop_nova[0][j_filho][k] = m_i_pop_anterior[0][j_pai][k1];
                    cked_2 = false;
                    break;
                }
            }
        }
    } while(!(cked_1 && cked_2));
}

void mutapais(int j) {
    if(((double)rand() / RAND_MAX) <= TAXAMUTA) {
	    int pt_mt_1 = rand() % TAMCROMO;
	    int pt_mt_2 = rand() % TAMCROMO;
	    posicao* aux = m_i_pop_nova[0][j][pt_mt_1];
        m_i_pop_nova[0][j][pt_mt_1] = m_i_pop_nova[0][j][pt_mt_2];
        m_i_pop_nova[0][j][pt_mt_2] = aux;
	}
}

int checaparada(void){

	return (i_geraativa == QTGERA - 1);
}

void mostrapop(void) {
    printf("\n\n**PRINT POPULACAO***\n\n");
    int j, k;
    for(j = 0; j < TAMPOP; j++) {
        for(k = 0; k < TAMCROMO; k++) {
            printf("G.%d, %d.%d dado= %d, [%d][%d]\n",
                    i_geraativa + 1, j + 1, k + 1, m_i_pop_anterior[0][j][k]->dado,
                     m_i_pop_anterior[0][j][k]->linha, m_i_pop_anterior[0][j][k]->col);
        }
        printf("\n");
    }
	return;
}

void print_arq_rota()
{
    unsigned long long int k;
    unsigned long long int j;
    FILE *ptr_arq;
    ptr_arq = fopen("result_rota.json","w");
    fprintf(ptr_arq , "%s", "[");

    pos_cromo(indice_notas[TAMPOP - 1], 0, &j);
    for (k = 0; m_i_pop_anterior[0][j][k]->dado != final.dado; k++) {
        if (m_i_pop_anterior[0][j][k]->dado == inicio.dado)
            fprintf(ptr_arq , "%d,", m_i_pop_anterior[0][j][k]->dado);
        else fprintf(ptr_arq , "%d,", m_i_pop_anterior[0][j][k]->dado);
    }
    fprintf(ptr_arq , "%d", m_i_pop_anterior[0][j][k]->dado);
    fprintf(ptr_arq , "%s", "]\n");
    fclose(ptr_arq);
}
