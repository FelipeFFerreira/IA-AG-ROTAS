/*
 * AlgoritmoGenetico.c
 *
 * Programa que simula um algoritmo gen�tico.
 *
 * Felipe Ferreira Nascimento (Ci�ncia da Computa��o)
 * Gabriel Romano Godoi Pereira (Ci�ncia da Computa��o)
 * Jaime Mathias de Lara Bueno (Ci�ncia da Computa��o)
 * Marcus Vinicius de Souza Olimpio da Silva (Ci�ncia da Computa��o)
 * Willy Pestana Filho (Ci�ncia da Computa��o)
 *
 * Disciplina: Intelig�ncia Artificial II
 *
 * Professor: Marcio Luiz Piva
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "rotas.h"

int main(void)
{
srand((unsigned long long)time(NULL) );
	init_mapa();
	print_mapa();
	criapop();
	mostrapop();

	while(TRUE) {
		avaliapop();
		system("pause");
		reproduzpop();
		//mostrapop();
		if(checaparada()){
            avaliapop();
            break;
		};
	}
	mostrapop();
    //bests_cromo();
    print_arq_rota();
	return 0;
}
