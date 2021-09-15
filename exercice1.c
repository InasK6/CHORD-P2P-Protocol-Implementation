/*
AUTEURS: Inas KACI 3873337
	Amandine Ta 3300320		

*/

/*
Commandes de compilation et d'exécution:
mpicc -o exo1 exercice1.c -lm
mpirun -np 9  --oversubscribe ./exo1
*/

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <time.h>
#include <math.h>


#include <stdbool.h>
#include <limits.h>
/*
M est une constante qui définit le nombre de bits d'un identifiant de pair ou une clé d'une donnée
*/

#define M 4
/*
N est une constante qui définit le nombre de paires dans le réseau pair à pair 
*/
#define N 8
#define KEY 0
#define INITIATOR 1


#define TAGINIT 0 // tag du message d'initialisation
#define LOOKUP 1  // tag des messages de recherche du noeud responsable d'une clé
#define END 2	// tag signalant la fin de la recherche d'un noeud
#define LASTCHANCE 3 // cas de test si le successeur est responsable de la clé
/******************************************************************************/
/*
variables et structures de données utiles:

*/

int range; // borne supérieure des valeurs possibles des clés

struct Identity{
	int MPI_ID;
	int CHORD_ID;
}
;
/*chaque pair dispose de :*/
struct Identity rang; // son identifiant CHORD et MPI
int chord_finger[M]; //identifiants CHORD fingers du processus courant 
int mpi_finger[M]; //identifiants MPI fingers du processus courant 
int successor; // id CHORD du successeur

/******************************************************************************/

/******************************************************************************/
/* 
 * Fonctions Outils
 */

/**
fonction de hashage aléatoire  garantissant l'unicité des valeurs tirés 

entrée:
	id: identifiant MPI d'une paire ou identifiant d'une donnée
retour: 
	identifiant unique dans l'ensemble des paires, ou dans l'ensemble des données
	cette ce valeur est dans {0, ..., 2^M -1 }
*/

bool exists(int k, int tab[], int NB ){
	bool res= false;
	int i;
	for (i=0; i<NB; i++){
		res=(res || (tab[i]==k));
	
	}
	return res;
}

bool test_exists(){
	int tab[3]={1, 2, 3};
	printf("%d exists in tab: %d\n", 1, exists(1, tab, 3));
	printf("%d doesn't exist in tab: %d\n", 4, exists(4, tab, 3));
}
int hash_func(int id, int tab[], int NB){
	// tirage d'une valeur aléatoire
	
	int res=rand()%range;
	// boucle utile pour garantir l'unicité dans l'ensemble de clés associées
	while(exists(res, tab, NB)){
		res=rand()%range;
	}
	
	return res;
}
/*
returns true if k \in ]a,b]
*/
bool app(int k,int a,int b){

	if(a<b)
		return k>a && k<=b;
		
	else
		return (k>a && k<range ) || (k>=0 && k<=b);

}
/* retourne la valeur positive associée au modulo*/
int modulo (int val, int Z){
	int mod=val%Z;
	if(mod<0){
		return mod+Z;
	}	
	return mod;
}

void affichage(int identifiers[], int successors[] ){
	printf("--------------------------------------------------------------------------------------------------\n");
	int i;
	
	printf("identifiants CHORD :{");
	for( i=0; i<N+1; i++) {
		printf("%d , ", identifiers[i]);
	}
	printf("\n");
	
	printf("successeurs CHORD :{");
	for( i=0; i<N+1; i++) {
		printf("%d , ", successors[i]);
	}
	printf("}\n");
	

	

}

/******************************************************************************/
/* Tests des fonctions outils
/******************************************************************************/

void test_hash_func(int tab[], int NB){
	int val=pow(2,M);
	printf("tirage d'une valeur dans {0, %d -1}\n", val);
	printf("%d\n", hash_func(0, tab, NB));
}




void test_app(){
	printf("expected 1, res=%d\n", app(3, 2, 5));
	printf("expected 0, res=%d\n", app(3, 4, 5));
	printf("expected 1, res=%d\n", app(3, 2,1));
	printf("expected 0, res=%d\n", app(3, 4,1));
}

int test_modulo(){
	printf("%d\n", modulo(10, 8));
	printf("%d\n", modulo(2, 8));
	printf("%d\n", modulo(-1, 8));
	printf("%d\n", modulo(-3, 8));

}
/******************************************************************************/


/******************************************************************************/
/**
Entrée:
	- node: id CHORD dont on est en train de remplir la table finger
	- k : le numéro du finger [0, M-1]

Retour: 
	- numéro du noeud associé au finger[k]
**/
struct Identity getFinger(int node, int k, int identifiers[]){
	struct Identity succ;
	succ.CHORD_ID=-1;
	succ.MPI_ID=-1;
	double val=(double)k;
	int tmp= pow(2,val);
	//int tmp=0;
	int value= (node+ tmp) % range;
	
	int distance=INT_MAX;
	
	// la valeur du finger est l'identifiant du 1er noeud qui succède à value
	//chercher cet identifiant dans le tableau
	//on ignore le noeud 0 qui corresoon au simumateur
	int i;
	for (i=1; i<=N; i++){
		//if( identifiers[i]!=node){
			int diff=modulo(identifiers[i]-value, range);
			
			if(diff<distance){
				succ.CHORD_ID=identifiers[i];
				succ.MPI_ID=i;
				distance=diff;
			}
		//}
	
	}
	return succ;

}
/* 
Initialisation de la DHT par le processus 0
   pour chaque processus correspondant à un noeud:
   - lui attribuer un identifiant unique dans Chord
   - allouer et remplir son tableau finger
   - initialiser la variable next: le noeud successeur dans la DHT 

Remarque: le processus 0 se charge de lancer la simulation

 */

void simulation(){
	int identifiers[N+1];
	memset(identifiers, -1, N+1);
	int successors[N+1];
	int i=0;
	// initialisation des identifiants des paires
	identifiers[0]=-1;
	for(i=1; i<=N; i++){
		identifiers[i]=hash_func(i, identifiers, N+1);
		// envoyer la valeur de l'id CHORD au paire associé
		MPI_Send(&identifiers[i], 1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); 
	}
	// noeud initiant la simulation n'a pas de successeur, ne fait pas partie du réseau de paires
	successors[0]=-1;

	int fingers_MPI[N+1][M];
	int fingers_CHORD[N+1][M];
	int j=0;
	MPI_Status status;
	
	struct Identity tmp;
	tmp.MPI_ID=-1;
	tmp.CHORD_ID=-1;
	for(j=0; j<M; j++){
		fingers_MPI[0][j]=-1;
		fingers_CHORD[0][j]=-1;
	}
	for(i=1; i<=N; i++){
		//printf("-------%d-------- \n", identifiers[i]);
		for(j=0; j<M; j++){
			tmp=getFinger(identifiers[i], j, identifiers);
			
			fingers_MPI[i][j]=tmp.MPI_ID;
			fingers_CHORD[i][j]=tmp.CHORD_ID;
			
		}
	}

	
	for(i=1; i<=N; i++){
		successors[i]=fingers_CHORD[i][0];
		MPI_Send(&successors[i], 1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); 
	}
	
	
	affichage(identifiers, successors);
	// affichage des tables des fingers
	for(i=0; i<=N; i++){
		printf("finger %d { " , identifiers[i]);
		for(j=0; j<M; j++){
			printf("%d, ", fingers_CHORD[i][ j]);
		}
		printf("}\n");
	}
	printf("--------------------------------------------------------------------------------------------------\n");
	
	// envoyer les fingers associés pour chaque noeud
	for (i=1; i<=N; i++){
		MPI_Send(fingers_MPI[i], M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); 
		// j'utilise un synchronous send pour eviter d'avoir à renvoyer un message par 
		// processus qui signale que son initialisation a bien été faite
		MPI_Send(fingers_CHORD[i], M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); 
	}
	
	// après la bonne initialisation 
	// Récéption d'un message de tous les processus signalant qu'ils ont bien été initialisés
	int done=0;
	for(i=1; i<=N; i++){
		MPI_Recv(&done, 1, MPI_INT,MPI_ANY_SOURCE, TAGINIT, MPI_COMM_WORLD, &status);
	}
	
	/*  Initiate Lookup */
	// tirer de manière aléatoire un identifiant de paire 
	// tirage selon l'identifiant MPI, donc pas besoin de vérifier l'existence du pair
	// étant donné qu'il y en a N et que leurs identifiants MPI se succédent
	int random=(rand()%N)+1;
	//identifiant CHORD associé:
	int chord_random=identifiers[random];
	// tirage aléatoire d'une clé à chercher
	int key=rand()%range;
	// contenu du message
	int msg[2];
	msg[KEY]=key;
	msg[INITIATOR]=chord_random;
	// envoi de la demande de recherche au pair
	printf("**************************************************************************\n");
	printf("déclenchement d'une recherche de la clé %d à partir du noeud chord %d\n", key, chord_random);
	printf("**************************************************************************\n");
	MPI_Send(msg, 2, MPI_INT, random, LOOKUP, MPI_COMM_WORLD); 
	
	
	//attente de response avec le résultat de la recherche
	int result_node;
	
	MPI_Recv(&result_node, 1, MPI_INT,MPI_ANY_SOURCE, END, MPI_COMM_WORLD, &status);
	printf("le resultat %d reçu par le simulateur\n", result_node);
	//propagation du message de terminaison à tous les processus
	for(i=1; i<=N; i++){
		MPI_Send(&result_node, 1, MPI_INT,i, END, MPI_COMM_WORLD);
	}
	
	
	
}


/******************************************************************************/
/* fonctions exécutées par les Paire*/

/*

initialisation des données d'une paire
*/

void initialisation(){
	MPI_Status status;
	// recevoir son identifiant CHORD du simlateur d'id 0
	MPI_Recv(&rang.CHORD_ID, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
	// reception de l'id CHORD du successeur
	MPI_Recv(&successor, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
	// recevoir la liste de fingers associés à la paire courante
	MPI_Recv(mpi_finger, M , MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
	MPI_Recv(chord_finger, M , MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
	
	// envoi de message signalant la bonne initialisation
	int done=1;
	MPI_Send(&done, 1, MPI_INT,0,TAGINIT, MPI_COMM_WORLD); 
	

}
/*  étant donné une clé, retourne le noeud de plus grande valeur*/ 
int findNextIndex(int k, int initiateur){
	int i;
	bool res=false;
	for(i=M-1; i>=0; i--){
	// faire un tableau pour pouvoir passer 2 valeurs
		if((res=app(k,chord_finger[i], initiateur)) && chord_finger[i]!= initiateur){
			return i;
		}
		printf("res: %d\n", res);
	}
	return -1;

}
#define RES 0
int receive(){
	int value[2];
	MPI_Status status;
	MPI_Recv(value, 2, MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if(status.MPI_TAG==END){
		// value : id CHORD du noeud responsable de la clé
		printf("Noeud {%d, %d} : fin de la recherche de la clé, res: %d\n", rang.CHORD_ID, rang.MPI_ID, value[RES]);
		
		return END;
	}
	if(status.MPI_TAG==LOOKUP){
		/*printf("lookup sur %d , fingers  {",rang.CHORD_ID);
		int i;
		for (i=0; i<M; i++){
			printf("%d, ", chord_finger[i]);
		}
		printf("\n");*/
		
		// value : clé recherchée 
		if (app(value[KEY], rang.CHORD_ID, chord_finger[0])){
			// le successeur est responsable de la clé: 
			// envoyer ce résultat: noeud responsable au simulateur
			printf("demande de lookup reçu par %d, le responsable de la clé est le successeur %d\n", rang.CHORD_ID, successor);
			MPI_Send(&chord_finger[0], 1, MPI_INT, 0, END, MPI_COMM_WORLD);
			
			 
		}
		else{
			printf("key %d, initiateur %d\n", value[KEY], value[INITIATOR]);
			//int next=findNextIndex(value[KEY], value[INITIATOR]);
			int next=findNextIndex(value[KEY], rang.CHORD_ID);
			//printf("next %d\n", next);
			
			if ( next !=-1){
				printf("demande de lookup reçu par %d, propagée au finger %d\n", rang.CHORD_ID, chord_finger[next]);		//sleep(1);
				MPI_Send(value, 2, MPI_INT, mpi_finger[next], LOOKUP, MPI_COMM_WORLD); 
			}
			else{
			// ne sert à rien , car le successeur n'a pas l'info des clés gérées
			//MPI_Send(&key, 1, MPI_INT, mpi_finger[0], LASTCHANCE, MPI_COMM_WORLD); 
			
			// rien trouvé, donc aucun noeud responsable
			int fail=-1;
			printf("pas de responsable trouvé dans le réseau\n");
			MPI_Send(&fail, 1, MPI_INT, 0, END, MPI_COMM_WORLD);
			}
		
		}
		return LOOKUP;
	}
	
	/*if(status.MPI_TAG==LASTCHANCE){
		if 	
	
	}*/

}

/* Question 2: Algo de recherche du pair responsable d'une clé CHORD par la paire courante
*/

void lookup(){

	int result_node;
	MPI_Status status;

	initialisation();
	int res=-1;
	while(res!=END){
		res=receive();
	}

}

/******************************************************************************/

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	srand(time(NULL));
	

  	
  	int nb_proc;
  	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

	if (nb_proc != N+1) {
      		printf("Nombre de processus incorrect , il doit être égale à %d , sinon changez la valeur de N dans le programme !\n", N+1);
      		MPI_Finalize();
      		exit(2);
   	}
  	
   	MPI_Comm_rank(MPI_COMM_WORLD, &rang.MPI_ID);
  	
  	range=pow(2,M);
  	
   	if (rang.MPI_ID == 0) {
   		printf("----------------------------------------INFOS DU PROGRAMME----------------------------------------\n");
   		printf("nombre de bits d'une clé M=%d\n", M);
		printf("nombre de pairs dans le réseau pair à pair %d\n", N);
		
  		printf("valeurs des clés dans {0, %d-1}\n", range);
  		printf("--------------------------------------------------------------------------------------------------\n");
      		simulation();
      	
   	} else {
   		
   		
      		lookup();
      		//printf("id du processus courant MPI: %d, CHORD: %d\n", rang.MPI_ID, rang.CHORD_ID);
   	}
   	
   	//test_modulo();
   	//test_app();
  	//test_exists();
  	//printf("test app %d in %d\n", app(7,4, 0), rang.MPI_ID);

	MPI_Finalize(); 
	return 0;
}
