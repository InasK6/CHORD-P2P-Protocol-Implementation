/*
AUTEUR: Inas KACI 3873337
	Amandine Ta 3300320
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

#define M 3 // nomble de bits d'un identifiant de pair
#define N 2 // nombre de paire définit dans le réseau
#define KEY 0
#define INITIATOR 1

/* Tags des différent message */
#define TAGINIT 0 // tag du message d'initialisation
#define LOOKUP 1  // tag des messages de recherche du noeud responsable d'une clé
#define END 2	// tag signalant la fin de la recherche d'un noeud
#define ADDID 3

/*
Propositon de l'algorithme :
*/
/*
Hypothèse de base : 
- Communications fiables et asychrones
- Exécution sans fautes
- Aucun noeud ne connait la taille de l'anneau
- Le nouveau noeud est capable initialement de n’envoyer des messages qu’à un unique pair de la DHT ( choisi arbitrairement ). Pour envoyer un message à d’autres pairs, il devra être informé de leur existence par un pair déjà présent dans la DHT. Donc il peut communiquer avec n'importe quel pair à condition de savoir qu'il existe
- DHT CHORD bien initialisée
- toute paire a une liste inverse_p qui donne la liste des identifiants qui l'utilisent.
- On connaît l'identifiant que l'on veut rentrer ( ajouter ) dans la DHT.
- On sait qu'il n'est pas encore utilisé ( déjà présent ) dans la DHT donc pas besoin de vérifier cette propriété.

- chaque noeud possède:
		* un identifiant CHORD et MPI unique
		* une référence vers succ[i] son successeur 
		* une table de finger CHORD et MPI
		* Une table inverse qui contient les identifiants dont la table finger contient le noeud courant
		
Etat global que l'on cherche à garantir:
- Nouveau noeud bien inséré
- Les fingers tables de tous les processus du réseau sont mis à jour en fonction de cette insertion, donc cohérentes avec le nouvel état de la DHT


1- Construction des tables inverse
- Le nouveau noeud contacte le noeud de la DHT choisi aléatoirement
- Ce dernier fait un lookup pour chercher le noeud responsable de la clé du nouveau noeud ( son successeur ) et renvoyer ses identifiants au nouveau noeud
  Ceci se fait en O(log N) 
- Le nouveau noeud envoie un message REQUEST pour demander les tables inverse et finger associées de son successeur
- Le nouveau noeud traverse la table reverse pour trouver le prédecesseur de son successeur qui deviendra son prédecesseur à lui 
- Le nouveau noeud envoie un message REQUEST à son predecesseur pour obtenir sa table finger
- Le nouveau noeud construit sa finger table en fonction des deux finger tables de son successeur et son prédecesseur et envoie O(log N) messages UPDATE vers les noeuds contenus dans sa finger tables afin qu'ils mettent à jour leurs tables inverse O(N)
- Le nouveau noeud envoie un message UPDATE à chaque noeud de la table inverse de son successeur afin de signaler sa création et lancer une mise à jour de la table finger sur chacun des noeuds destination. Pire cas: O(N-1) et récurisvement mettre à jour les tables inverse des noeuds concernés ( le nouveau noeud et le successeur ) O(N) 


On a donc une compléxité O(logN + N+ N) =O(N) ceci étant le pire cas, nous pouvons considérer que la compléxité de l'algorithme en moyenne est < N , donc sous linéaire 


Pour résumer:

---- 1er étape ----

Trouvéer la clé i à scindér
Puis insérer la valeur x tel que 
succ_x = succ_i;
succ_i = x;
initialisation des finger de x;

---- 2ème étape ----

Parcourire la liste des inverses de i :
 - envoyer un message a chaque inverse de i pour leur signaler la création de la clé x.
 - chaque processus p qui reçoit le message doit :
    - parcourir sa liste de finger et vérifier s'il y a i :
        - S'il y a i, verifier si la clé doit être mise à jour à x ou garder la valeur i :
            - Si mise à jour: modifier la valeur dans la table puis rajouter la valeur p dans la liste des inverses de x.
        - A la fin, Si la table de finger ne contient plus aucun i, enlever p de la liste des inverses de i.

---- Fin algo ----

Justification de la correction de l'algorithme:

	Terminaison:
	L'algorithme de recherche du successeur se termine car le nombre de noeuds lors de la recherche est fixe
	Les modifications suite à l'ajout d'un noeud concernent un nombre fini de noeuds ( successeurs, prédecesseur, noeuds dans finger et inverse de successeurs et prédecesseur ). Ces modifications étant constantes, elles se terminent aussi.
	 Le nombre de messages envoyés est fini, donc l'algorithme d'ajout de noeud se termine 
	
	
	Vivacité:
	Propriété de vivacité assurée, un noeud finit par trouver l'endroit ou il doit s'insérer en un temps fini
	
	Sureté: Les noeuds ayant des identifiants unique, on est assuré qu'il n'y aura qu'un seul noeud successeur correspondant 
	donc, c'est bien le nouveau noeud est bien ajouté au bon endroit
	Chaque noeud amenée à subir modification est au courant de l'arrivée de la nouvelle paire
Et chaque processus qui a reçu une modification en informe ceux pour qui les modification apportées induisent une mise à jour. Ceci assure donc la cohérence de l'état global après l'insertion du noeud
	
	Dépendance:
	L'algorithme de recherche garantie qu'il n'y a d'ajout de noeud tant que tous les noeuds n'ont pas été parcourus. On a donc bien la décision qui est précédée causalement par un événement de chaque processus
	


**Exemple du déroulement de l'algo dans une DHT**

Inicialisation de la DHT CHORD :

nombre de bits d'une clé M=3
nombre de pairs dans le réseau pair à pair 4
valeurs des clés dans {0, 7}

pair présent : {2,4,5,6}
id: 2, succ_2 : 4, finger_2: { 4, 4, 6 } inverse_2: {5, 4, 6}
id: 4, succ_4 : 5, finger_4: { 5, 6, 2 } inverse_4: {2, 5}
id: 5, succ_5 : 6, finger_5: { 6, 2, 2 } inverse_5: {4}
id: 6, succ_6 : 2, finger_6: { 2, 2, 2 } inverse_6: {2, 4, 5}

On veut rajoutée le pair d'id : 0
création de la paire :
id: 0, succ_0 : 2, finger_0: {?} inverse_0:{?}

Calcule du paire qui a pour succeseur le future succeseur de 0
-> 0 envoie un message à 6 New_SUCC 0 de la DHT (msg : 1)

La paire 6 reçoit New_SUCC 0 :
- Met à jour succ_6 = 0
- Vérifie si elle à des 2 dans finger_6 -> oui
	finger_6[0] -> 0
	finger_6[1] -> 0
- On a rajoutée des 0 dans finger_6 donc envoie un message fin à 0 pour rajoutée 6 à la liste des inverse_0 (msg:2)
id: 6, succ_6 : 0, finger_6: { 0, 0, 2 } inverse_6: {2, 4, 5}
-> 6 envoie un messge à 2 INSER 0 (msg : 3)

La paire 0 reçoit La terminaison de 6:
-> 0 met a jour sont finger_0 et inverse_0;

La paire 2 reçoit INSER 0 de 6:
-> 2 envoie un message à 4 et 5 CHECK 0 (msg:4&5)
-> 2 envoie un message de terminaison à 0 (msg:6)

La paire 0 reçoit La terminaison de 2:
-> 0 met a jour sont finger_0 et inverse_0;

La paire 4 reçoit CHECK 0 de 2:
vérifie pour tout les finger_4 qui sont égale à 2 si il change:
	Ne change pas.
-> 4 envoie un message de terminaison à 0 (msg:7)

La paire 0 reçoit La terminaison de 4:
-> 0 met a jour sont finger_0 et inverse_0;

La paire 5 reçoit CHECK 0 de 2:
vérifie pour tout les finger_5 qui sont égale à 2 si il change:
	finger_5[1] -> 0
- On a rajoutée des 0 dans finger_5 donc envoie un message fin à 0 pour rajoutée 5 à la liste des inverse_0 (msg:8)
La paire 0 reçoit La terminaison de 5:
-> 0 met a jour sont finger_0 et inverse_0;

Paires Final : {0,2,4,5,6}
id: 0, succ_0 : 2, finger_0: { 2, 2, 4 } inverse_0:{5,6}
id: 2, succ_2 : 4, finger_2: { 4, 4, 6 } inverse_2: {5, 4, 6}
id: 4, succ_4 : 5, finger_4: { 5, 6, 2 } inverse_4: {2, 5}
id: 5, succ_5 : 6, finger_5: { 6, 1, 2 } inverse_5: {4}
id: 6, succ_6 : 2, finger_6: { 0, 0, 2 } inverse_6: {2, 4, 5}
*/


/*
Code: 

- initialisation avec ajout de inverse 
- piocher du nouvelle paire à insérer 
- Recherche du successeur de la nouvelle paire
il manque l'insertion
*/
/*

Commandes d'exécution et compilation:
mpicc -o exo3 exercice3.c -lm
mpirun -np 3  --oversubscribe ./exo3
*/
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
// Ajout exercice 3
int inverse_mpi[8]; // liste des identifiant mpi des inverse
/*
 * Fonctions Outils
 */

/*
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

int hash_func(int id, int tab[], int NB){
	// tirage d'une valeur aléatoire

	int res=rand()%range;
	// boucle utile pour garantir l'unicité dans l'ensemble de clés associées
	while(exists(res, tab, NB)){
		res=rand()%range;
	}
	return res;
}
/* returns true if k \in ]a,b] */
bool app(int k,int a,int b){

	if(a<b)
		return k>a && k<=b;

	else
		return (k>a && k<range ) || (k>=0 && k<=b);

}

/* retourne la valeur positive associée au modulo */
int modulo (int val, int Z){
	int mod=val%Z;
	if(mod<0){
		return mod+Z;
	}
	return mod;
}

void affichage(int identifiers[], int successors[] ){
	printf("------------------------------------\n");
	int i;
	printf("identifiants CHORD :{");
	for( i=0; i<N+1; i++) {
		if(N==i) printf("%d ", identifiers[i]);
		else printf("%d , ", identifiers[i]);
	}
	printf("}\n");

	printf("successeurs CHORD :{");
	for( i=0; i<N+1; i++) {
		if(N == i ) printf("%d ", successors[i]);
		else  printf("%d, ", successors[i]);
	}
	printf("}\n");
}

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
	for(i=1; i<=N; i++){
		identifiers[i]=hash_func(i, identifiers, N+1);
		// envoyer la valeur de l'id CHORD au paire associé
		MPI_Send(&identifiers[i], 1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
	}
	// noeud initiant la simulation n'a pas de successeur, ne fait pas partie du réseau de paires

	successors[0]=-1;
	int fingers_MPI[N+1][M];
	int fingers_CHORD[N+1][M];
	int p = pow(2,M);
	int inverse_MPI[N+1][p];
	int j=0;
	MPI_Status status;

	//inicialisation de inverse_MPI
	for(int num_l = 0;num_l<N+1;num_l++){
		for(int num_c = 0;num_c<p;num_c++){
			inverse_MPI[num_l][num_c] = 0;
		}
	}

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

	for(int num_p=1; num_p<=N; num_p++){
		int id_chord = identifiers[num_p];
		for(int num_c=0; num_c<M; num_c++){
			int id_proc = fingers_MPI[num_p][num_c];
			inverse_MPI[id_proc][id_chord]=id_chord;
		}
	}

	affichage(identifiers, successors);
	// affichage des tables des fingers
	for(i=0; i<=N; i++){
		printf("finger %d { " , identifiers[i]);
		for(j=0; j<M; j++){
			printf("%d, ", fingers_CHORD[i][j]);
		}
		printf("}\n");

	}
	//printf("i = 3,2, %d\n",fingers_CHORD[3][2] );


	printf("Matrice inverse : \n");
	for(int num_l = 0;num_l<N+1;num_l++){
		printf("inverse  %d : [ ",identifiers[num_l] );
		for(int num_c = 0;num_c<p;num_c++){
			if(inverse_MPI[num_l][num_c] !=0)
				printf("%d ,",inverse_MPI[num_l][num_c] );
		}
		printf("] \n");
	}
	printf("--------------------------------------------------------------------------------------------------\n");

	// envoyer les fingers associés pour chaque noeud
	for (i=1; i<=N; i++){
		MPI_Send(fingers_MPI[i], M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
		// j'utilise un synchronous send pour eviter d'avoir à renvoyer un message par
		// processus qui signale que son initialisation a bien été faite
		MPI_Send(fingers_CHORD[i], M, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
		MPI_Send(inverse_MPI[i],8, MPI_INT, i, TAGINIT, MPI_COMM_WORLD);
	}

	// après la bonne initialisation
	// Récéption d'un message de tous les processus signalant qu'ils ont bien été initialisés
	int done=0;
	for(i=1; i<=N; i++){

		MPI_Recv(&done, 1, MPI_INT,MPI_ANY_SOURCE, TAGINIT, MPI_COMM_WORLD, &status);
	}

	/* Initiate add paire*/
	int new_id = hash_func(i, identifiers, N+1);

	printf("**************************************************************************\n");
	printf("Déclanchement de l'ajout d'un nouveau pair d'id : %d \n", new_id);
	printf("**************************************************************************\n");

	// tirer de manière aléatoire un identifiant de paire
	// tirage selon l'identifiant MPI, donc pas besoin de vérifier l'existence du pair
	// étant donné qu'il y en a N et que leurs identifiants MPI se succédent
	int random=(rand()%N)+1;
	//identifiant CHORD associé:
	int chord_random=identifiers[random];
	int key=new_id;
	// contenu du message
	int msg[2];
	msg[KEY]=key;
	msg[INITIATOR]=chord_random;
	// envoi de la demande de recherche au pair
	printf("**************************************************************************\n");
	printf("Recherche à quel paire est associer la future paire %d à partir du noeud chord %d\n", key, chord_random);
	printf("**************************************************************************\n");
	MPI_Send(msg, 2, MPI_INT, random, LOOKUP, MPI_COMM_WORLD);


	//attente de response avec le résultat de la recherche
	int result_node;

	MPI_Recv(&result_node, 1, MPI_INT,MPI_ANY_SOURCE, END, MPI_COMM_WORLD, &status);
	printf("le resultat %d reçu par le simulateur\n", result_node);
	//MPI_Send(&new_id,1, MPI_INT, result_node,ADDID,MPI_COMM_WORLD);
	//propagation du message de terminaison à tous les processus
	for(i=1; i<=N; i++){
		MPI_Send(&result_node, 1, MPI_INT,i, END, MPI_COMM_WORLD);
	}

}



/*
	fonctions exécutées par les Paire
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
	// Rajout de la liste des inverse
	MPI_Recv(inverse_mpi, 8, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
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
			int next=findNextIndex(value[KEY], rang.CHORD_ID);
			//printf("next %d\n", next);

			if ( next !=-1){
				printf("demande de lookup reçu par %d, propagée au finger %d\n", rang.CHORD_ID, chord_finger[next]);
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
}


void lookup(){

	int result_node;
	MPI_Status status;

	int res=-1;
	while(res!=END){
		res=receive();
	}

}


void insertion (){
	initialisation();
	lookup();
	//int succ_id = findNexIndex(id, value[INITIATOR]);


}

/******************************************************************************/

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	srand(time(NULL));

  int nb_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

	if (nb_proc != N+1) {
      		printf("Nombre de processus incorrect !\n");
      		MPI_Finalize();
      		exit(2);
   	}

  MPI_Comm_rank(MPI_COMM_WORLD, &rang.MPI_ID);

  range=pow(2,M);

  if (rang.MPI_ID == 0) {
  // Inicialisation de da DHT CHORD
  	printf("---------------INFOS DU PROGRAMME-----------------\n");
  	printf("nombre de bits d'une clé M=%d\n", M);
		printf("nombre de pairs dans le réseau pair à pair %d\n", N);
  	printf("valeurs des clés dans {0, %d-1}\n", range);
  	printf("--------------------------------------------------\n");
    simulation();

	}else{
		insertion();
	}
	MPI_Finalize();
	return 0;
}
