/*
AUTEURS: Inas KACI 3873337
	Amandine Ta 3300320		

*/
/*
Commandes d'éxecution et de compilation
mpicc -o exo2 exercice2.c -lm
mpirun -np 9  --oversubscribe ./exo2

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

#define N 8
#define M 4

#define NON_INITIATEUR 0
#define INITIATEUR 1
#define BATTU 2
#define ELU 3


// tags de messages
#define TAGINIT 0
#define  TAGSTATE 1 //permet de signaler à un processus s'il est initiateur ou pas 
#define OUT 2
#define IN 3
#define END 4
/**
hypthèses:
	- Anneau bidirectionnel
	- chaque noeud possède:
		* un identifiant CHORD et MPI unique
		* une référence vers succ[i] son successeur et pred[i] son prédecesseur
		
	- nombre non nul d'initiateurs ( plusieurs candidats simultanées possibles )
	- Aucun noeud ne connait la taille de l'anneau
	- communications fiables et asynchrones 
	- Exécution sans fautes
	- compexité <O(N²)
	- Reseau statique en cours de l'execution de l'algorithme ( tant qu'il n'a pas finit de s'executer, un noeud ne peux pas rejoindre ou quitter le réseau )
	
	En partant de ces hypothèses, l'algorithme qui les satisfait est l'algorithme d'élection 
	étudié en TD pour la configurantion Bidirectionnelle asynchrone: *Hirschberg & Sinclair*
description de l'algorithme
	
	Dans un reseau pair à pair, pour prendre toute décision, il faut d'abord élire un leader afin d'éviter toute incohérence, d'ou le fait que l'on va se base sur l'algorithme Hirschberg & Sinclair. Nous allons le modifier afin que le noeud leader contienne à la fin la liste de tous les noeuds puis qu'il la diffuse aux reste des noeuds pour que chacun puisse construire sa propre table finger. Ceci revient donc à utiliser un algorithme à vague pour rassembler des informations du réseau puis les diffuser
	Etant donné les valeurs d'identifiants CHORD sont prises dans l'intervalle [0, 2^M -1] et qu'elles sont uniques, notre tableau de liste de noeuds résultat est au plus de taille 2^M
	Ainsi, le message envoyé par l'initiateur contient un tableau contenant en case 0 son identifiant.
	 
	 À l'étape k, le site envoie dans l'anneau dans les deux directions un jeton out avec un tabelau de taille 2^M contenant son identité à la case 0
	 Ces jetons parcourent une distance 2^k avant de revenir à l'emetteur avec in comme type associé
	 à la réception des deux jetons in, l'initiateur passe à l'étape suivante (k+1)
	 
	 Un site j voit passer un jeton out :
	 si tab[INITIATEUR]>j ou j n'est pas initiateur, j transmet le jeton à son voisin en rajoutant son noeud à la prochaine case vide ( qui contient -1 ) dans le tableau 
	 si tab[INITIATEUR]<j et j est initiateur, le jeton est détruit( non transmis )
	 si tab[INITIATEUR]==j l'initiateur tab[INITIATEUR]  contient toutes les noeuds du reseau car il est le noeud leader et a fait le tour complet de l'anneau
	 
	 Il transmets donc ce tableau complet des identifiants CHORD à son voisin avec le TAG END ,
	 le voisin propage le tableau complet à son voisin et ainsi de suite jusqu'à faire un tour complet
	 Ainsi tous les noeuds connaissent le contenu du tableau complet donc connaissent tous tous les identifiants présents sur le réseua
	 Donc chaque noeud est capable de calculer sa propre table finger selon l'algorithme classique:
	 si p=CHOTD_ID du processus courant
	 table de taille M et pour tout i dans {0, M-1} : finger[i]= (p+2^i)[2^M]
	
	 	
	
justification de sa correction

	Terminaison:
	L'algorithme Hirschberg & Sinclair se termine garantit l'élection d'un leader unique,
	Une fois ce leader élu, il envoi un unique message un à son successeur qui le diffuse récursivement à son successeur jusqu'à revenir sur le noeud leader. Par hypthèse, Etant donné que le nombre de noeud du réseau est fixe en cours d'execution de cet algorithme, il y aura N envoi de message avant de retomber sur le leader, donc cette partie de l'algorithme est en O(N)
	
	Décision:
	Une décision est prise par tous les processus ( tous créent leurs tables finger ). Propriété de vivacité assurée 
	
	Sureté: L'algorithme Hirschberg & Sinclair garantit qu'il y a un unique leader élu, 
	donc, c'est bien le même tableau qui est renvoyé à tous les noeuds 
	
	Dépendance:
	L'élection de leader avec L'algorithme Hirschberg & Sinclair garantie qu'il n'y a d'élection tant que tous les noeuds n'ont pas été parcourus, de même la diffusion du tableau se fait en parcourant tous les noeuds. Puis, le calcul de la table des fingers ( la décision ) se fait sur chaque noeud. On a donc bien la décision qui est précédée causalement par un événement de chaque processus
	
justification de sa complexité
	Algorithme de Hirschberg & Sinclair est en O(NlogN) car il y au plus N initateurs et pour chaque initiateur au plus logN rounds avec l'envoi de deux messages chacun
	L'étape finale de diffusion du tableau résultat est en O(N)
**/

/**
Code 
*/

/*************************************************************************
Variables
**************************************************************************/

struct Identity{
	int MPI_ID;
	int CHORD_ID;
}
;
/*chaque pair dispose de :*/
struct Identity rang; // son identifiant CHORD et MPI
int range;
int state;
int mpi_finger[M];
int chord_finger[M];

/*************************************************************************
Tool functions
**************************************************************************/
bool exists(int k, int tab[], int NB ){
	//bool res= false;
	int i;
	for (i=0; i<NB; i++){
		if(tab[i]==k){
			return true;
		}
	
	}
	return false;
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

/* retourne la valeur positive associée au modulo*/
int modulo (int val, int Z){
	int mod=val%Z;
	if(mod<0){
		return mod+Z;
	}	
	return mod;
}

void simulation(){
	// tirer aléatoirement les id CHORD
	int i;
	int identifiers[N+1];
	identifiers[0]=-1;
	// le processus 0 ne participe pas à l'overlay de la DHT
	// les valeurs de sa table finger prennent -1
	for (i=0; i<M; i++){
		mpi_finger[i]=-1;
		chord_finger[i]=-1;
	}
	
	for(i=1; i<=N; i++){
		identifiers[i]=hash_func(i, identifiers, N+1);
		// envoyer la valeur de l'id CHORD au paire associé
		MPI_Send(&identifiers[i], 1, MPI_INT, i, TAGINIT, MPI_COMM_WORLD); 
	} //O(N)
	
	//afficher la liste des valeurs 
	printf(" chord identifiers : {");
	for(i=0; i<=N; i++){
		printf("%d , ", identifiers[i]);
		
	}
	printf(" }\n");
	

	
	// construction d'un anneau bidirectinnel selon l'ordre des des rangs MPI 
	// cela revient juste à autoriser à n noeud d'envoyer un message à son successeur et son précedecesseur
	// c'est à dire si le pair d'identifiant MPI p, il peut communiquer pendant la contruction de la DHT avec les noeuds p-1 et p+1
	 
	 

	 int cpt=0;
	 int states[N+1];
	 
	
	
	 for( i=1; i<=N;i++){
	 	states[i]=rand()%2;
	 	if (states[i]==INITIATEUR){
	 		cpt++;
	 	}
	 	
	  }
	 if (cpt==0){
	 	// aucun processus n'a été choisi en tant qu'initiateur
		// il faut avoir au moins un initiateur
		// en choisir un aléatoirement
		int initiator = 1+(rand()%N);
		states[initiator]=INITIATEUR;
	 }
	 for(i=1; i<=N; i++){
	 	MPI_Send(&states[i], 1, MPI_INT, i, TAGSTATE, MPI_COMM_WORLD); 
	 }
	 
	 // chaque finger affiche sa table de finger et se termine 
	 
	 
	 

}

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
	// N n'est pas connu dans le programme courant par les paires
	// le calculer
	int NB=0;
	for(i=1; i<=range; i++){
		if( identifiers[i]!=-1){
			NB++;
		}
	}
	for (i=1; i<=NB; i++){
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

//int round=0;
int nb_in;

#define TTL 0
int K=0;
void init_round(int k){
	if(k==0){
	
	 }
	 nb_in=0;
	 // envoi jeton à ses deux voisins 
	 /* identifiers[0] contient le l'identifiant de l'initiateur
	 identifiers[1] contient 2^k
	 pas la peine d'envoyer la valeur de l'emetteur, on sait que en identifiants MPI c'est soit i-1, soit i+1
	 à partir de identifers[2] nous avons les autres noeuds du réseau
	 
	 */
	 int identifiers[range+1];
	 identifiers[INITIATEUR]=rang.CHORD_ID;
	 identifiers[TTL]=pow(2,k);
	 int i;
	 for(i=2; i<range+2; i++){
	 
	 	identifiers[i]=-1;
	 }
	 int next=modulo(rang.MPI_ID+1, N+1);
	 if ( next==0){
	 	// incrémentation car le processus 0 fait office de simulateur et ne fait pas partie du réseau 
	 	
	 	next++;
	 }
	 
	MPI_Send(identifiers, range+1, MPI_INT, next, OUT, MPI_COMM_WORLD); 
	//printf("OUT: from %d to %d, Initiateur %d\n", rang.CHORD_ID, next, identifiers[INITIATEUR]);
	int previous=modulo(rang.MPI_ID-1, N+1);
	if ( previous==0){
	 	// incrémentation car le processus 0 fait office de simulateur et ne fait pas partie du réseau 
	 	
	 	previous = N;
	 }
	MPI_Send(identifiers, range+1, MPI_INT, previous, OUT, MPI_COMM_WORLD);
	//printf("OUT: from %d to %d, Initiateur %d\n", rang.CHORD_ID, previous, identifiers[INITIATEUR]);
	K=k+1;
	 
}

/**
écrire value à la première case contenant -1
**/
void addId(int tab[], int NB, int value){
	int i;
	for (i=2; i<NB; i++){
		if(tab[i]==-1){
			tab[i]=value;
			return;
		}
	}

}
/**
	returns true if val is left to ref in the ring 
*/
bool isLeft(int val, int ref){
	int diff1= modulo(ref-val, N+1);
	int diff2= modulo(val-ref, N+1);
	
	if(diff1<diff2){
		return true;
	}
	return false;
}

bool test_isLeft(){
	printf("expected 1, res: %d\n", isLeft(5, 1));
	printf("expected 1, res: %d\n", isLeft(1, 2));
	printf("expected 0, res: %d\n", isLeft(3, 1));
	printf("expected 0, res: %d\n", isLeft(2, 5));
}
int receive(){
	int identifiers[range+1];
	MPI_Status status;
	MPI_Recv(identifiers, range+1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (status.MPI_TAG==OUT){
		if(state!=INITIATEUR || rang.CHORD_ID< identifiers[INITIATEUR]){
			state=BATTU;
			int k=identifiers[TTL];
			addId(identifiers, range+1, rang.CHORD_ID);
			//printf(" add %d starting from %d\n", rang.CHORD_ID,identifiers[INITIATEUR] );
			
			//printf("intermediate identifiers at INITIATOR %d { ", identifiers[INITIATEUR] );
			int j;
			for(j=0; j<=range; j++){
			
				printf(" %d, ", identifiers[j]);
			}
			printf("\n"); 
			
			k--;
			if( k==0){
				identifiers[TTL]=k;
				MPI_Send(identifiers, range+1, MPI_INT,status.MPI_SOURCE, IN, MPI_COMM_WORLD);
				//printf("IN: from %d to %d, Initiateur %d\n", rang.CHORD_ID,status.MPI_SOURCE, identifiers[INITIATEUR]);
			}
			else{
				int next;
				if(isLeft(status.MPI_SOURCE,rang.MPI_ID)){
					next=modulo(rang.MPI_ID+1, N+1);
					if ( next==0){
	 	
	 				next++;
					 }
				}
				else{
					next=modulo(rang.MPI_ID-1, N+1);
					if ( next==0){
	 	
	 					next = N;
					 }
				}
				identifiers[TTL]=k;
				MPI_Send(identifiers, range+1, MPI_INT,next, OUT, MPI_COMM_WORLD);
				//printf("OUT: from %d to %d, Initiateur %d\n", rang.CHORD_ID, next, identifiers[INITIATEUR]);
			}
		}
		if( rang.CHORD_ID==identifiers[INITIATEUR]){
			// le message a fait le tour, donc l'initiateur est leader et le tableau contient tous les noeuds .
			state=ELU;
			// envoyer le tableau complet à tous les noeud en commençant par le noeud voisin qui lui va le diffuser à son tour
			int next;
			next=modulo(rang.MPI_ID+1, N+1);
			if ( next==0){
	 			next++;
			}
			MPI_Send(identifiers, range+1, MPI_INT,next, END, MPI_COMM_WORLD);
			//printf("END: from %d to %d, Initiateur %d\n", rang.CHORD_ID, next, identifiers[INITIATEUR]);
			/*printf("result identifiers at INITIATOR %d { ", rang.CHORD_ID);
			int j;
			for(j=0; j<=range; j++){
			
				printf(" %d, ", identifiers[j]);
			}
			printf("\n"); */
			
		}
	}
	if (status.MPI_TAG==IN){
		/*printf("IN sur %d avec identifiers of initiator %d {",rang.CHORD_ID, identifiers[INITIATEUR ]);
		
		int j;
		for(j=0; j<=range; j++){
		
			printf(" %d, ", identifiers[j]);
		}
		printf("\n"); */
		
		
		if( identifiers[INITIATEUR]!=rang.CHORD_ID){
			int next;
				if(isLeft(status.MPI_SOURCE,rang.MPI_ID)){
					next=modulo(rang.MPI_ID+1, N+1);
					if ( next==0){
	 	
	 				next++;
					 }
				}
				else{
					next=modulo(rang.MPI_ID-1, N+1);
					if ( next==0){
	 	
	 					next = N;
					 }
				}
				identifiers[TTL]=0;
				MPI_Send(identifiers, range+1, MPI_INT,next, IN, MPI_COMM_WORLD);
				//printf("IN: from %d to %d, Initiateur %d\n", rang.CHORD_ID, next, identifiers[INITIATEUR]);
		}
		else{
			nb_in++;
			if(nb_in==2){
				printf("----------------------------------------- STARTING ROUND %d -----------------------------------------\n" , K);
				init_round(K);
			}
		}
	}
	if( status.MPI_TAG==END){
		// calculer la finger table
		
		//printf("-------%d-------- \n", identifiers[i]);
		int j;
		struct Identity tmp;
		for(j=0; j<M; j++){
			tmp=getFinger(rang.CHORD_ID, j, identifiers);
			
			mpi_finger[j]=tmp.MPI_ID;
			chord_finger[j]=tmp.CHORD_ID;
			
		}
		// l'afficher
		printf("finger %d { " ,rang.CHORD_ID);
		for(j=0; j<M; j++){
			printf("%d, ", chord_finger[ j]);
		}
		printf("}\n");
		
		// transmettre le message si pas leaer 
		if( state!=ELU){
			int next=modulo(rang.MPI_ID+1, N+1);
			if ( next==0){
	 			next++;
			}
			MPI_Send(identifiers, range+1, MPI_INT,next, END, MPI_COMM_WORLD);
			//printf("END: from %d to %d, Initiateur %d\n", rang.CHORD_ID, next, identifiers[INITIATEUR]);
			
		}
	}
		
	return status.MPI_TAG;
}
void initialisation(){
	//chaque proc récupère son numéro d'ID CHORD
	MPI_Status status;
	// recevoir son identifiant CHORD du simlateur d'id 0
	MPI_Recv(&rang.CHORD_ID, 1, MPI_INT, 0, TAGINIT, MPI_COMM_WORLD, &status);
	// réception de l'état du procressus courant
	MPI_Recv(&state, 1, MPI_INT, 0, TAGSTATE, MPI_COMM_WORLD, &status);
	if( state==INITIATEUR){
		printf("le pair chord: %d, MPI: %d est initiateur\n", rang.CHORD_ID, rang.MPI_ID);
		printf("----------------------------------------- STARTING ROUND %d -----------------------------------------\n" , K);
		init_round(0);
		
	}
	
	
}



void compute_finger_table(){
	initialisation();
	int res=receive();
	while (res!=END){
		res=receive();
	}
	//receive();	
	

}


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	srand(time(NULL));
	
	int nb_proc;
  	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

	if (nb_proc != N+1) {
      		printf("Nombre de processus incorrect !, il doit être égale soit à %d ou être changé dans le programme \n", N+1);
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
   	
   		compute_finger_table();
   	
   	}
   	//test_isLeft();
	
	MPI_Finalize(); 
	return 0;
}
