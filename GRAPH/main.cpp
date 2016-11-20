#include "problems.h"

int main (int argc, char *argv[]){

    //N : # of Nodes
    //M : # of Factors
    //p : # of variables that enter in a factor
    //q : # of colors each variable can assume
    
    int N, M, seed;
    int p=3;
    int q=2;

    
    if(argc==4){
        int i = 1;
        N     = atoi(argv[i++]);
        M     = atoi(argv[i++]);
        seed  = atoi(argv[i++]);
    }
    else{
        cout << "argument: N, M, seed" << endl;
    }
    
    
    srand(seed);
    
    //SatInstance G(N,p);
    //f_plantedSatInstance(G,N,M,p,q);

    //XorsatInstance G(N,p);
    //f_XorsatRandomDecimation(G, N, M, p, q);
    
    SatInstance G(N,p);
    f_BPguidedDecimation(G,N,M,p,q);
    
    return 1;
}
