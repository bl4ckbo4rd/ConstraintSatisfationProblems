#include "problems.h"

void f_XorsatGraph(XorsatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

void f_plantedXorsatGraph(XorsatInstance& G, int N, int M, int p, int q){
    
    //here we explicitely define a planted configuration
    //vector<int> ps = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;
    //or we can create it randomly
    vector<int> ps;
    ps.reserve(N);
    for (int i=0; i<N; i++){
        if (2*(double)rand()/RAND_MAX-1 > 0)
            ps.push_back(1);
        else
            ps.push_back(0);
    }
    
    //we build an ErdosRenyi factor graph for which ps is a solution
    G.plantedErdosRenyi(M,ps);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    cout << "print 1 if the planted solution satisfy the formula: " << G.check(ps) << endl;
    
    //we run a BP algorithm starting from the planted solution to verify that it is a solution:
    cout << "running BP on the graph starting from the planted solution to verify that it is a solution: " << endl;
    
    //all the nodes are biased towards the planted solutions
    vector<int> v_bias;
    for (int i=0; i<N; i++){
        v_bias.push_back(i);
    }
    //and these are the color indices to which we set the nodes
    vector<int> v_q = ps;
    
    int verbose=1;
    BP It1(q,G,verbose);
    It1.initDecimation(v_bias,v_q);
    
    It1.BPsweep();
 
    return;
}

void f_leafRemoval(XorsatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    //run the LR algorithm to find the 2-core structure in the verbose mode so to print the values d of each node.
    //d=1 on node i means that i is part of the 2-core structure.
    
    int verbose=1;
    G.LeafRemoval(verbose);
    
    return;
}

void f_XorsatRandomDecimation(XorsatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();

    int verbose=1;
    BP It(q,G,verbose);
    
    //this is the initial number of biased spins
    int N_prime = (int)(0.2*N);
    
    //and these are the color indices to which we set the nodes
    vector<int> v_q;
    v_q.reserve(N_prime);
    for (int i=0; i<N_prime; i++){
        if (2*(double)rand()/RAND_MAX-1 > 0)
            v_q.push_back(1);
        else
            v_q.push_back(0);
    }
    
    //-------------------------------------------------------------------------------------------------
    //here we give information on the structure of the graph
    
    G.graphStructure();
    
    //-------------------------------------------------------------------------------------------------
    //we run the decimation method starting from the biased nodes:
    cout << "running decimation on the graph where the first " << N_prime << " nodes are biased towards: " << endl;
    
    for (int i=0; i<N_prime; i++)
        cout << "bias: " << v_q[i] << endl;
    
    vector<int> v_bias;
    for (int i=0; i<N_prime; i++){
        v_bias.push_back(i);
    }
    
    It.initDecimation(v_bias,v_q);
    It.randomDecimation();
    
}

void f_SatGraph(SatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

void f_plantedSatGraph(SatInstance& G, int N, int M, int p, int q){
    
    //here we explicitely define a planted configuration
    //vector<int> ps = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;
    //or we can create it randomly
    vector<int> ps;
    ps.reserve(N);
    for (int i=0; i<N; i++){
        if (2*(double)rand()/RAND_MAX-1 > 0)
            ps.push_back(1);
        else
            ps.push_back(0);
    }
    
    //we build an ErdosRenyi factor graph for which ps is a solution
    G.plantedErdosRenyi(M,ps);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    cout << "print 1 if the planted solution satisfy the formula: " << G.check(ps) << endl;
    
    //we run a BP algorithm starting from the planted solution to verify that it is a solution:
    cout << "running BP on the graph starting from the planted solution to verify that it is a solution: " << endl;
    
    //all the nodes are biased towards the planted solutions
    vector<int> v_bias;
    for (int i=0; i<N; i++){
        v_bias.push_back(i);
    }
    //and these are the color indices to which we set the nodes
    vector<int> v_q = ps;
    
    int verbose=1;
    BP It1(q,G,verbose);
    It1.initDecimation(v_bias,v_q);
    
    It1.BPsweep();
    
    return;
    
}

void f_BPsweepSat(SatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    int verbose=1;
    BP It(q,G,verbose);
    
    It.BPsweep();
    
}

void f_BPguidedDecimation(SatInstance& G, int N, int M, int p, int q){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    bool verbose=1;
    BP It(q,G,verbose);
    
    vector <int> v_q;
    vector <int> v_bias;
    
    It.initDecimation(v_bias,v_q);
    It.BPguidedDecimation(N);
    
}

