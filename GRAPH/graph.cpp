#include "graph.h"




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Node
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

int Node::numberOfFactors(){
  return v_fac.size();
};




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Factor
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//input variables:
//p_f : index of the factor
//p_p : # of the variables entering in a factor

Factor::Factor(int p_f, int p_p) : f ( p_f ), p(p_p) {
    v_node.reserve(p);
};

int Factor::numberOfNodes(){
    return v_node.size();
};

//--------------------------------------------------------------------------------------------------------------------------------------------// methods of the derived class FactorXorsat

//input variables:
//i,j,k : indeces of the nodes entering in the factor

double FactorXorsat::clause(int i, int j, int k){
    
    /*
    //this is the finite T clause
    int s0   = (int)pow(-1.,i);
    int s1   = (int)pow(-1.,j);
    int s2   = (int)pow(-1.,k);
    double J = (int)pow(-1.,type);
    double beta=2;
    return exp(-beta*(1-J*s0*s1*s2));
    */

    //this is the T=0 clause
    if(type)
        return i^j^k;
    else
        return !(i^j^k);
    
};

//input variables:
//i,j,k : indeces of the nodes entering in the factor

//given a triplet of variables entering in a clause, we can construct only one "planted" clause.
//a clause of the XORSAT problem is specified by the value type being 0 or 1, drawn at randon in the non-planted setting.
//for instance, consider the clause x1 ^ x2 ^ x3.
//it is easy to see that the only choice of type for which clause is satisfied is type = x1 ^ x2 ^ x3.

void FactorXorsat::plantedClause(int i, int j, int k){
    
    //this is the T=0 clause
    type=i^j^k;
    
};

//-----------------------------------------------------------------------------------------------------------------------------------------------// methods of the derived class FactorSat

//input variables:
//i,j,k : indeces of the nodes entering in the factor

double FactorSat::clause(int i, int j, int k){

    bool J1 = v_J[0];
    bool J2 = v_J[1];
    bool J3 = v_J[2];
    
    //this is the T=0 clause
    return (2*J1 - 1) * (i - (1-J1)) | (2*J2 - 1) * (j - (1-J2)) | (2*J3 - 1) * (k - (1-J3));

};

//input variables:
//i,j,k : indeces of the nodes entering in the factor

//given a triplet of variables entering in a clause, we can construct several planted "clauses".
//a clause of the SAT problem is specified by three values J1, J2, J3, being 0 or 1 if the corresponding variables appear negated or unnegated in the clause,
//drawn at random in the non-planted setting.
//for instance, consider the clause x1 | !x2 | !x3 where ! indicates the negation. in this case we set J1=1, J2=0, J3=0.
//it is easy to see that the only choice of J's for which the triplet do not satisfy the clause is 0, 1, 1, and more generally, for a particular realization of (x1,x2,x3) we need
//to exclude J1=!x1, J2=!x2, J3=!x3.

void FactorSat::plantedClause(int i, int j, int k){
    
    vector<bool> temp_v_J, p_J;
    bool J1, J2, J3;
    int flag=1;
    
    while (flag){

        if (2*(double)rand()/RAND_MAX-1 > 0)
            J1=1;
        else
            J1=0;
    
        if (2*(double)rand()/RAND_MAX-1 > 0)
            J2=1;
        else
            J2=0;
    
        if (2*(double)rand()/RAND_MAX-1 > 0)
            J3=1;
        else
            J3=0;

        p_J = make_vector<bool>() << 1-i << 1-j << 1-k;
        temp_v_J = make_vector<bool>() << J1  << J2  << J3;
    
        if (p_J != temp_v_J)
            flag=0;
    }

    //this is the T=0 clause
    v_J = temp_v_J;
    
};




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Graph
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//input variables:
//p_N : # of nodes
//p_p : # of variables entering in a factor

Graph::Graph(int p_N, int p_p) : N ( p_N ), p ( p_p ) {
  v.reserve(N);
  for (int i = 0; i < N; ++i){
    v.push_back (Node(i));
  }
};

int Graph::numberOfTotalFactors(){
  return F.size();
};

int Graph::numberOfTotalNodes(){
  return N;
};

void Graph::factorsOfNode(int i){
  cout << "node " << i << " has " << v[i].numberOfFactors() << " factors: " << endl;
  for (vector<int>::iterator it = v[i].v_fac.begin() ; it != v[i].v_fac.end(); ++it)
    cout << *it << endl;
};

void Graph::nodesOfFactor(int a){
  cout << "factor " << a << " has " << F[a]->numberOfNodes() << " nodes: " << endl;
  for (vector<int>::iterator it = F[a]->v_node.begin() ; it != F[a]->v_node.end(); ++it)
    cout << *it << endl;
};

//------------------------------------------------------------------------------------------------------------------------------------------// methods of the derived class XorsatInstance

//input variables:
//p_a    : factor index;
//p_type : factor type (+1 or 0);
//v_da   : nodes attached to the factor

int XorsatInstance::addFactor(int p_a, bool p_type, vector<int> v_da){
    
    vector <int> v1 = v_da;
    std::sort(v1.begin(), v1.end());
    
    int flag=1;
    
    //before adding a factor we check that the factor does not already exist.
    //if the factor already exists, flag is set to 0.
    //to this aim it is sufficient to check that the first node of v_da does not appear in a factor with the two other nodes of v_da.
    
    for(vector<int>::iterator it_a = v[v_da[0]].v_fac.begin(); it_a != v[v_da[0]].v_fac.end(); it_a++){
        
        vector <int> v2 = F[*it_a]->v_node;
        std::sort(v2.begin(), v2.end());
        
        if (v1 == v2){
            flag=0;
        }
    }
    
    if(flag){
        
        //a is a pointer to the derived class FactorXorsat
        FactorXorsat* a = new FactorXorsat(p_a,p_type,p);
        
        for (int i=0; i<v_da.size(); i++)
            a->v_node.push_back(v_da[i]);
        
        //F, defined as vector of pointer to Factor, can be filled also with pointer to the derived classes such as FactorXorsat
        F.push_back(a);
        
        for(int i=0; i<v_da.size(); i++)
            v[v_da[i]].v_fac.push_back(p_a);
    }
    
    return flag;
    
};

void XorsatInstance::graphStructure(){
    
    cout << endl;
    cout << endl;
    
    cout << "structure of the graph: " << endl;
    cout << endl;
    
    for (int i=0; i<N; i++)
        factorsOfNode(i);
    
    for (int a=0; a<M; a++)
        nodesOfFactor(a);
    
    cout << endl;
    cout << endl;

    cout << "factors types: " << endl;
    for (int a=0; a<M; a++)
        //F is defined as vector of pointer to Factor but we filled it with pointer to the derived class FactorXorsat.
        //In order to use members of the derived class we need to use the cast operator
        cout << a << " type: " << ((FactorXorsat*)F[a])->type << endl;

    cout << endl;
    cout << endl;
    
};

//input variables:
//p_M: # of factors

void XorsatInstance::ErdosRenyi(int p_M){

    M = p_M;
    
    vector<int> v;
    int i,j,k;
    int a = 0;
    bool type;
    int flag;

    while (a<M){
        i = rand() % N ;
        j = rand() % N ;
        k = rand() % N ;

        if (i!=j && i!=k && j!=k){
            v = make_vector<int>() << i << j << k;
            if (2*(double)rand()/RAND_MAX-1 > 0)
                type=1;
            else
                type=0;
        
            flag=addFactor(a, type, v);
            if (flag) a++;
        }
    }

};

//input variables:
//p_M: # of factors
//ps : the planted configuration

void XorsatInstance::plantedErdosRenyi(int p_M, vector<int> &ps){
    
    
    
    M = p_M;
    
    vector<int> v;
    int i,j,k;
    int a = 0;
    int type;
    int flag;
    
    while (a<M){
        i=rand() % N ;
        j=rand() % N ;
        k=rand() % N ;
        
        if (i!=j && i!=k && j!=k){
            v = make_vector<int>() << i << j << k;
            type=0;
            flag=addFactor(a, type, v);
            
            if (flag){
                a++;
                //ptr is the pointer to the last added factor
                //please remind that F is a vector of pointers to Factor's.
                //in order to use this pointer as a pointer to XorsatFactor we need to use the cast operator
                Factor* ptr = F.back();
                ((FactorXorsat*)F[ptr->f])->plantedClause(ps[i],ps[j],ps[k]);
                
                //cout << ps[i] << " " << ps[j] << " " << ps[k] << " " << ((ps[i]+ps[j]+ps[k]) % 2 == ((FactorXorsat*)F[ptr->f])->type) << endl;
            }
        }
    }
};

//input variables:
//ps : the planted configuration

bool XorsatInstance::check(vector<int> &ps){
    
    bool prod=1;
    int i, j, k;
    for (int a=0; a<M; a++){
        i = F[a]->v_node[0];
        j = F[a]->v_node[1];
        k = F[a]->v_node[2];
        
        prod*=((ps[i]+ps[j]+ps[k]) % 2 == ((FactorXorsat*)F[a])->type);

    }
    
    return prod;
}

//input variable:
//verbose : 1/0 for the verbose/non verbose version of the method

void XorsatInstance::LeafRemoval(int verbose){
    
    int count=1, flag, errorflag;
    int l;
    
    //v_d[i]=0 on variables that are removed from the graph in a leaf removal strategy and 1 otherwise
    //F_d[i]=0 on factors that are removed from the graph in a leaf removal strategy and 1 otherwise
    
    vector <int> v_d(N,1);
    vector <int> F_d(M,1);
    
    while (count>0){
        count=0;
        for (int i=0; i<N; i++){
            //here we select a node and we try to eliminate it, with the corresponding constraint to which it is still attached
            if(v_d[i]){
                
                flag=0;
                
                //here we check how many constraints are still attached to the node i
                l=0;
                for (vector<int>::iterator it_a = v[i].v_fac.begin(); it_a != v[i].v_fac.end(); it_a ++){
                    if (F_d[*it_a]==1)
                        l++;
                }
                
                if (l <= 1){
                    flag=1;
                    v_d[i]=0;
                    
                    int errorflag = 0;
                    for(vector<int>::iterator it_a = v[i].v_fac.begin(); it_a != v[i].v_fac.end(); it_a++){
                        if(F_d[*it_a]){
                            F_d[*it_a]=0;
                            errorflag++;
                        }
                    }
                    //this is a check. errorflag has to be 1.
                    if(errorflag > 1)
                        cout << "error in leaf removal" << endl;
                }
                
                count += flag;
                //cout << i << " " << flag << endl;
            }
        }
    }
    
    for(int i=0; i<N; i++)
        v[i].d=v_d[i];
    
    for(int a=0; a<M; a++)
        ((FactorXorsat*)F[a])->d=F_d[a];
    
    if (verbose){
        cout << endl;
        cout << endl;
        cout << "2-core structure (d=1 on frozen variables) : " << endl;
        for (int i=0; i<N; i++)
            cout << "----- d: " << v_d[i] << endl;
        
        cout << endl;
        cout << endl;
    }
    
};

//--------------------------------------------------------------------------------------------------------------------------------------------// methods of the derived class SatInstance

//input variables:
//p_a    : factor index;
//p_v_J  : p=3-component vector whose components are 0 or 1 depending if the corresponding variable is negated or not
//v_da   : nodes attached to the factor

int SatInstance::addFactor(int p_a, vector<bool> p_v_J, vector<int> v_da){
    
    vector <int> v1 = v_da;
    std::sort(v1.begin(), v1.end());
    
    int flag=1;
    
    //before adding a factor we check that the factor does not already exist.
    //if the factor already exists, flag is set to 0.
    //to this aim it is sufficient to check that the first node of v_da does not appear in a factor with the two other nodes of v_da.
    
    for(vector<int>::iterator it_a = v[v_da[0]].v_fac.begin(); it_a != v[v_da[0]].v_fac.end(); it_a++){
        
        vector <int> v2 = F[*it_a]->v_node;
        std::sort(v2.begin(), v2.end());
        
        if (v1 == v2){
            flag=0;
        }
    }
    
    if(flag){
        
        //a is a pointer to the derived class FactorSat
        FactorSat* a = new FactorSat(p_a,p_v_J,p);
        
        for (int i=0; i<v_da.size(); i++)
            a->v_node.push_back(v_da[i]);
        
        //F, defined as vector of pointer to Factor, can be filled also with pointer to the derived classes such as FactorSat
        F.push_back(a);
        
        for(int i=0; i<v_da.size(); i++)
            v[v_da[i]].v_fac.push_back(p_a);
    }
    
    return flag;
    
};

void SatInstance::graphStructure(){
    
    cout << endl;
    cout << endl;
    
    cout << "structure of the graph: " << endl;
    cout << endl;
    
    for (int i=0; i<N; i++)
        factorsOfNode(i);
    
    for (int a=0; a<M; a++)
        nodesOfFactor(a);
    
    cout << endl;
    cout << endl;
    
    cout << "factors types: " << endl;
    for (int a=0; a<M; a++)
        //F is defined as vector of pointer to Factor but we filled it with pointer to the derived class FactorSat.
        //In order to use members of the derived class we need to use the cast operator
        cout << a << " J's: " << ((FactorSat*)F[a])->v_J[0] << " " << ((FactorSat*)F[a])->v_J[1] << " " << ((FactorSat*)F[a])->v_J[2] << endl;
    
    cout << endl;
    cout << endl;
    
};

//input variables:
//p_M: # of factors

void SatInstance::ErdosRenyi(int p_M){
    
    M = p_M;
    
    vector<int> v;
    vector <bool> v_J;
    int i,j,k;
    bool J1, J2, J3;
    int a = 0;
    int flag;
    
    while (a<M){
        i = rand() % N ;
        j = rand() % N ;
        k = rand() % N ;
        
        if (i!=j && i!=k && j!=k){
            v = make_vector<int>() << i << j << k;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J1=1;
            else
                J1=0;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J2=1;
            else
                J2=0;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J3=1;
            else
                J3=0;
            
            v_J = make_vector<bool>() << J1 << J2 << J3;
            
            flag=addFactor(a, v_J, v);
            if (flag) a++;
        }
    }
    
};

//input variables:
//p_M : # of factors
//ps  : planted configuration

void SatInstance::plantedErdosRenyi(int p_M, vector<int> &ps){
    
    M = p_M;
    
    //ps is the planted solution;
    
    vector<int> v;
    vector <bool> v_J;
    int i,j,k;
    int a = 0;
    bool J1, J2, J3;
    int flag;
    
    while (a<M){
        i=rand() % N ;
        j=rand() % N ;
        k=rand() % N ;
        
        if (i!=j && i!=k && j!=k){
            v   = make_vector<int>() << i  << j  << k;
            v_J = make_vector<bool>() << 0 << 0 << 0;
            flag=addFactor(a, v_J, v);
            
            if (flag){
                a++;
                //ptr is the pointer to the last added factor
                //please remind that F is a vector of pointers to Factor's.
                //in order to use this pointer as a pointer to XorsatFactor we need to use the cast operator
                Factor* ptr = F.back();
                ((FactorSat*)F[ptr->f])->plantedClause(ps[i],ps[j],ps[k]);
                
                //J1 = ((FactorSat*)F[ptr->f])->v_J[0];
                //J2 = ((FactorSat*)F[ptr->f])->v_J[1];
                //J3 = ((FactorSat*)F[ptr->f])->v_J[2];
                //cout << ps[i] << " " << ps[j] << " " << ps[k]  << " " << ((2*J1 - 1) * (ps[i] - (1-J1)) | (2*J2 - 1) * (ps[j] - (1-J2)) | (2*J3 - 1) * (ps[k] - (1-J3))) << endl;
               
            }
        }
    }
};

//input variables:
//ps  : planted configuration

bool SatInstance::check(vector<int> &ps){
    
    bool prod=1;
    int i, j, k;
    bool J1, J2, J3;
    for (int a=0; a<M; a++){
        i = F[a]->v_node[0];
        j = F[a]->v_node[1];
        k = F[a]->v_node[2];
        
        J1 = ((FactorSat*)F[a])->v_J[0];
        J2 = ((FactorSat*)F[a])->v_J[1];
        J3 = ((FactorSat*)F[a])->v_J[2];
        
        prod*=((2*J1 - 1) * (ps[i] - (1-J1)) | (2*J2 - 1) * (ps[j] - (1-J2)) | (2*J3 - 1) * (ps[k] - (1-J3)));
        
    }
    
    return prod;
}




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Messages
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


Messages::Messages(int p_q, Graph& p_G) : q (p_q), G (p_G) {
    
    N=G.numberOfTotalNodes();
    M=G.numberOfTotalFactors();

    eta_NodeToFac.resize(N);
    nu_FacToNode.resize(M);
    marginal.resize(N);
    bias.resize(N);
    
    initMessages();
    
};


void Messages::initMessages(){

    //Init messages eta to uniform distributions
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        int size_di=G.v[it_i->n].numberOfFactors();  //di is the set of factors attached to node i. size_di is the number of such factors
        eta_NodeToFac[it_i->n].resize(size_di);
      
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
        
            vector <double> eta(q,1./q);
            eta_NodeToFac[it_i->n][index_a].resize(q);
            eta_NodeToFac[it_i->n][index_a]=eta;

        
        }
    }

    //init messages nu to 0.
    //we don't need to do that because we usually start running BP equations from the eta messages (and thus we only need to initialize them)
    //and the nu are computed afterwords. by the way we need to store memory for these vectors.
    for(vector<Factor*>::iterator it_a=G.F.begin(); it_a !=G.F.end(); ++it_a){
        int size_da=G.F[(*it_a)->f]->numberOfNodes();  //da is the set of nodes attached to factor a. size_da is the number of such nodes and should be p
        nu_FacToNode[(*it_a)->f].resize(size_da);
        for (vector<int>::iterator it_i = G.F[(*it_a)->f]->v_node.begin(); it_i != G.F[(*it_a)->f]->v_node.end(); ++it_i){
            int index_i = distance (G.F[(*it_a)->f]->v_node.begin(), it_i);
            vector <double> nu(q,0.);
            nu_FacToNode[(*it_a)->f][index_i].resize(q);
            nu_FacToNode[(*it_a)->f][index_i]=nu;
        }
    }
    
    //set no bias, i.e set bias to the uniform distribution
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        bias[it_i->n].resize(q);
        vector <double> b(q,1./q);
        bias[it_i->n]=b;
    }

    //init marginals to uniform distributions
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        marginal[it_i->n].resize(q);
        marginal[it_i->n]=bias[it_i->n];
    }

};

void Messages::nuUpdate (){

    for (vector<Factor*>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
        for (vector<int>::iterator it_i = G.F[(*it_b)->f]->v_node.begin() ; it_i != G.F[(*it_b)->f]->v_node.end(); ++it_i){
            //cout << "message from factor " << it_b->f << " to node " << *it_i << endl;
            //compute message nu_b_to_i
            int index_i = distance (G.F[(*it_b)->f]->v_node.begin(), it_i);
            vector <vector <double> > in_eta;
            for(vector<int>::iterator it_j = G.F[(*it_b)->f]->v_node.begin(); it_j != G.F[(*it_b)->f]->v_node.end(); ++it_j){
                if(*it_j!=*it_i){
                    //cout << "j: " << *it_j << endl;
                    vector<int>::iterator it = find(G.v[*it_j].v_fac.begin(), G.v[*it_j].v_fac.end(), (*it_b)->f);
                    int index_b = distance (G.v[*it_j].v_fac.begin(), it);
                    in_eta.push_back(eta_NodeToFac[*it_j][index_b]);
                }
            }
            
            vector <double> nu(q,0.);
            //this part needs to be generalized to general kind of clauses.
            //now it is tailed to describe a 3-body interaction.
            for (int k0=0; k0<q; ++k0){
                for (int k1=0; k1<q; ++k1){
                    for (int k2=0; k2<q; ++k2){
                        nu[k0] += G.F[(*it_b)->f]->clause(k0,k1,k2) * in_eta[0][k1] * in_eta[1][k2];
                    }
                }
            }
            
            //message nu_b_to_i
            nu_FacToNode[(*it_b)->f][index_i] = nu;
        }
    }

};

void Messages::etaUpdate (){

    for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            //cout << "message from node " << it_i->n << " to factor " << *it_a << endl;
            //compute message eta_i_to_a
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
            vector <double> fin_eta(q,1.);
            for(vector<int>::iterator it_b = G.v[it_i->n].v_fac.begin(); it_b != G.v[it_i->n].v_fac.end(); ++it_b){
                if(*it_b!=*it_a){
                    //cout << "b: " << *it_b << endl;
                    vector<int>::iterator it = find(G.F[*it_b]->v_node.begin(), G.F[*it_b]->v_node.end(), it_i->n);
                    int index_i = distance (G.F[*it_b]->v_node.begin(), it);
                    for(int k=0; k<q; ++k){
                        fin_eta[k] *= nu_FacToNode[*it_b][index_i][k];
                    }
                }
            }
            for(int k=0; k<q; ++k)
                fin_eta[k]*=bias[it_i->n][k];
            
            //message eta_i_to_a
            eta_NodeToFac[it_i->n][index_a] = fin_eta;
        }
    }

};

int Messages::etaNormalize(){

    int errorflag=1;

    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
            double sum=0.;
            for (int k=0; k<q; k++){
                sum+=eta_NodeToFac[it_i->n][index_a][k];
            }
            if(sum!=0)
                for (int k=0; k<q; k++){
                    eta_NodeToFac[it_i->n][index_a][k]/=sum;
                }
            else
                errorflag=0;
            
        }
    }
    return errorflag;
};

void Messages::nodeMarginals(){
    //for each node i
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        vector <double> marg(q,1.);
        //we compute the product of the messages nu_FacToNode[a][i] coming from the factor attached to i
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            int index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            for (int k=0; k<q; k++){
                marg[k]*=nu_FacToNode[*it_a][index_i][k];
            }
        }
        //and we take into account the bias on the node i
        for (int k=0; k<q; k++){
            marg[k]*=bias[it_i->n][k];
        }
        
        //finally we normalize the marginal
        double sum=0.;
        for (int k=0; k<q; k++){
            sum+=marg[k];
        }
        
        if (!sum)
            cout << "node " << it_i->n << " receives conflicting messages, see below: " << endl;
        else{
            for (int k=0; k<q; k++){
                marg[k]/=sum;
            }
            marginal[it_i->n]=marg;
        }
        
    }

};

void Messages::setHardBias(vector<int>& v_bias, vector<int>& v_q){
    
    //v_bias contains the indices of nodes that have to be biased
    //v_q contains the color towards which they have to be biased
    
    int size=v_bias.size();
    if (size != v_q.size()){
        cout << "error: v_q and v_bias have to have the same size!" << endl;
        return;
    }
    else{
        for (int i=0; i<size; i++){
            int n=v_bias[i];
            vector <double> b(q,0);
            //we set the bias towards the index color specified in v_q[i]
            b[v_q[i]]=1;
            bias[n]=b;
        }
    }
    
};

//the three following methods of the class Messages print the state of the BP algorithm
void Messages::etaState(){

  cout << endl;
  cout << "---*---*---*---printing messages from nodes to factors---*---*---*---" << endl;

  for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
    for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
      cout << "message from node " << it_i->n << " to factor " << *it_a << endl;
      int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
      cout << eta_NodeToFac[it_i->n][index_a][0] << " " << eta_NodeToFac[it_i->n][index_a][1] << endl;
    }
  }

};

void Messages::nuState(){

  cout << endl;
  cout << "---*---*---*---printing messages from factors to nodes---*---*---*---" << endl;

  for (vector<Factor*>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
    for (vector<int>::iterator it_i = G.F[(*it_b)->f]->v_node.begin() ; it_i != G.F[(*it_b)->f]->v_node.end(); ++it_i){
      cout << "message from factor " << (*it_b)->f << " to node " << *it_i << endl;
      int index_i = distance (G.F[(*it_b)->f]->v_node.begin(), it_i);
      cout << nu_FacToNode[(*it_b)->f][index_i][0] << " " << nu_FacToNode[(*it_b)->f][index_i][1] << endl;
    }
  }

};

void Messages::marginalState(){

  cout << endl;
  cout << "---*---*---*---printing marginals---*---*---*---" << endl;

  for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
    cout << "marginal of node " << it_i->n << endl;
    cout << marginal[it_i->n][0] << " " << marginal[it_i->n][1] << endl;
  }

};




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class BP
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


BP::BP(int p_q, Graph& p_G) : q (p_q), G (p_G), mess (q,G) { N=G.N; };

void BP::initDecimation(vector<int>& v_bias, vector<int>& v_q){
    
    for(int i=0; i<N; i++)
        notFixedSpins.push_back(i);
    
    if(v_bias.size()){

        int size=v_bias.size();
        if (size != v_q.size()){
            cout << "error: v_q and v_bias have to have the same size!" << endl;
            return;
        }
        else{
            for (int i=0; i<size; i++){
                int n=v_bias[i];
                fixedSpins.push_back(n);
                fixedValues.push_back(v_q[i]);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
                
                //we bias eta's towards the index color specified in v_q[i]
 
                //for (vector<int>::iterator it_a = G.v[n].v_fac.begin(); it_a != G.v[n].v_fac.end(); ++it_a){
                //    int index_a = distance (G.v[n].v_fac.begin(), it_a);
                //    vector <double> eta(q,0);
                //    eta[v_q[i]]=1;
                //    mess.eta_NodeToFac[n][index_a]=eta;
                //}
                
            }
        }

        mess.setHardBias(v_bias, v_q);

    }

};

int BP::BPsweep(int verbose){

    mess.nuUpdate();
    mess.etaUpdate();
    int flag = mess.etaNormalize();
    if(verbose){
        mess.nodeMarginals();
        BPprint();
    }
    return flag;
};

void BP::BPprint(){

    mess.nuState();
    mess.etaState();
    mess.marginalState();

};

void BP::BPiteration(int T, int verbose=0){
    
    for (int t=0; t<T; t++)
        int flag = BPsweep(verbose);
    
}

void BP::randomDecimation(int verbose=0){

    //Pictorially this is how this method works:
    
    //time t-1:
    //set bias and fix a group G of variables
    //time t:
    //update nu  (using eta at time t-1)
    //update eta (using nu  at time t and bias at time t-1)
    //at this time step, no variables get fixed because of the set G, but some other variables
    //may get fixed because of the variables that we fixed at time t-2.
    //time t+1:
    //update nu  (using eta at time t)
    //use these nu's to set bias and fix other variables
    
    //it is clear then that the variables fixed at time t-1 have some effect on the others after 2 iterations, not 1.
    
    //g is the number of spins that became frozen (decimated) at each iteration.
    //g_past is the number of spins that became frozen (decimated) at the previous iteration.
    //we do keep track of both for the reason explained above.
    //f is a boolean variable that is set to 1 when frustrated variables are found
    
    int  g, g_past=1;
    int  t=0;
    bool f=0;
    
    while (f==0 && fixedSpins.size()<N){
        cout << "--------------------------------------time: " << t << "------------------------------------ " << endl;
        cout << endl;
        if (verbose){
            cout << "frozen variables:" << " (size=" << fixedSpins.size()    << ")" << endl;
            vec_print(fixedSpins);
            cout << "free variables:"   << " (size=" << notFixedSpins.size() << ")" << endl;
            vec_print(notFixedSpins);
        }
        
        int flag = BPsweep(verbose);

        if(!flag)
            break;

        f=findFrustratedSpins();
        g=fixSpins(verbose);
        
        
        //if no contraddiction are found and no spins gets fixed, we randomly pick one spin and we fix it
        if (f == 0 && (g+g_past==0)){
            fixRandomSpin(verbose);
        }
        
        
        t++;
        g_past = g;
    }
    
    if (verbose){
        cout << endl;
        cout << "----------------------------------END OF ITERATION----------------------------------------" << endl;
        cout << "(maybe partial) solution: " << endl;
        cout << "size of the solution: "     << fixedSpins.size() << endl;
        vec_print(fixedSpins);
        vec_print(fixedValues);
    }
    
    cout << endl;

};

void BP::BPguidedDecimation(int T, int verbose=0.){
    
    reservePreviousMarginals();
    
    double eps=1.;
    
    BPsweep(verbose);
    storePreviousMarginals();
    
    while (eps>pow(10.,-3)){
        BPsweep(verbose);
        storePreviousMarginals();
        eps = compareMarginals();
    }
    
};

void BP::reservePreviousMarginals(){
    
    prev_marginal.resize(N);
    
    for(int i=0; i<N; i++)
        prev_marginal[i].resize(q-1);
    
};

void BP::storePreviousMarginals(){
    
    for(int i=0; i<N; i++){
        for (int k=0; k<q-1; k++)
            prev_marginal[i][k] = mess.marginal[i][k];
    }
    
};

double BP::compareMarginals(){
    
    double tmp, max = 0.;
    
    for(int i=0; i<N; i++){
        for (int k=0; k<q-1; k++)
            tmp = abs(prev_marginal[i][k] - mess.marginal[i][k]);
            if (tmp > max)
                max = tmp;
    }
    
    return max;
    
};
    
void BP::findMostBiased(){
    
    int i_max=0;
    double tmp, max=0.;
    
    for (int i=0; i<N; i++){
        tmp = abs(mess.marginal[i][0]-mess.marginal[i][1]);
        if (tmp > max){
            i_max=i;
            max=tmp;
        }
    }
    
    return;
};




    
void BP::fixRandomSpin(int verbose){
    
    vector<int> v_bias;
    vector<int> v_q;
    bool col;
    
    int N_free = notFixedSpins.size();
    int i = rand() % N_free ;
    if (2*(double)rand()/RAND_MAX-1 > 0)
        col = 1;
    else
        col = 0;
    
    int n = notFixedSpins[i];
    
    v_bias.push_back(n);
    v_q.push_back(col);
    fixedSpins.push_back(n);
    fixedValues.push_back(col);
    notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
    
    if(v_bias.size() != 0)
        mess.setHardBias(v_bias,v_q);
    
    if (verbose){
        cout << "*************** printing the nodes that gets randomly frozen: " << endl;
        for(int i=0; i<v_bias.size(); i++)
            cout << v_bias[i] << " " << v_q[i] << endl;
    }
    
};


int BP::findFrustratedSpins(){
    
    vector<int>().swap(frustratedSpins);
    
    bool f=0;
    
    //for each node i
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        vector <double> marg(q,1.);
        //we compute the product of the messages nu_FacToNode[a][i] coming from the factor attached to i
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            int index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            for (int k=0; k<q; k++){
                marg[k]*=mess.nu_FacToNode[*it_a][index_i][k];
            }
        }
        //and we take into account the bias on the node i
        for (int k=0; k<q; k++){
            marg[k]*=mess.bias[it_i->n][k];
        }
        
        //finally we normalize the marginal
        double sum=0.;
        for (int k=0; k<q; k++){
            sum+=marg[k];
        }
        
        if (!sum){
            frustratedSpins.push_back(it_i->n);
            f=1;
        }
    }
    
    return f;
    
};

int BP::fixSpins(int verbose){

    int q=1;

    //v_bias contains the indices of nodes that have to be biased
    //v_q contains the color towards which they have to be biased

    vector<int> v_bias;
    vector<int> v_q;


    //for each factor a
    for (vector<Factor*>::iterator it_a = G.F.begin() ; it_a != G.F.end(); ++it_a){
        //we look at all the nodes to which a is sending messages
        for (vector<int>::iterator it_i = G.F[(*it_a)->f]->v_node.begin() ; it_i != G.F[(*it_a)->f]->v_node.end(); ++it_i){
            //if this not is not fustrated (i.e. it is not receinving conflicting messages from factors), we fix it
            bool flag = find(frustratedSpins.begin(), frustratedSpins.end(), *it_i) != frustratedSpins.end();
            if(!flag){
                
                int index_i = distance (G.F[(*it_a)->f]->v_node.begin(), it_i);
                //if this message contains a clear indication of the color towards which
                //the node should be biased, we add the node to v_bias and the color to v_q,
                //after having checked that the node has not been frozen yet.
                if (mess.nu_FacToNode[(*it_a)->f][index_i][q]==1.){
                    bool contains = find(fixedSpins.begin(), fixedSpins.end(), *it_i) != fixedSpins.end();
                    if(!contains){
                        v_bias.push_back(*it_i);
                        v_q.push_back(1);
                        fixedSpins.push_back(*it_i);
                        fixedValues.push_back(1);
                        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), *it_i), notFixedSpins.end());

                    }
                    continue;
                }
                if (mess.nu_FacToNode[(*it_a)->f][index_i][q]==0.){
                    bool contains = find(fixedSpins.begin(), fixedSpins.end(), *it_i) != fixedSpins.end();
                    if(!contains){
                        v_bias.push_back(*it_i);
                        v_q.push_back(0);
                        fixedSpins.push_back(*it_i);
                        fixedValues.push_back(0);
                        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), *it_i), notFixedSpins.end());
                    }
                    continue;
                }
            }
        }
    }
    
    if(v_bias.size() != 0)
        mess.setHardBias(v_bias,v_q);
    

    if (verbose){
        cout << "*************** printing the nodes that gets frozen at this time step: " << endl;
        for(int i=0; i<v_bias.size(); i++)
            cout << v_bias[i] << " " << v_q[i] << endl;

        cout << endl;
        cout << endl;
    }
    
    return v_bias.size();

};


void vec_print(vector<int>& vec){
    for (int i=0; i<vec.size(); i++)
        cout << vec[i] << ' ';
    cout << endl;
}

