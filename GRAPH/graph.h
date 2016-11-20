#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>

using namespace std;



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a node
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Node{
public:

    Node(int p_n) : n ( p_n ) {};
  
    int n;                                                  //n is the index of the node
    vector <int> v_fac;                                     //v_fac contains the indices of factors attached to n
    
    int d;                                                  //this value is determined by the Leaf Removal Algorithm.
                                                            //it is 1 if the variable belongs to the 2-core structure, 0 otherwise
    
    int numberOfFactors();                                  //this method returns the number of factors to which the node n is attached
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Factor{
public:
    
    Factor(int, int);

    int f;                                                  //f is the index of the factor
    int p;                                                  //p is the number of variables entering in a factor. p=3 in the 3-XORSAT problem

    vector <int> v_node;                                    //v_node contains the indices of nodes attached to f
    
    int numberOfNodes();                                    //this method returns the numer of variables to which a factor is attached. this number is p by definition
    
    virtual double clause(int, int, int) { return 0; };     //this method specifies the clause represented by a factor. it is specified in the specific
                                                            //factors implmented in the derived classes.
    
};

//------------------------------------------------------------------------------------------------------------------------------------------// derived class representing a XORSAT clause

class FactorXorsat: public Factor{
public:
    
    FactorXorsat(int p_f, bool p_type, int p_p) : Factor(p_f,p_p) { type = p_type; };
    
    bool type;                                              //type specifies the kind the XORSAT factor we are considering.
                                                            //it can be 1 or 0.
                                                            //in the first case we consider the XOR of p variables.
                                                            //in the second case we consider the negation of the XOR of p variables.
    
    int d;                                                  //this value is determined by the Leaf Removal Algorithm. It is 1 if factor belongs to the reduced graph, 0 otherwise.
    
    double clause(int, int, int);                           //this method implements the XORSAT clause represented by a factor.
    
    void plantedClause(int, int, int);                      //this method looks at the input variables and sets the type of the factor accordingly.
                                                            //input: the p=3 variables involved in each clause
    
};

//----------------------------------------------------------------------------------------------------------------------------------------------// derived class representing a SAT clause

class FactorSat: public Factor{
public:
    
    FactorSat(int p_f, vector<bool> p_v_J, int p_p) : Factor(p_f,p_p) {v_J.resize(p); v_J = p_v_J; };
    
    vector <bool> v_J;                                      //v_J is a set of p binary variables specifying the SAT clause.
                                                            //consider the clause x1 | !x2 | !x3 where ! indicates the negation.
                                                            //be Jk the k-component of the vector v_J.
                                                            //in this case we set J1=1, J2=0, J3=0 and we write the clause as
                                                            //(2 J1 - 1) [x1 - (1 - J1)] | (2 J2 - 1) [x2 - (1 - J2)] | (2 J3 - 1) [x3 - (1 - J3)]
                                                            //in other words, we just send xk to wk = (2 Jk - 1) [xk - (1 - Jk)]
                                                            //and we see that if Jk is 1, then wk = xk while if Jk=0, then wk = 1 - xk = ! xk

    double clause(int, int, int);                           //this method implements the SAT clause represented by a factor. input: the p=3 variables involved in each clause
    
    void plantedClause(int, int, int);                      //this method looks at the input variables and sets the vector v_J accordingly.
                                                            //input: the p=3 variables involved in each clause

};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor graph
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Graph{
public:
    
    Graph(int, int);
    
    int N;
    int p;
    int M;
    
    vector <Node> v;                                        //v contains all the nodes of the graph
                                                            //this vector is filled by the constructor.
    
    vector <Factor*> F;                                     //F contains pointers to all the factors of the graph.
                                                            //this vector is filled in the classes derived from Graph, that specify the model under consideration.
                                                            //since the model is specified by the kind of factor, this vector will be filled by pointers to the classes
                                                            //derived from Factor.
                                                            //thus this vector, that should contains pointers to Factors, will eventually contain pointers to derived Factors.
                                                            //this is the reason why pointers are used. Without using them, we would have had "object slicing".

    
    int numberOfTotalFactors();                             //this method returns the total number of factors in the graph
    int numberOfTotalNodes();                               //this method returns the total number or nodes in the graph

    void factorsOfNode(int);                                //this method returns the factors attached to the input node
    void nodesOfFactor(int);                                //this method returns the nodes attached to the input factor
    
    //the following methods are virtual
    virtual void ErdosRenyi(int) {};                        //this method build an ER graph
    virtual void plantedErdosRenyi(int, vector<int> &) {};

    virtual void graphStructure() {};                       //this method prints the structure of the graph
    

    
};

//-------------------------------------------------------------------------------------------------------------------------------// derived class representing a XORSAT instance (formula)

class XorsatInstance : public Graph{
public:
    
    XorsatInstance(int p_N, int p_p) : Graph(p_N,p_p) {};
    
    int addFactor(int, bool, vector<int>);                  //the addFactor method implement the operations needed when adding a factor in the graph
                                                            //namely one needs to create a FactorXorsat object a, store its neighour variables and for each of them
                                                            //add the index of a as a neghbouring factor.
                                                            //input variables: factor index, factor type (+1 or 0), vector of nodes attached to the factor.
                                                            //output: is 1 if the operation could be done, 0 if not (the factor already exists).
    
    void ErdosRenyi(int);                                   //ErdosRenyi method generates an ER graph.
                                                            //input: number of links
    
    void plantedErdosRenyi(int, vector<int> &);             //this method generates a planted ER graph, starting from an initial configuration
                                                            //input: number of constraints, input configuration
    
    
    void LeafRemoval(int);                                  //this function implement the LR strategy to identity the 2-core structure of the initial formula
                                                            //if a node appears only in one factor, it is removed, together with this factor.
                                                            //the elimination of this factor, may make nodes that appeared in two factors, appear only in one.
                                                            //they will be thus eliminated when encountered.
                                                            //we iterate on all the nodes of the graph that have not been removed yet. Let's call this process a graph sweep.
                                                            //after each graph sweep we look at how many new nodes have been removed.
                                                            //we stop when there are no more nodes eliminate.
    
    void graphStructure();                                  //this method prints the structure of the graph
    
    bool check(vector<int> &);                              //this method checks that the planted solution verifies all the clauses.
                                                            //output: 1 if it does, 0 if it does not.
    
};

//----------------------------------------------------------------------------------------------------------------------------------// derived class representing a SAT instance (formula)

class SatInstance : public Graph{
public:
    
    SatInstance(int p_N, int p_p) : Graph(p_N,p_p) {};
    
    int addFactor(int, vector<bool>, vector<int>);          //the addFactor method implement the operations needed when adding a factor in the graph
                                                            //namely one needs to create a FactorSat object a, store its neighour variables and for each of them
                                                            //add the index of a as a neghbouring factor.
                                                            //input variables: factor index, p=3-component vector whose components are 0 or 1 if the corresponding variable
                                                            //appear negated or unnegated in the clause, vector of nodes attached to the factor.
                                                            //output: is 1 if the operation could be done, 0 if not (the factor already exists).
    
    
    void ErdosRenyi(int);                                   //ErdosRenyi method generates an ER graph.
                                                            //input: number of links
    
    void plantedErdosRenyi(int, vector<int> &);             //this method generates a planted ER graph, starting from an initial configuration
                                                            //input: number of constraints, input configuration
    
    void graphStructure();                                  //this method prints the structure of the graph
    
    bool check(vector<int> &);                              //this method checks that the planted solution verifies all the clauses.
                                                            //output: 1 if it does, 0 if it does not.
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------------------------------------// class reprensenting the BP messages
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Messages{
public:
    
    Messages(int, Graph &);

    int N;
    int M;
    int q;
    Graph G;

    vector < vector < vector <double> > > eta_NodeToFac;        //eta_NodeToFac contains the messages from nodes to factors
    vector < vector < vector <double> > > nu_FacToNode;         //nu_FacToNode contains the messages from factors to node
    //Loosely speaking
    //eta_NodeToFac[i][a][k] is the k-th component of the message from node i to factor a
    //while
    //nu_FacToNode[b][i][k]  is the k-th component of the message from node b to factor i

    vector < vector <double> > marginal;                        //marginal probability for each spin
    vector < vector <double> > bias;                            //bias contains the biased distribution for each spin. it is uniform for unbiased spins
  
    void initMessages();                                        //this method initializes eta's to umbiased values 1/q and the nu's to 0. bias are set
                                                                //to the uniform distribution and the marginals are set to the uniform distribution as well.
                                                                //it is invoked inside the constructor

    //nuUpdate and etaUpdate are the two functions that implement a BP sweep on the whole graph:
    void nuUpdate();                                            //this method updates the nu messages
    void etaUpdate();                                           //this method updates the eta messages

    int  etaNormalize();                                        //this method normalizes the eta's. it returns 0 if at least one node receives conflicting messages, 1 otherwise
    void nodeMarginals();                                       //this method computes (do not print) node marginals

    void setHardBias(vector<int> & , vector<int> &);            //this method defines hard biases towards color q for each node contained
                                                                //in the first input vector and according to the color index contained in the second input vector.
                                                                //default input vectors are empty vectors.
                                                                //input variables: vector of nodes to be biased, vector of colors to towards which node biases have to be biased
    
    
    //these are the printing functions
    void etaState();                                            //this method returns the state of the eta's
    void nuState();                                             //this method returns the state of the nu's
    void marginalState();                                       //this method returns the state of the marginals

};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------// class reprensenting BP algorithms
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class BP{
public:
    
    BP(int, Graph &);
    
    int N;
    int q;
    Graph G;
    Messages mess;
    
    void initDecimation(vector<int> &, vector<int> &);          //this method set an hard bias on some specified variables
                                                                //the first vector contains the biased variables
                                                                //the second vector contais the colours toward which the biased nodes are biased.
    
    vector <int> fixedSpins;                                    //this vector contains the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <int> fixedValues;                                   //this vector contains the values of the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the constructor and by the method fixSpins during the decimation process.
    
    vector <int> notFixedSpins;                                 //this vector contains the indices of the spin that are not frozen.
                                                                //at t=0, it is formed by all the spins. As t increases and the decimation process continues, its size
                                                                //decreases.
                                                                //it is filled by the constructor and by the method fixSpins during the decimation process.
    
    vector <int> frustratedSpins;                               //this vector contains the spins that receive conflicting messages from the factors.
                                                                //it is filled by the method findFrustratedSpins().
    
    vector < vector <double> > prev_marginal;
    
    
    int findFrustratedSpins();                                  //this method evaluates, for each node, all the messages from its factors (like nodeMarginals) and seeks if they contain
                                                                //contraddictory information. In this case, this spin is frustrated.
                                                                //when frustrated spins are found the method returns 1.
    
    int fixSpins(int);                                          //this method fix spins values, by calling setHardBias, for spins that receive a clear indication
                                                                //towards one color from at least one of its factors.
                                                                //it is invoked in the method BP_decimation,
                                                                //when a node receives a clear message from at least one if its factors.
                                                                //it also fills the vector fixedSpins and, when doing this, erase nodes from the vector notFixedSpins.
                                                                //input variable:
                                                                //verbose: set it to 1 to print the nodes that get frozen
    
    void fixRandomSpin(int);                                    //this method fix one spin at random at a random value, by calling setHardBias.
                                                                //it is invoked in the method BP_decimation, no spins get frozen because of messages.
                                                                //it also fills the vector fixedSpins and, when doing this, erase nodes from the vector notFixedSpins.
                                                                //input variable:
                                                                //verbose: set it to 1 to print the node that get frozen

    int BPsweep(int);                                           //this method updates all the nu's and all the eta's in the graph.
                                                                //it returns 0 if at least one node receives conflicting messages, 1 otherwise
                                                                //input variable:
                                                                //verbose: set it to 1 to print the messages, 0 otherwise
    
    void randomDecimation(int);                                 //this method propagates warnings at each time steps and, at each time step,
                                                                //it fixes spins by calling fixSpins. when this is not possible it calls fixRandomSpins.
                                                                //it can stop earlier if contraddicting messages are founded.
                                                                //input variable:
                                                                //verbose: set it to 1 to print the messages, 0 otherwise


    void BPguidedDecimation(int, int);

    void reservePreviousMarginals();

    void storePreviousMarginals();
    
    double compareMarginals();

    void findMostBiased();
    
    void BPiteration(int, int);                                 //this method iterates BP equations by calling BP_sweep
                                                                //input variables:
                                                                //T: time of BP_sweep;
                                                                //verbose: set it to 1 to print the messages, 0 otherwise
    
    void BPprint();                                             //this method prints the BP messages and the marginals

};


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------// useful functions

//this function defines allows to fill a vector in one single line
//it has been downloaded from https://gist.github.com/pablomtz/5577626
//examples:
//vector<int> v = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;

template <typename T>
class make_vector {
 public:
  typedef make_vector<T> my_type;
  my_type& operator<< (const T& val) {
    data_.push_back(val);
    return *this;
  }
  operator std::vector<T>() const {
    return data_;
  }
 private:
  std::vector<T> data_;
};

//this function print all the elements of a vector.
void vec_print(vector<int>& vec);


