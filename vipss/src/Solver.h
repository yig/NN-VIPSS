#ifndef SOLVER_H
#define SOLVER_H


#include <iostream>
#include <vector>
#include <armadillo>
#include <nlopt.hpp>


// using namespace std;

enum ConstraintMethod{
    CM_EQUAL,
    CM_LESS,
    CM_GREATER,
    CM_RANGE
};


struct triple{
    int i;
    int j;
    double val;
    triple(int i,int j,double val):i(i),j(j),val(val){}
};

struct LinearVec{

    std::vector<int>ts;
    std::vector<double>ws;
    double a;
    double b;
    int label;
    ConstraintMethod mt;


    LinearVec(){}
    LinearVec(int label):label(label){}
    void clear();
    void push_back(int t, double w);
    void set_label(int label);
    void set_ab(double a, double b);
};

struct QuadraticVec{

    std::vector<triple>ts;
    double b;
    ConstraintMethod mt;

    QuadraticVec(){}
    QuadraticVec(std::vector<triple>&ts, double b):ts(ts),b(b){}

    void push_back(int i,int j,double val);
    void set_b(double b);
    void clear();
};

struct NLOptConstraint{
    nlopt::vfunc constraintfunc;
    ConstraintMethod mt;
    void *data;

};


struct Solution_Struct{
    int Statue;
    double init_energy;
    double energy;
    double time;
    std::vector<double>solveval;
    arma::vec solvec_arma;

    Solution_Struct():init_energy(-1),energy(-1){}
    void init(int n);
};

class Solver{


public:

//    GRBEnv objenv;
//    GRBModel objmodel;

//    std::vector<GRBVar>c_vars;

public:
//    Solver():objenv(GRBEnv()),objmodel(GRBModel(objenv)){

//        objmodel.getEnv().set(GRB_IntParam_OutputFlag, 0);

//    }

    Solver(){}


public:


//    static int solveQuadraticProgramming(std::vector<triple> &M, std::vector<triple> &Ab,int n, std::vector<double>&solveval);

//    static int solveQuadraticProgramming(std::vector<triple> &M, std::vector<LinearVec> &Ab,int n, std::vector<double>&solveval);
//    static int solveQuadraticProgramming(arma::mat &M, std::vector<LinearVec> &Ab,int n, Solution_Struct &sol);


//    static int solveLinearLS(std::vector<triple> &M, std::vector<double> &b, int n, std::vector<double>&solveval);

//    static int solveQCQP(arma::mat &M, std::vector<LinearVec> &Ab, std::vector<QuadraticVec> &Cb, int n, Solution_Struct &sol);

    //int solveQP_multiupdata_main(Solution_Struct &sol, std::vector<LinearVec> &Ab, int n, bool suppressinfo = true);
    //int solveQP_multiupdata_init(arma::mat &M, int n);



    static int nloptwrapper(std::vector<double>&lowerbound,
                   std::vector<double>&upperbound,
                   nlopt::vfunc optfunc,
                   void *funcPara,
                   std::vector<NLOptConstraint>&nl_constraints,
                   double tor,
                   int maxIter,
                   Solution_Struct &sol
                   );

    static int nloptwrapper(std::vector<double>&lowerbound,
                   std::vector<double>&upperbound,
                   nlopt::vfunc optfunc,
                   void *funcPara,
                   double tor,
                   int maxIter,
                   Solution_Struct &sol
                   );

    static int nloptwrapperDirect(std::vector<double>&lowerbound,
                   std::vector<double>&upperbound,
                   nlopt::vfunc optfunc,
                   void *funcPara,
                   double tor,
                   int maxIter,
                   Solution_Struct &sol
                   );

};




















#endif // SOLVER_H
