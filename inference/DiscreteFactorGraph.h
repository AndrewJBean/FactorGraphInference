#pragma once

#include "FloatType.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



class Variable
{
public:
    int cardinality;
    vector< int > fns;
};

class Function
{
public:
    vector< int > vars;
    int values;
    int ValuesSize;
};



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



class FactorGraph
{
public:
    int ReadGraph();                                                        // DONE
    int WriteGraph();                                                       // DONE
    static int MyGet(int &);                                                // DONE
    static int MyGet(FloatType &);                                          // DONE
    static void MyPut(int);                                                 // DONE
    static void MyPut(FloatType);                                           // DONE
    static int MyGet(int &,ifstream &);                                     // DONE
    static int MyGet(FloatType &,ifstream &);                               // DONE
    static void MyPut(int,ofstream &);                                      // DONE
    static void MyPut(FloatType,ofstream &);                                // DONE
    int ReadGraphBin();                                                     // DONE
    int WriteGraphBin();                                                    // DONE
    int ReadGraphBin(string const &);                                       // DONE
    int WriteGraphBin(string const &);                                      // DONE
    int ReadGraphBin(ifstream &);                                           // DONE
    int WriteGraphBin(ofstream &);                                          // DONE
    int GetNumVars(){return variables.size();}                              // DONE
    int GetNumFns(){return functions.size();}                               // DONE
    FloatType LogValue(vector< int > const &);

    vector< Variable > variables;
    vector< Function > functions;
    vector< vector< FloatType > > values;
};


