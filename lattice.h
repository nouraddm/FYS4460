#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <string>
#include <../../../Desktop/FYS4460/ProjectOne/include/armadillo>
#include "atom.h"
//#include "generatequantities.h"
#include "lib.h"

using namespace std;
using namespace arma;

class Lattice
{
public:
    Lattice();
    vec getCell();
    void setCell(const vec &inCell);
protected:
    vec cells;
};

#endif // LATTICE_H
