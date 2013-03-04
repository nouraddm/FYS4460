#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <iostream>
#include <string>
#include <../../../Desktop/FYS4460/ProjectOne/include/armadillo>
#include "atom.h"
#include "generatequantities.h"
#include "lib.h"

using namespace std;
using namespace arma;

class Potentials
{
public:
    Potentials();
    vec Zero_potential(const vec3 &dr);
    vec Lennard_Jones_potential(const vec3 &dr);
protected:
    vec3 f;
    vec3 zeroVector;
    double innerProduct;
};

#endif // POTENTIALS_H
