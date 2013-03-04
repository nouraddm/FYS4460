#ifndef GENERATEQUANTITIES_H
#define GENERATEQUANTITIES_H

#include <iostream>
#include <string>
#include <../../../Desktop/FYS4460/ProjectOne/include/armadillo>
#include "atom.h"
#include "cell.h"
#include "lib.h"

class Potentials;

using namespace std;
using namespace arma;

class GenerateQuantities
{
public:
    GenerateQuantities(string nameOfThermostat);
    double sampleNormal(long *idum);
    void printToFile();
    void printVelocity();
    void generatePosition();
    void generateVelocity();
    void generateForce();
    void genTimeDevelepment();
    void generateCells();
    void acceleration(Atom *atom1, Atom *atom2);
    double radialDF(Atom *atom1, Atom *atom2);
    double BerendsenThermostat(double temp);
    void AndersenThermostat(Atom *atom);
protected:
    Cell* cells;
    vector<Atom*> atoms;
    Potentials* potential;
    int nsteps;
    double density;
    double volume;
    int N;
    int Nx, Ny, Nz;
    double sigma;
    double b, tau, k_b;
    double T0, T, T_bath, m, eV, epsilon;
    double standardDeviation;
    double L, dt;
    int Lc;
    double cellSize;
    int nCells;
    double r_cut;
    long idum;
    string atom_Name, thermostat;
};

#endif // GENERATEQUANTITIES_H
