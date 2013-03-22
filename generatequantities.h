#ifndef GENERATEQUANTITIES_H
#define GENERATEQUANTITIES_H

#include <iostream>
#include <string>
#include <../../../Desktop/FYS4460/ProjectOne/include/armadillo>
#include "atom.h"
#include "cell.h"
#include "pressurecells.h"
#include "lib.h"

class Potentials;

using namespace std;
using namespace arma;

class GenerateQuantities
{
public:
    GenerateQuantities(string nameOfThermostat, string half_Density, string Flow, string pore);
    double sampleNormal(long *idum);
    void printToFile();
    void printVelocity();
    void printPressure(const vec &inPressure);
    void printFlowProfile(const vec &inPosition);
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
    PressureCells* pressCells;
    vector<Atom*> atoms;
    Potentials* potential;
    int nsteps;
    double density;
    double volume, volumePressCells;
    int N;
    int Nx, Ny, Nz;
    double sigma;
    double b, tau, k_b;
    double T0, T, T_bath, m, eV, epsilon;
    double standardDeviation;
    double L, dt;
    int Lc, nLpc;
    double cellSize, Lpc;
    int nCells, nPressCells;
    double r_cut, ;
    long idum;
    string atom_Name, thermostat, half_density, flow, shapeOfPores;
    double pore_radius;
    vec3 cell_center;
    bool who_moves;
    double porosity;
    double pi, F;
    int numberOfPores, flowDirection;
};

#endif // GENERATEQUANTITIES_H
