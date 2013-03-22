#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <cmath>
#include "math.h"
#include "lib.h"
#include "atom.h"
#include <../../../Desktop/FYS4460/ProjectOne/include/armadillo>
#include "generatequantities.h"
#include "potentials.h"
#include <time.h>
//#include <mpi.h>


using namespace std;
using namespace arma;


int main(int argc, char** argv)
{
    //int numprocs;
    //MPI_Init(&argc, &argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //cout << "Number of processes " << numprocs << endl;


    time_t start,end;
    time (&start);

    string command = "rm Argon.xyz";
    system(command.c_str());

    string thermostat = "noThermostat";
    string thermostat1 = "Berendsen";
    string thermostat2 = "Andersen";
    string half_density = "half density";
    string full_density = "full density";
    string flow = "withFlow";
    string noFlow = "noFlow";
    string sphericalPore = "spherical";
    string cylindricalPore = "cylindrical";

    GenerateQuantities app(thermostat1, half_density, noFlow, sphericalPore);

    app.generatePosition();
    app.generateVelocity();
    //app.printVelocity();
    //app.genTimeDevelepment();
    app.generateCells();

    time (&end);
    double dif = difftime (end,start);
    printf ("Elasped time is %.2lf minutes.", dif/60 );
    //MPI_Finalize();
}
