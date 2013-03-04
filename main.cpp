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


using namespace std;
using namespace arma;


int main()
{
    time_t start,end;
    time (&start);

    string command = "rm Argon6.xyz";
    system(command.c_str());

    string thermostat = "noThermostat";
    string thermostat1 = "Berendsen";
    string thermostat2 = "Andersen";

    GenerateQuantities app(thermostat1);

    app.generatePosition();
    app.generateVelocity();
    //app.printVelocity();
    //app.genTimeDevelepment();
    app.generateCells();

    time (&end);
    double dif = difftime (end,start);
    printf ("Elasped time is %.2lf minutes.", dif/60 );
}
