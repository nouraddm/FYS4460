#include "generatequantities.h"
#include "potentials.h"

GenerateQuantities::GenerateQuantities(string nameOfThermostat, string half_Density, string Flow, string pore)
{
    potential = new Potentials();

    Nx = 8;
    Ny = 8;
    Nz = 8;
    N  = Nx*Ny*Nz*4;

    atom_Name = "Ar";
    half_density = half_Density;
    flow = Flow;
    shapeOfPores = pore;

    m     = 39.948*1.66053886*pow(10.0,-27.0);
    k_b   = 1.3806503*pow(10.0,-23.0);
    sigma = 3.405*pow(10.0,-10.0);
    if(half_density == "full density")
    {
        b = 5.72*pow(10.0,-10.0);

    }else if(half_density == "half density")
    {
        b = pow(2.0,1/3)*5.72*pow(10.0,-10.0);
    }
    //b     = 5.72*pow(10.0,-10.0);
    eV    = 1.60217646*pow(10.0,-19.0);

    epsilon = 1.0318*pow(10.0,-2.0) * eV;
    tau     = sigma * sqrt(m / epsilon);
    T0      = 119.74;
    T       = 1.5;
    T_bath  = 1.5;
    L = ((Nx) * (b/sigma));
    dt = 0.005;
    nsteps = 2000;
    nLpc = 10;
    Lpc = L / nLpc;
    nPressCells = nLpc * nLpc * nLpc;
    volumePressCells = Lpc * Lpc * Lpc;
    Lc = (int) (L / 3);
    volume = L*L*L;
    density = N / volume;
    r_cut = L / (1.0*Lc) + 0.0001;
    nCells = Lc * Lc * Lc;
    idum = -1;
    pore_radius = (1.0 * pow(10.0,-9.0)) / sigma;
    thermostat = nameOfThermostat;
    standardDeviation = sqrt(T);
    cell_center << L/2 << L/2 << L/2;
    who_moves = false;
    pi = acos(-1);
    F = 0.1 * epsilon / sigma;
    numberOfPores = 5;
    flowDirection = 2;
}

double GenerateQuantities::sampleNormal(long *idum)
{
    double u = ran0(idum) * 2 - 1;
    double v = ran0(idum) * 2 - 1;
    double r = u * u + v * v;
    if (r == 0 || r > 1) return sampleNormal(idum);
    double c = sqrt(-2 * log(r) / r);
    return u * c;
}

void GenerateQuantities::printToFile()
{
    string atom_name;
    ofstream myfile ("Argon.xyz", ios_base::app);
    if (myfile.is_open())
    {
        myfile << N << endl;
        myfile << "This is the comment line.\n";

        for(int i = 0; i < atoms.size(); i++)
        {
            Atom *A = atoms[i];
            if (A->inMotion == false){
                atom_name = "Ar";
            }
            if (A->inMotion == true){
                atom_name = "C";
            }
            myfile << atom_name << " " << A->getPosition()(0) << " " << A->getPosition()(1) << " " << A->getPosition()(2) << endl;
        }
        myfile.close();
    }
    else cout << "Unable to open file";
}

void GenerateQuantities::printVelocity()
{
    ofstream myfile ("Velocity.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < atoms.size(); i++)
        {
            Atom *A = atoms[i];
            myfile << A->getVelocity()(0) * (sigma/tau)<< " " << A->getVelocity()(1)  * (sigma/tau) << " " << A->getVelocity()(2)  * (sigma/tau)<< " "<< norm(A->getVelocity(),2)  * (sigma/tau)<<endl;
            }
        myfile.close();
    }
    else cout << "Unable to open file";
}

void GenerateQuantities::printPressure(const vec &inPressure)
{
    ofstream myfile ("Pressure.txt", ios_base::app);
    if (myfile.is_open())
    {
        for(int i = 0; i < nPressCells; i++)
        {
            myfile << inPressure(i)<< " ";
        }
        myfile << endl;
        myfile.close();
    }
    else cout << "Unable to open file";
}

void GenerateQuantities::printFlowProfile(const vec &inPosition)
{
    ofstream myfile ("Flow.txt", ios_base::app);
    if (myfile.is_open())
    {
        myfile << inPosition(0) << " " <<  inPosition(1) << " " << inPosition(2) << endl;
        myfile.close();
    }
    else cout << "Unable to open file";
}

void GenerateQuantities::generatePosition()
{
    Atom *A;

    vec R(3);
    vec r(3);

    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            for(int k = 0; k < Nz; k++)
            {
                R(0) = k*b/sigma;
                R(1) = j*b/sigma;
                R(2) = i*b/sigma;

                A = new Atom(R);
                atoms.push_back(A);

                r(0) = R(0) + 0.5*b/sigma;
                r(1) = R(1) + 0.5*b/sigma;
                r(2) = R(2);

                A = new Atom(r);
                atoms.push_back(A);

                r(0) = R(0);
                r(1) = R(1) + 0.5*b/sigma;
                r(2) = R(2) + 0.5*b/sigma;

                A = new Atom(r);
                atoms.push_back(A);

                r(0) = R(0) + 0.5*b/sigma;
                r(1) = R(1);
                r(2) = R(2) + 0.5*b/sigma;

                A = new Atom(r);
                atoms.push_back(A);
            }
        }
    }
}

void GenerateQuantities::generateVelocity()
{
    vec velocity(3);
    vec3 velocitySum = zeros(3);

    for(int i = 0; i < atoms.size(); i++)
    {
        // Get randomly distributed velocity,
//        velocity(0) = ran0(&idum);
//        velocity(1) = ran0(&idum);//(double) rand()/RAND_MAX;//
//        velocity(2) = ran0(&idum);
        // or, Gaussian distributed velocity
        velocity(0) = sampleNormal(&idum);
        velocity(1) = sampleNormal(&idum);
        velocity(2) = sampleNormal(&idum);
        //v /= norm(v,2);
        velocity *= standardDeviation;

        Atom* A = atoms.at(i);
        A->setVelocity(velocity);

        velocitySum += A->getVelocity();
    }
    velocitySum /= (1.0 * N);

    for(int i = 0; i < atoms.size(); i++)
    {
        Atom* A = atoms.at(i);
        velocity = A->getVelocity() - velocitySum;
        A->setVelocity(velocity);
        A->setKinetic(0.5*dot(velocity,velocity));
    }
    printVelocity();
}

void GenerateQuantities::generateForce()
{
    double normPosition;
    double Epot;
    vec force(3);
    vec DeltaR_ij(3);
    double dWork;

    for(int i = 0; i < atoms.size(); i++)
    {
        dWork = 0.0;
        for(int j = i+1; j < atoms.size(); j++)
        {
            DeltaR_ij = atoms[i]->getPosition() - atoms[j]->getPosition();
            for(int k = 0; k<3; k++)
            {
                if(abs(DeltaR_ij(k)-L) < abs(DeltaR_ij(k)))
                {
                     DeltaR_ij(k) = DeltaR_ij(k) - L;
                }
                if(abs(DeltaR_ij(k)+L) < abs(DeltaR_ij(k)))
                {
                    DeltaR_ij(k) =  DeltaR_ij(k) + L;
                }
            }
            normPosition = dot(DeltaR_ij, DeltaR_ij);
            force = 24.0 * (2.0 * pow(normPosition,-3.0) - 1.0 ) * DeltaR_ij * pow(normPosition,-4.0);
            Epot = 4 * (pow(normPosition,-6.0) - pow(normPosition,-3.0));
            dWork += dot(force,DeltaR_ij);

            atoms[i]->setForce(force);
            atoms[j]->setForce(-force);
            atoms[i]->setPotential(Epot);
            atoms[j]->setPotential(Epot);
        }
        atoms[i]->setPressure(dWork);
    }
}

// Verlet algorithm without cells
void GenerateQuantities::genTimeDevelepment()
{
    double dotPosition;
    vec force(3);
    vec acceleration(3);
    vec DeltaR_ij(3);
    vec oldPosition(3);
    vec newPosition(3);
    vec tempVelocity(3);
    vec newVelocity(3);

    generateForce();
    printToFile();

    for(int t = 1; t < nsteps; t++)
    {
        // Calculate temp. velocities and new positions for every atom
        for(int i = 0; i < atoms.size(); i++)
        {
            acceleration = atoms[i]->getForce();
            oldPosition = atoms[i]->getPosition();
            tempVelocity = atoms[i]->getVelocity() + 0.5 * acceleration * dt;
            newPosition = oldPosition + tempVelocity * dt;

            int k = 0;
            while(k<3)
            {
                if(newPosition(k) > L )
                {
                    newPosition(k) -= L;
                }
                if(newPosition(k) < 0 )
                {
                    newPosition(k) += L;
                }k++;
            }
            atoms[i]->setPosition(newPosition);
            atoms[i]->setVelocity(tempVelocity);
            atoms[i]->setForce(-acceleration);
        }

        // Calculate forces for all atoms in their new positions
        for(int i = 0; i < atoms.size(); i++)
        {
            for(int j = i+1; j < atoms.size(); j++)
            {
                DeltaR_ij = atoms[i]->getPosition() - atoms[j]->getPosition();
                for(int k = 0; k<3; k++)
                {
                    if(abs(DeltaR_ij(k)-L) < abs(DeltaR_ij(k)))
                    {
                        DeltaR_ij(k) = DeltaR_ij(k) - L;
                    }
                    if(abs(DeltaR_ij(k)+L) < abs(DeltaR_ij(k)))
                    {
                        DeltaR_ij(k) =  DeltaR_ij(k) + L;
                    }
                }
                //force = potential->Lennard_Jones_potential(DeltaR_ij);
                dotPosition = dot(DeltaR_ij, DeltaR_ij);
                force = 24.0 * (2.0 * pow(dotPosition,-3.0) - 1.0 ) * DeltaR_ij * pow(dotPosition,-4.0);
                atoms[i]->setForce(force);
                atoms[j]->setForce(-force);
            }
        }
        //Then calculate new velocities
        for(int i = 0; i < atoms.size(); i++)
        {
            newVelocity = atoms[i]->getVelocity() + 0.5 * atoms[i]->getForce() * dt;
            atoms[i]->setVelocity(newVelocity);
        }

        printToFile();
    }

}

void GenerateQuantities::generateCells()
{
    int ix, iy, iz;
    int nBins = 100;
    double Epot;
    double kinetic;
    double press;
    double r_sphereSum, r_cylinderSum;

    cells = new Cell[nCells];
    pressCells = new PressureCells[nPressCells];

    vec temp(2);
    vec r(3);
    vec d(3);
    vec r_center(3);
    vec cellIndices(3);
    vec indexes(3);
    vec3 newVelocity;
    vec3 oldPosition;
    vec3 newPosition;
    vec3 tempVelocity;
    vec3 force;
    vec kinetic_energy(nsteps);
    vec potential_energy(nsteps);
    vec total_energy(nsteps);
    vec temperature(nsteps);
    vec pressure(nsteps);
    vec displacement(nsteps);
    vec bins(nBins);
    vec sphere_radius(numberOfPores);
    vec cylinder_radius(numberOfPores);
    vec local_kinetic_energy(nPressCells);
    vec local_pressure(nPressCells);
    vec cellPressure(nPressCells);
    vec3 wFlow;
    wFlow << 0.0 << 0.0 << F;

    mat Index_Matrix;
    Index_Matrix = zeros<mat>(3,27);
    mat pore_Matrix;
    pore_Matrix = zeros<mat>(3,numberOfPores);
    mat pore_MatrixC;
    pore_MatrixC = zeros<mat>(2,numberOfPores);

    // Place every atom in the same position
    for (int i = 0; i < atoms.size(); i++)
    {
        atoms[i]->inMotion = false;
     }

    // Generate a matrix with 20 pores; the position of each pore is randomly genersated
    // The radius of each sphere is also randomly generated
    // The spheres do not overlap
    if(shapeOfPores == "spherical")
    {
        for(int j = 0; j < numberOfPores; j++)
        {
            //sphere_radius(j) = (ran0(&idum) * (pow(10.0,-9.0)/sigma)) + pore_radius;
            sphere_radius(j) = (((double) rand()/RAND_MAX) * (pow(10.0,-9.0)/sigma)) + pore_radius;

            for(int i = 0; i < 3; i++)
            {
                //pore_Matrix(i,j) = ran0(&idum) * L;
                pore_Matrix(i,j) = ((double) rand()/RAND_MAX) * L;
            }
            for(int k = 0; k < j; k++)
            {
                if( norm((pore_Matrix.col(k) - pore_Matrix.col(j)),2) <= sphere_radius(j) + sphere_radius(k))
                {
                    k = j;
                    j -= 1;
                }
            }
        }

        r_sphereSum = 0.0;
        for (int i = 0; i < numberOfPores; i++)
        {
            r_sphereSum += pow(sphere_radius(i),3.0);
        }
        porosity = 1.0 - (4/3) * (pi / volume) * r_sphereSum;
        cout << "Porosity: " << porosity << endl;

        for (int i = 0; i < atoms.size(); i++)
        {
            Atom *atom = atoms[i];
            newPosition = atom->getPosition();
            for(int l = 0; l < pore_Matrix.n_cols; l++)
            {
                r_center = newPosition - pore_Matrix.col(l);

                for(int k = 0; k < 3; k++)
                {
                    if(abs(r_center(k)-L) < abs(r_center(k)))
                    {
                        r_center(k) = r_center(k) - L;
                    }
                    if(abs(r_center(k)+L) < abs(r_center(k)))
                    {
                        r_center(k) = r_center(k) + L;
                    }
                }

                if(dot(r_center,r_center) <= sphere_radius(l) * sphere_radius(l))
                {
                    atom->inMotion = true;
                }
            }
        }
    }else if(shapeOfPores == "cylindrical")
    {
        cout << "cylinder" << endl;
        for(int j = 0; j < numberOfPores; j++)
        {
            //cylinder_radius(j) = (ran0(&idum) * pow(10.0,-9.0)/sigma) + pore_radius;
            cylinder_radius(j) = pore_radius;//(((double) rand()/RAND_MAX) * pow(10.0,-9.0)/sigma) + pore_radius;
            for(int i = 0; i < 2; i++)
            {
                //pore_MatrixC(i,j) = ran0(&idum) * L;
                pore_MatrixC(i,j) = 0.5 * L;//((double) rand()/RAND_MAX) * L;
            }
            for(int k = 0; k < j; k++)
            {
                if( norm((pore_MatrixC.col(k) - pore_MatrixC.col(j)),2) <= cylinder_radius(j) + cylinder_radius(k))
                {
                    cout << " k : "<<k<<endl;
                    k = j;
                    j -= 1;
                }
            }
        }
        r_cylinderSum = 0.0;
        for (int i = 0; i < numberOfPores; i++)
        {
            r_cylinderSum += pow(cylinder_radius(i),2.0);
        }
        porosity = 1.0 - (pi / volume) * r_cylinderSum * L;
        cout << "Porosity: " << porosity << endl;


        for (int i = 0; i < atoms.size(); i++)
        {
            Atom *atom = atoms[i];
            newPosition = atom->getPosition();
            for(int l = 0; l < pore_MatrixC.n_cols; l++)
            {
                temp << newPosition(0)<< newPosition(1);
                r_center = temp - pore_MatrixC.col(l);

                for(int k = 0; k < 2; k++)
                {
                    if(abs(r_center(k)-L) < abs(r_center(k)))
                    {
                        r_center(k) = r_center(k) - L;
                    }
                    if(abs(r_center(k)+L) < abs(r_center(k)))
                    {
                        r_center(k) = r_center(k) + L;
                    }
                }

                if(dot(r_center,r_center) > cylinder_radius(l) * cylinder_radius(l))
                {
                    atom->inMotion = true;
                }
            }
        }
    }

/* Generate the cells: Every cell has a number,
*  and is associated with three indices (i,j,k)
*
*  At the same time generate the pressure cells
*/
    int cellNumber = 0;
    for(int i = 0; i < Lc; i++)
    {
        for(int j = 0; j < Lc; j++)
        {
            for(int k = 0; k < Lc; k++)
            {
                cells[cellNumber].cellIndices << i << j << k;
                cellNumber++;
            }
        }
    }

    // Generate the pressure cells
    int pressCellNumber = 0;
    for(int i = 0; i < nLpc; i++)
    {
        for(int j = 0; j < nLpc; j++)
        {
            for(int k = 0; k < nLpc; k++)
            {
                pressCells[pressCellNumber].cell_index << i << j << k;

                pressCellNumber++;
            }
        }
    }

    // Construct an index matrix to use later to find
    // neighbor cells.
    int l= 0;
    for(int i = -1; i < 2; i++)
    {
        for(int j = -1; j < 2; j++)
        {
            for(int k = -1; k < 2; k++)
            {
                Index_Matrix(0,l) = i;
                Index_Matrix(1,l) = j;
                Index_Matrix(2,l) = k;

                l +=1;
            }
        }
    }

    // Find neighbor cells and add them to neighbor list
    for(int i = 0; i < nCells; i++)
    {
        for(int t = 0; t < Index_Matrix.n_cols; t++)
        {
            indexes = Index_Matrix.col(t) + cells[i].cellIndices;
            for(int k = 0; k < 3; k++)
            {
                if (indexes(k) > Lc-1)
                {
                    indexes(k) -= Lc;
                }
                else if (indexes(k) < 0)
                {
                    indexes(k) += Lc;
                }
            }
            for(int j = 0; j < nCells; j++)
            {
                if(cells[j].cellIndices(0)==indexes(0) &&
                        cells[j].cellIndices(1)==indexes(1) &&
                        cells[j].cellIndices(2)==indexes(2) && i!=j)
                {
                    cells[i].addNeighbor(&cells[j]);
                }
            }
        }
    }

    // Set atoms in cells
    for(int i = 0; i < atoms.size(); i++)
    {
        r = atoms[i]->getPosition();

        ix = int (r(0)/r_cut);
        iy = int (r(1)/r_cut);
        iz = int (r(2)/r_cut);

        cellIndices << ix << iy << iz;

        for(int j = 0; j < nCells; j++)
        {
            if(cells[j].cellIndices(0)==cellIndices(0) &&
                    cells[j].cellIndices(1)==cellIndices(1) &&
                    cells[j].cellIndices(2)==cellIndices(2))
            {
                cells[j].addAtoms(atoms[i]);
            }
        }

    }

    // Set atoms in pressure cells
    for(int i = 0; i < atoms.size(); i++)
    {
        r = atoms[i]->getPosition();

        ix = int (r(0)/Lpc);
        iy = int (r(1)/Lpc);
        iz = int (r(2)/Lpc);

        cellIndices << ix << iy << iz;

        for(int j = 0; j < nPressCells; j++)
        {
            if(pressCells[j].cell_index(0)==cellIndices(0) &&
                    pressCells[j].cell_index(1)==cellIndices(1) &&
                    pressCells[j].cell_index(2)==cellIndices(2))
            {
                pressCells[j].addAtoms(atoms[i]);
            }
        }
    }

//------------------------------------------------------------------------------
//     For every atom with a given position find the corrosponding indices,
//     then, add them to the cell list.
//------------------------------------------------------------------------------
    generateForce();
    printToFile();

    // Intialize bins for g(r)
    for (int i = 0; i < nBins; i++)
    {
        bins(i) = 0;
    }


    // Calculate the kinetic energy and thereby the temperature at t = 0.
    kinetic_energy(0)   = 0.0;
    potential_energy(0) = 0.0;
    displacement(0)     = 0.0;
    pressure(0)         = 0.0;
    for(int i = 0; i < atoms.size(); i++)
    {
        kinetic_energy(0) += atoms[i]->getKinetic();
        potential_energy(0) += atoms[i]->getPotential();
        pressure(0) += atoms[i]->getPressure();
    }
    temperature(0) = T;//(2.0 * kinetic_energy(0)) / (3.0 * N);
    total_energy(0) = kinetic_energy(0) + potential_energy(0);
    pressure(0) = density * temperature(0) + pressure(0) / (3.0 * volume);


    // Calculate the pressure in every cell
    for(int cellNumber = 0; cellNumber < nPressCells; cellNumber++)
    {
        local_kinetic_energy(cellNumber) = 0.0;
        local_pressure(cellNumber) = 0.0;
        for(int i = 0; i < (pressCells[cellNumber].atoms).size(); i++)
        {
            local_kinetic_energy(cellNumber) += pressCells[cellNumber].atoms[i]->getKinetic();
            local_pressure(cellNumber) += pressCells[cellNumber].atoms[i]->getPressure();
        }
        cellPressure(cellNumber) = ( 2 * local_kinetic_energy(cellNumber) + local_pressure(cellNumber));
        cellPressure(cellNumber) /= (3 * volumePressCells);
    }

    printPressure(cellPressure);


//---------------------------------------------------------------------------------------
//                        ** Time development **
//---------------------------------------------------------------------------------------
    for(int t = 1; t < nsteps; t++)
    {
        cout << "t: " << t << endl;

        // Calculate temp. velocities and new positions for every atom
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++)
        {
            for(int j = 0; j < (cells[cellNumber].atoms).size(); j++)
            {
                force = cells[cellNumber].atoms[j]->getForce();
                tempVelocity = cells[cellNumber].atoms[j]->getVelocity() + 0.5 * force * dt;
                oldPosition = cells[cellNumber].atoms[j]->getPosition();
                newPosition = oldPosition + tempVelocity * dt;
                d = newPosition - oldPosition;

                if (cells[cellNumber].atoms[j]->inMotion == who_moves)
                {
                    int k = 0;
                    while(k < 3)
                    {
                        if(newPosition(k) > L )
                        {
                            newPosition(k) -= L;
                            if (k == flowDirection)
                            {
                                printFlowProfile(newPosition);
                            }
                        }
                        if(newPosition(k) < 0 )
                        {
                            newPosition(k) += L;
                            if (k == flowDirection)
                            {
                                printFlowProfile(newPosition);
                            }
                        }k++;
                    }
                    cells[cellNumber].atoms[j]->setPosition(newPosition);
                    cells[cellNumber].atoms[j]->setDisplacement(d);
                    cells[cellNumber].atoms[j]->setVelocity(tempVelocity);
                }
            }
        }

        //clear Atoms In Cells();
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++){
              cells[cellNumber].atoms.clear();
          }

        //clear Atoms in pressure cells();
        for(int cellNumber = 0; cellNumber < nPressCells; cellNumber++){
              pressCells[cellNumber].atoms.clear();
          }

        // Set atoms in pressure cells
        for(int i = 0; i < atoms.size(); i++)
        {
            r = atoms[i]->getPosition();

            ix = int (r(0)/Lpc);
            iy = int (r(1)/Lpc);
            iz = int (r(2)/Lpc);

            cellIndices << ix << iy << iz;

            for(int j = 0; j < nPressCells; j++)
            {
                if(pressCells[j].cell_index(0)==cellIndices(0) &&
                        pressCells[j].cell_index(1)==cellIndices(1) &&
                        pressCells[j].cell_index(2)==cellIndices(2))
                {
                    pressCells[j].addAtoms(atoms[i]);
                }
            }
        }

        // Set atoms in the cells
        for(int i = 0; i < atoms.size(); i++)
        {
            r = atoms[i]->getPosition();

            ix = int (r(0)/r_cut);
            iy = int (r(1)/r_cut);
            iz = int (r(2)/r_cut);

            cellIndices << ix << iy << iz;

            for(int j = 0; j < nCells; j++)
            {
                if(cells[j].cellIndices(0)==cellIndices(0) &&
                   cells[j].cellIndices(1)==cellIndices(1) &&
                   cells[j].cellIndices(2)==cellIndices(2))
                {
                    cells[j].addAtoms(atoms[i]);
                }
            }
        }

        // Reset the force, the potential and pressure for every atom and then calculate the force between atoms
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++)
        {
            cells[cellNumber].hasCalculated = false;
            for(int j = 0; j < (cells[cellNumber].atoms).size(); j++)
            {
                Atom *atom = cells[cellNumber].atoms[j];
                force = atom->getForce();
                Epot = atom->getPotential();
                kinetic = 0.0;
                press = 0.0;

                if(flow == "withFlow")
                {
                    cells[cellNumber].atoms[j]->setForce(-force);
                    cells[cellNumber].atoms[j]->setForce(wFlow);
                }else if (flow == "noFlow")
                {
                    cells[cellNumber].atoms[j]->setForce(-force);
                }
                cells[cellNumber].atoms[j]->setPotential(-Epot);
                cells[cellNumber].atoms[j]->setKinetic(kinetic);
                cells[cellNumber].atoms[j]->setPressure(press);
            }

        }


        // Calculate the forces between atoms in the same cell
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++)
        {
            for(int i = 0; i < (cells[cellNumber].atoms).size(); i++)
            {
                for(int j = i+1; j < (cells[cellNumber].atoms).size(); j++)
                {
                    Atom *atom1 = cells[cellNumber].atoms[i];
                    Atom *atom2 = cells[cellNumber].atoms[j];
                    if(atom1->inMotion != who_moves && atom2->inMotion != who_moves)
                    {
                    }else{
                        acceleration(atom1, atom2);
                    }
                }
            }
            // Calculate the forces between atoms in the neighboring cells
            for(int neighborNumber = 0; neighborNumber < (cells[cellNumber].neighbor).size(); neighborNumber++)
            {
                Cell* neighbor = cells[cellNumber].neighbor[neighborNumber];
                if(!neighbor->hasCalculated==false)
                {
                    for(int i = 0; i < (cells[cellNumber].atoms).size(); i++)
                    {
                        for(int j = 0; j < (neighbor->atoms).size(); j++)
                        {
                            Atom *atom1 = cells[cellNumber].atoms[i];
                            Atom *atom2 = neighbor->atoms[j];
                            if(atom1->inMotion != who_moves && atom2->inMotion != who_moves )
                            {
                            }else{
                                acceleration(atom1, atom2);
                            }
                        }
                    }
                }
            }
            cells[cellNumber].hasCalculated = true;
        }


        // Calculate new velocities
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++)
        {
            for(int j = 0; j < (cells[cellNumber].atoms).size(); j++)
            {
                Atom *atom = cells[cellNumber].atoms[j];

                if(atom->inMotion == who_moves)
                {
                    newVelocity = atom->getVelocity() + 0.5 * atom->getForce() * dt;
                    cells[cellNumber].atoms[j]->setVelocity(newVelocity);
                    cells[cellNumber].atoms[j]->setKinetic(dot(newVelocity,newVelocity)*0.5);
                }

            }
        }
        // Calculate g(r)
        for (int i = 0; i < atoms.size(); i++)
        {
            for(int j = 0; j < atoms.size(); j++)
            {
                Atom *atom1 = atoms[i];
                Atom *atom2 = atoms[j];
                double dr = radialDF(atom1,atom2);
                int BinIndex = nBins * (int) 2*dr/L;
                if (BinIndex >= nBins)
                {

                }else{
                    bins(BinIndex) += 1;
                }
            }
        }

        // Calculate system dynamics
        kinetic_energy(t)   = 0.0;
        potential_energy(t) = 0.0;
        pressure(t)         = 0.0;
        displacement(t)     = 0.0;
        for(int i = 0; i < atoms.size(); i++)
        {
            kinetic_energy(t) += atoms[i]->getKinetic();
            potential_energy(t) += atoms[i]->getPotential();
            pressure(t) += atoms[i]->getPressure();
            displacement(t) += dot(atoms[i]->getDisplacement(), atoms[i]->getDisplacement());
        }

        total_energy(t) = kinetic_energy(t) + potential_energy(t);
        temperature(t) = (2.0 * kinetic_energy(t)) / (3.0 * N);
        pressure(t) = density * temperature(t) + pressure(t) / (3.0 * volume);
        displacement(t) /= N;

        // Calculate the pressure in every cell
        for(int cellNumber = 0; cellNumber < nPressCells; cellNumber++)
        {
            local_kinetic_energy(cellNumber) = 0.0;
            local_pressure(cellNumber) = 0.0;
            for(int i = 0; i < (pressCells[cellNumber].atoms).size(); i++)
            {
                local_kinetic_energy(cellNumber) += pressCells[cellNumber].atoms[i]->getKinetic();
                local_pressure(cellNumber) += pressCells[cellNumber].atoms[i]->getPressure();
            }
            cellPressure(cellNumber) = ( 2 * local_kinetic_energy(cellNumber) + local_pressure(cellNumber));
            cellPressure(cellNumber) /= (3 * volumePressCells);
        }

        printPressure(cellPressure);


        //Then calculate new velocities
        for(int cellNumber = 0; cellNumber < nCells; cellNumber++)
        {
            for(int j = 0; j < (cells[cellNumber].atoms).size(); j++)
            {
                Atom *atom = cells[cellNumber].atoms[j];
                if(atom->inMotion == who_moves)
                {
                    newVelocity = atom->getVelocity();

                    if(thermostat == "noThermostat")
                    {
                        cells[cellNumber].atoms[j]->setVelocity(newVelocity);
                    }else if(thermostat == "Berendsen")
                    {
                        if(t==1500)
                        {
                            thermostat = "noThermostat";
                        }else{
                            double tempBerendsen = BerendsenThermostat(temperature(t));
                            newVelocity *= tempBerendsen;
                            cells[cellNumber].atoms[j]->setVelocity(newVelocity);
                        }

                    }else if(thermostat == "Andersen")
                    {
                        Atom *atom1 = cells[cellNumber].atoms[j];
                        AndersenThermostat(atom1);
                    }
                }
            }
        }
         printToFile();

    }
//-------------------------------------------------------------------------------------------
//                        The end of time development
//-------------------------------------------------------------------------------------------

    // Radial distribution function g(r):
    for (int i = 0; i < nBins; i++)
    {
        bins(i) /= nsteps;
    }

    // Write dynamical properties to file
    ofstream myfile ("Energy.txt");
    if (myfile.is_open())
    {
        for(int t = 0; t < nsteps; t++)
        {
            myfile << t << " " << total_energy(t) << " " << temperature(t) << " "<< kinetic_energy(t)
                   << " " << potential_energy(t) << " " << pressure(t)<< " " << displacement(t) << endl;
//            myfile << t << " " << temperature(t) << " "<< pressure(t) << endl;
        }
        myfile.close();
    }
    else cout << "Unable to open file";

    ofstream myfil ("radial.txt");
    if (myfil.is_open())
    {
        for(int t = 0; t < nBins; t++)
        {
            myfil << t << "  " << bins(t) << endl;
        }
        myfil.close();
    }
    else cout << "Unable to open file";

}

void GenerateQuantities::acceleration(Atom *atom1, Atom *atom2)
{
    double Ep;
    double p;
    double innerProduct, drSixthPower, drTwelvethPower;
    vec r1(3);
    vec r2(3);
    vec dr(3);
    vec force(3);

    r1 = atom1->getPosition();
    r2 = atom2->getPosition();

    dr = r1 - r2;
    for(int k = 0; k<3; k++)
    {
        if(abs(dr(k)-L) < abs(dr(k)))
        {
             dr(k) = dr(k) - L;
        }
        if(abs(dr(k)+L) < abs(dr(k)))
        {
            dr(k) = dr(k) + L;
        }
    }
    //Ep = 4.0 * (pow(dot(dr,dr),-6.0) - pow(dot(dr,dr),-3.0));
    innerProduct = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2);
    drSixthPower = innerProduct * innerProduct * innerProduct;
    drTwelvethPower = drSixthPower * drSixthPower;

    Ep = 4.0 * ((1.0 / drTwelvethPower) - (1.0 / drSixthPower));
    force = potential->Lennard_Jones_potential(dr);

    //p = dot(force, dr);
    p = force(0)*dr(0) + force(1)*dr(1) + force(2)*dr(2);

    atom1->setForce(force);
    atom2->setForce(-force);
    atom1->setPotential(Ep);
    atom2->setPotential(Ep);
    atom1->addPressure(p);

}

double GenerateQuantities::radialDF(Atom *atom1, Atom *atom2)
{
    vec r1(3);
    vec r2(3);
    vec dr(3);
    r1 = atom1->getPosition();
    r2 = atom2->getPosition();

    dr = r1 - r2;
    for(int k = 0; k<3; k++)
    {
        if(abs(dr(k)-L) < abs(dr(k)))
        {
             dr(k) = dr(k) - L;
        }
        if(abs(dr(k)+L) < abs(dr(k)))
        {
            dr(k) = dr(k) + L;
        }
    }
    return norm(dr,2);
}

double GenerateQuantities::BerendsenThermostat(double temp)
{
    double tau = 15 * dt;
    double gamma = sqrt( 1.0 + (dt /tau) * ((T_bath / temp) - 1.0));
    return gamma;
}


void GenerateQuantities::AndersenThermostat(Atom *atom)
{
    double random_number = (double) rand()/RAND_MAX;
    double tau = 15 * dt;
    double std = sqrt(T_bath);

    vec velocity(3);

    if(random_number < (dt/tau))
    {
        velocity(0) = sampleNormal(&idum);
        velocity(1) = sampleNormal(&idum);
        velocity(2) = sampleNormal(&idum);
        //velocity /= norm(velocity,2);
        velocity *= std;
        atom->setVelocity(velocity);
    }
}
