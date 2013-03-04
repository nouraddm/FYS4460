#include "lattice.h"

Lattice::Lattice()
{
}

vec Lattice::getCell()
{
    return cells;
}

void Lattice::setCell(const vec &inCell)
{
    cells = inCell;
}

//void Lattice::box()
//{
//    int ix, iy, iz;
//    // int Lc = (int) (L / r_cut);
//    // int nCells = pow(Lc,3);

//    vec r(3);

//    // "cells" gir boksene
//    // "pointer" gir lista med atomer i en gitt celle

//    for(int c = 0; c < nCells; c++)
//    {
//        cells[c] = EMPTY;
//    }

//    for(int i = 0; i < list.size(); i++)
//    {
//        r = list[i]->getPosition();

//        ix = int (r(0)/3.0);
//        iy = int (r(1)/3.0);
//        iz = int (r(2)/3.0);

//        pointer[i] = cells[boxNumber(ix,iy,iz)];
//        cells[boxNumber(ix,iy,iz)] = i;
//    }
//}


//int Lattice::boxNumber(int i, int j, int k)
//{
//    int indexNumber;
//    int lenght, height;
//    lenght = 8;
//    height = 8;

//    indexNumber = (i * lenght * height) + (j * height) + k;
//    return indexNumber;
//}
