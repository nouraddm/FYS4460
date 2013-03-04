#include "atom.h"

//Atom::Atom()
//{

//}

Atom::Atom(vec ri)
{
    r = ri;
}

vec3 Atom::getPosition()
{
    return r;
}

vec3 Atom::getVelocity()
{
    return v;
}

vec3 Atom::getForce()
{
    return f;
}

double Atom::getPotential()
{
    return potential;
}

double Atom::getKinetic()
{
    return kinetic;
}

double Atom::getPressure()
{
    return pressure;
}

void Atom::setPosition(const vec &inPosition)
{
    r = inPosition;
}

void Atom::setVelocity(const vec &inVelocity)
{
    v = inVelocity;
}

void Atom::setForce(const vec &inForce)
{
    f += inForce;
}

void Atom::addForce(const vec &force)
{
    f += force;
}

void Atom::setPotential(double inPotential)
{
    potential += inPotential;
}

void Atom::setKinetic(double inKinetic)
{
    kinetic = inKinetic;
}

void Atom::setPressure(double inPressure)
{
    pressure = inPressure;
}

void Atom::addPressure(double inPressure)
{
    pressure += inPressure;
}

vec3 Atom::getDisplacement()
{
    return d;
}

void Atom::setDisplacement(const vec &inDisplacement)
{
    d += inDisplacement;
}


