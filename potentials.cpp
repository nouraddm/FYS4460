#include "potentials.h"

Potentials::Potentials()
{
    zeroVector = zeros(3);
}

vec Potentials::Zero_potential(const vec3 &dr)
{
    return zeroVector;
}

vec Potentials::Lennard_Jones_potential(const vec3 &dr)
{
    //innerProduct = dot(dr,dr);
    //f = 24.0 * (2.0 * pow(innerProduct,-3.0) - 1.0 ) * dr * pow(innerProduct,-4.0);

    innerProduct = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2);
    double drSixthPower = innerProduct * innerProduct * innerProduct;
    double drTwelvethPower = drSixthPower * drSixthPower;
    f = 24.0 * ((2.0 / drTwelvethPower) - (1.0 / drSixthPower) ) * (dr / innerProduct);

    return f;
}
