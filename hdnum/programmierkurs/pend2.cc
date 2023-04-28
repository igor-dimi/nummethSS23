#include <iostream>
#include <cmath>

double l(1.34);
double phi(3.0); //note that the initial phi can be large
double u(0.0);
double dt(1E-4);
double T(10);
double t(0.0);

int main()
{

    while (t < T){
        std::cout << t << " " << phi << std::endl;
        double phiprev(phi);
        phi += dt * u;
        u -= dt * (9.81 / l) * sin(phiprev);
        t += dt;
    }
}
