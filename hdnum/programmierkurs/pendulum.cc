#include <iostream>
#include <cmath>

#include <iostream>
#include <cmath>

double l(1.34); //Length of the pendulum chord in meters
double phi0(0.2); //Amplitude i.e. the initial angle in radians
double dt(0.05); //Time-step in seconds
double T(5.0); //End-time in seconds
double t(0.0); //Initial time value
int main(){
    while (t <= T){
        std::cout << t << " "
                << phi0 * cos(sqrt(9.81/l) * t)
                << std::endl;
        t += dt;
    }
}