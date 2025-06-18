#include <iostream>
#include <vector>
#include <cmath>
#include "TimeIter.cpp"
using namespace std;

// Initialize constants
const double end_time = 0.2; // End time for the simulation

double CalcSoundSpeed(const State& state)
{
    return sqrt(gamma_val * state.p / state.rho);
}

void InitSodProblem(vector<State>& states, vector<Conserved>& U) {
    int n = states.size();
    for (int i = 0; i < n; ++i)
    {
        if (i < n / 2)
        {
            states[i] = {1.0, 0.0, 1.0};
        }
        else 
        {
            states[i] = {0.125, 0.0, 0.1};
        }
        U[i] = StateToConserved(states[i]);
    }
}

int main() {
    double dx, dt;
    cout << "Enter dx and dt: ";
    cin >> dx >> dt;
    int n = static_cast<int>(1.0 / dx) + 1;

    vector<State> states(n);
    vector<Conserved> U(n);
    InitSodProblem(states, U);

    double time = 0.0;
    int step = 0;
    while (time < end_time) {
        if (time + dt > end_time)
        {
            dt = end_time - time;
        }

        time += dt;
        RungeKutta3(U, dt, dx, 1);
        step++;
        
        for (int i = 0; i < n; ++i)
        {
            states[i] = ConservedToState(U[i]);
        }
        
        if (step % 50 == 0)
        {
            cout << "Step: " << step << " Time: " << time << " dt: " << dt << endl;
        }
    }

    cout << "\nFinal state at t = " << time << "s:" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << "x=" << i*dx << ": Density=" << U[i].rho << ", Velocity=" << states[i].u << ", Pressure=" << states[i].p << endl;
    }
    return 0;
}