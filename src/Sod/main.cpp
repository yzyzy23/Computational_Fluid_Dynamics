#include <iostream>
#include <vector>
#include <cmath>
#include "TimeIter.cpp"
using namespace std;

// Initialize constants
const double end_time = 0.14; // End time for the simulation

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
    }
    /*State temp = states[1];
    states[1] = states[0];
    states[0] = temp;
    states[n - 1] = states[n - 4];
    states[n - 2] = states[n - 3];*/
    for (int i = 0; i < n; ++i)
    {
        U[i] = StateToConserved(states[i]);
    }
}

int main() 
{
    double dx, dt = 0.0;
    cout << "Enter dx : ";
    cin >> dx;
    cout << "Select flux process type (0 for FDS, 1 for FVS): ";
    int flux_process_type;
    cin >> flux_process_type;
    cout << "Select difference calculation type (0 for TVD, 1 for WENO, 2 for GVC): ";
    int diff_calc_type;
    cin >> diff_calc_type;
    int n = static_cast<int>(1.0 / dx);
    //cout << "CFL condition: " << dt * sqrt(gamma_val) / dx << endl;

    vector<State> states(n);
    vector<Conserved> U(n);
    InitSodProblem(states, U);

    double max_speed = 0;
    double time = 0.0;
    while (time < end_time)
    {
        max_speed = 0;
        for (int i = 0; i < n; ++i)
        {
            max_speed = max(max_speed, abs(CalcSoundSpeed(states[i]) + states[i].u));
        }
        dt = cfl * dx/ max_speed;
        if (time + dt > end_time)
        {
            dt = end_time - time;
        }
        cout << dt << endl;

        time += dt;
        RungeKutta3(U, dt, dx, flux_process_type, diff_calc_type, n);
        
        for (int i = 0; i < n; ++i)
        {
            states[i] = ConservedToState(U[i]);
        }
    }

    cout << "\nFinal state at t = " << time << "s:" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << "x=" << i*dx << ": Density=" << states[i].rho << ", Velocity=" << states[i].u << ", Pressure=" << states[i].p << endl;
    }

    // Save results to file
    FILE* file = fopen("SodResult.txt", "w");
    if (file == NULL)
    {
        cout << "Error opening file: SodResult.txt!" << endl;
        return -1;
    }
    for (int i = 0; i < n; ++i)
    {
        fprintf(file, "%lf %lf %lf %lf\n", i * dx, states[i].rho, states[i].u, states[i].p);
    }
    fclose(file);

    return 0;
}