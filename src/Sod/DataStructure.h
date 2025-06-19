#include <iostream>
#include <cmath>

using namespace std;

// Initialize constants
const double gamma_val = 1.4;
const double R_val = 286.9;
const double cfl = 0.2;
const double x_min = -1.0;
const double x_max = 1.0;
const double eps = 1e-5; // Small value to avoid division by zero

struct State {double rho, u, p;};
struct Conserved {double rho, rho_u, E;};

Conserved StateToConserved(const State& state)
{
    Conserved U;
    U.rho = state.rho;
    U.rho_u = state.rho * state.u;
    U.E = state.p / (gamma_val - 1.0) + 0.5 * state.rho * state.u * state.u;
    return U;
}

State ConservedToState(const Conserved& U)
{
    State state;
    state.rho = U.rho;
    state.u = U.rho_u / U.rho;
    state.p = (U.E - 0.5 * U.rho * state.u * state.u) * (gamma_val - 1.0);
    return state;
}