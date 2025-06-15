#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// Initialize constants
const double gamma_val = 1.4;
const double end_time = 0.2;

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

double CalcSoundSpeed(const State& state)
{
    return sqrt(gamma_val * state.p / state.rho);
}

Conserved CalcRoeFlux(const State& left, const State& right)
{
    Conserved U_left = StateToConserved(left);
    Conserved U_right = StateToConserved(right);
    
    double sqrt_rho_left = sqrt(left.rho);
    double sqrt_rho_right = sqrt(right.rho);
    double rho_avg = sqrt_rho_left * sqrt_rho_right;
    double u_avg = (sqrt_rho_left * left.u + sqrt_rho_right * right.u) / (sqrt_rho_left + sqrt_rho_right);
    double H_left = (U_left.E + left.p) / left.rho;
    double H_right = (U_right.E + right.p) / right.rho;
    double H_avg = (sqrt_rho_left * H_left + sqrt_rho_right * H_right) / (sqrt_rho_left + sqrt_rho_right);
    double c_avg = sqrt((gamma_val - 1.0) * (H_avg - 0.5 * u_avg * u_avg));

    // Entropy modification
    double lambda1 = u_avg - c_avg;
    double lambda2 = u_avg;
    double lambda3 = u_avg + c_avg;
    double eps = 0.1 * c_avg;
    lambda1 = (abs(lambda1) < eps) ? (lambda1 * lambda1 + eps * eps) / (2.0 * eps) : abs(lambda1);
    lambda2 = (abs(lambda2) < eps) ? (lambda2 * lambda2 + eps * eps) / (2.0 * eps) : abs(lambda2);
    lambda3 = (abs(lambda3) < eps) ? (lambda3 * lambda3 + eps * eps) / (2.0 * eps) : abs(lambda3);

    double alpha[3] = {((right.p - left.p) - rho_avg * c_avg * (right.u - left.u)) / (2.0 * c_avg * c_avg),
                       (right.rho - left.rho) - (right.p - left.p) / (c_avg * c_avg),
                       ((right.p - left.p) + rho_avg * c_avg * (right.u - left.u)) / (2.0 * c_avg * c_avg)};

    // Calculate vectors for Roe average
    double r_vec1[3] = {1.0, u_avg - c_avg, H_avg - u_avg * c_avg};
    double r_vec2[3] = {1.0, u_avg, 0.5 * u_avg * u_avg};
    double r_vec3[3] = {1.0, u_avg + c_avg, H_avg + u_avg * c_avg};
    
    Conserved F;
    F.rho = lambda1 * alpha[0] * r_vec1[0] + lambda2 * alpha[1] * r_vec2[0] + lambda3 * alpha[2] * r_vec3[0];   
    F.rho_u = lambda1 * alpha[0] * r_vec1[1] + lambda2 * alpha[1] * r_vec2[1] + lambda3 * alpha[2] * r_vec3[1];   
    F.E = lambda1 * alpha[0] * r_vec1[2] + lambda2 * alpha[1] * r_vec2[2] + lambda3 * alpha[2] * r_vec3[2];

    // Calculate physical fluxes
    Conserved F_left, F_right;
    F_left.rho = left.rho * left.u;
    F_left.rho_u = left.rho * left.u * left.u + left.p;
    F_left.E = left.u * (U_left.E + left.p);
    F_right.rho = right.rho * right.u;
    F_right.rho_u = right.rho * right.u * right.u + right.p;
    F_right.E = right.u * (U_right.E + right.p);

    // Calculate Roe flux
    Conserved F_roe;
    F_roe.rho = 0.5 * (F_left.rho + F_right.rho) - 0.5 * F.rho;
    F_roe.rho_u = 0.5 * (F_left.rho_u + F_right.rho_u) - 0.5 * F.rho_u;
    F_roe.E = 0.5 * (F_left.E + F_right.E) - 0.5 * F.E;
    
    return F_roe;
}

void ComputeFluxes(const vector<State>& states, vector<Conserved>& fluxes, double dx)
{
    int n = states.size();
    fluxes.resize(n + 1);

    // Boundary conditions
    State left = states[0];
    left.u = -left.u;
    fluxes[0] = CalcRoeFlux(left, states[0]);
    State right = states[n-1];
    right.u = -right.u;
    fluxes[n] = CalcRoeFlux(states[n-1], right);

    // Calculate fluxes for interior states
    for (int i = 1; i < n; ++i)
    {
        fluxes[i] = CalcRoeFlux(states[i-1], states[i]);
    }
}

void RungeKutta3(vector<Conserved>& U, double dt, double dx)
{
    int n = U.size();
    vector<Conserved> U0 = U;
    vector<Conserved> U1(n), U2(n);
    vector<State> states(n);
    vector<Conserved> fluxes;

    for (int i = 0; i < n; ++i)
    {
        states[i] = ConservedToState(U0[i]);
    }
    ComputeFluxes(states, fluxes, dx);
    for (int i = 0; i < n; ++i)
    {
        double flux_factor = dt / dx;
        U1[i].rho = U0[i].rho - flux_factor * (fluxes[i+1].rho - fluxes[i].rho);
        U1[i].rho_u = U0[i].rho_u - flux_factor * (fluxes[i+1].rho_u - fluxes[i].rho_u);
        U1[i].E = U0[i].E - flux_factor * (fluxes[i+1].E - fluxes[i].E);
    }

    for (int i = 0; i < n; ++i)
    {
        states[i] = ConservedToState(U1[i]);
    }
    ComputeFluxes(states, fluxes, dx);
    for (int i = 0; i < n; ++i)
    {
        double flux_factor = dt / dx;
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - flux_factor * (fluxes[i+1].rho - fluxes[i].rho));
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - flux_factor * (fluxes[i+1].rho_u - fluxes[i].rho_u));
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - flux_factor * (fluxes[i+1].E - fluxes[i].E));
    }

    for (int i = 0; i < n; ++i)
    {
        states[i] = ConservedToState(U2[i]);
    }
    ComputeFluxes(states, fluxes, dx);
    for (int i = 0; i < n; ++i)
    {
        double flux_factor = dt / dx;
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - flux_factor * (fluxes[i+1].rho - fluxes[i].rho));
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - flux_factor * (fluxes[i+1].rho_u - fluxes[i].rho_u));
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - flux_factor * (fluxes[i+1].E - fluxes[i].E));
    }
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

        RungeKutta3(U, dt, dx);
        time += dt;
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
        cout << "x=" << i*dx << ": Density=" << states[i].rho << ", Velocity=" << states[i].u << ", Pressure=" << states[i].p << endl;
    }
    return 0;
}