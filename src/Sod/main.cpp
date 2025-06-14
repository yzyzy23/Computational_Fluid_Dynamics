#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// Initialize constants
const double gamma = 1.4;
const double end_time = 1.0;

struct State
{
    double rho, u, p;
};
struct Conserved
{
    double rho, rho_u, E;
};

Conserved CalcFlux(const Conserved& U_left, const Conserved& U_right)
{
    Conserved F;
    F.rho = 0.5 * (U_left.rho + U_right.rho);
    F.rho_u = 0.5 * (U_left.rho_u + U_right.rho_u);
    F.E = 0.5 * (U_left.E + U_right.E);
    return F;
}
void RungeKutta3(vector<Conserved>& U, double dt)
{
    int n = U.size();
    vector<Conserved> U0 = U;
    vector<Conserved> U1(n), U2(n), U3(n);

    for (int i = 0; i < n; ++i)
    {
        Conserved F = CalcFlux(U0[i], U0[i + 1]);
        U1[i].rho = U0[i].rho + dt * F.rho;
        U1[i].rho_u = U0[i].rho_u + dt * F.rho_u;
        U1[i].E = U0[i].E + dt * F.E;
    }
    for (int i = 0; i < n; ++i)
    {
        Conserved F = CalcFlux(U1[i], U1[i + 1]);
        U2[i].rho = (3.0 * U0[i].rho + U1[i].rho + dt * F.rho) / 4.0;
        U2[i].rho_u = (3.0 * U0[i].rho_u + U1[i].rho_u + dt * F.rho_u) / 4.0;
        U2[i].E = (3.0 * U0[i].E + U1[i].E + dt * F.E) / 4.0;
    }
    for (int i = 0; i < n; ++i)
    {
        Conserved F = CalcFlux(U2[i], U2[i + 1]);
        U3[i].rho = (U0[i].rho + 2.0 * U2[i].rho + 2.0 * dt * F.rho) / 3.0;
        U3[i].rho_u = (U0[i].rho_u + 2.0 * U2[i].rho_u + 2.0 * dt * F.rho_u) / 3.0;
        U3[i].E = (U0[i].E + 2.0 * U2[i].E + 2.0 * dt * F.E) / 3.0;
    }
    for (int i = 0; i < n; ++i)
    {
        U[i] = U3[i];
    }
}

int main()
{
    // Initialize parameters
    double rho_left = 1.0;
    double u_left = 0.0;
    double p_left = 1.0;
    double rho_right = 0.125;
    double u_right = 0.0;
    double p_right = 0.1;
    double Cv = 4.138;
    double dx, dt;
    cin >> dx >> dt;
    double sigma = dt / dx;
}