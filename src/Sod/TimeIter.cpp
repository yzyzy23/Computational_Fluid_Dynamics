#include <iostream>
#include <vector>
#include "CalcFlux.cpp"

using namespace std;

void RungeKutta3FDS(vector<Conserved>& U, double dt, double dx)
{
    int n = U.size();
    vector<Conserved> U0 = U;
    vector<Conserved> U1(n), U2(n);
    vector<Conserved> fluxes(n);
    vector<Conserved> f(n);
    for (int i = 0; i < n; ++i)
    {
        f[i] = {0.0, 0.0, 0.0};
        fluxes[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFDS(U, dx, fluxes, f);
    for (int i = 0; i < n; ++i)
    {
        U1[i].rho = U0[i].rho - dt * f[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * f[i].rho_u;
        U1[i].E = U0[i].E - dt * f[i].E;
    }
    CalcFluxFDS(U1, dx, fluxes, f);
    for (int i = 0; i < n; ++i)
    {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * f[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * f[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * f[i].E);
    }
    CalcFluxFDS(U2, dx, fluxes, f);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * f[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * f[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * f[i].E);
    }
}
void RungeKutta3FVS(vector<Conserved>& U, double dt, double dx)
{
    int n = U.size();
    vector<Conserved> U0 = U;
    vector<Conserved> U1(n), U2(n);
    vector<Conserved> f_left(n), f_right(n), diff(n), temp_left(n), temp_right(n);
    for (int i = 0; i < n; ++i)
    {
        f_left[i] = {0.0, 0.0, 0.0};
        f_right[i] = {0.0, 0.0, 0.0};
        diff[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFVS(U, f_left, f_right);
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, 1);
    for (int i = 0; i < n; ++i)
    {
        U1[i].rho = U0[i].rho - dt * diff[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * diff[i].rho_u;
        U1[i].E = U0[i].E - dt * diff[i].E;
    }
    CalcFluxFVS(U1, f_left, f_right);
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, 1);
    for (int i = 0; i < n; ++i)
    {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * diff[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * diff[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * diff[i].E);
    }
    CalcFluxFVS(U2, f_left, f_right);
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, 1);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * diff[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * diff[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * diff[i].E);
    }
}
void RungeKutta3(vector<Conserved>& U, double dt, double dx, int flux_process_type)
{
    if (flux_process_type == 0)
    {
        // Use FDS flux processing
        RungeKutta3FDS(U, dt, dx);
    }
    else if (flux_process_type == 1)
    {
        // Use FVS flux processing
        RungeKutta3FVS(U, dt, dx);
    }
    else
    {
        // Handle other flux processing types if needed
        cout << "Unsupported flux process type: " << flux_process_type << endl;
    }
}