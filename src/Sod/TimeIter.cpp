#include <iostream>
#include <vector>
#include "CalcFlux.cpp"

using namespace std;

void ApplyBoundaryConditions(vector<Conserved>& U)
{
    U[0].rho = U[1].rho;
    U[0].rho_u = -U[1].rho_u;
    U[0].E = U[1].E;

    int n = U.size();
    U[n - 1].rho = U[n - 2].rho;
    U[n - 1].rho_u = U[n - 2].rho_u;
    U[n - 1].E = U[n - 2].E;
}
void add_ghost_cells(const vector<Conserved>& u, vector<Conserved>& u_ext)
{
    u_ext.push_back(u[0]);
    u_ext.push_back(u[0]);
    u_ext.insert(u_ext.end(), u.begin(), u.end());
    u_ext.push_back(u.back());
    u_ext.push_back(u.back());
}
void compute_dF(const vector<Conserved>& u_ext, vector<Conserved>& dF, int N, double dx)
{
    vector<Conserved> fluxes;
    CalcFluxFDS(u_ext, fluxes);
    dF.resize(N + 1);

    for (int i = 2; i < N + 3; ++i)
    {
    int idx_dF = i - 2;
    dF[idx_dF].rho = (fluxes[i].rho - fluxes[i-1].rho) / dx;
    dF[idx_dF].rho_u = (fluxes[i].rho_u - fluxes[i-1].rho_u) / dx;
    dF[idx_dF].E = (fluxes[i].E - fluxes[i-1].E) / dx;
    }
}
void RungeKutta3FDS(vector<Conserved>& u_conservation, double dt, int n, double dx)
{
    vector<Conserved> u_ext;
    add_ghost_cells(u_conservation, u_ext);
    vector<Conserved> dF1(u_conservation.size());
    compute_dF(u_ext, dF1, n, dx);
    
    vector<Conserved> u1 = u_conservation;
    for (int i = 0; i < u1.size(); ++i)
    {
        u1[i].rho -= dt * dF1[i].rho;
        u1[i].rho_u -= dt * dF1[i].rho_u;
        u1[i].E -= dt * dF1[i].E;
    }
    
    vector<Conserved> u1_ext;
    add_ghost_cells(u1, u1_ext);
    vector<Conserved> dF2(u_conservation.size());
    compute_dF(u1_ext, dF2, n, dx);
    
    vector<Conserved> u2 = u_conservation;
    for (int i = 0; i < u2.size(); ++i)
    {
        u2[i].rho = 0.75 * u_conservation[i].rho + 
                    0.25 * (u1[i].rho - dt * dF2[i].rho);
        u2[i].rho_u = 0.75 * u_conservation[i].rho_u + 
                      0.25 * (u1[i].rho_u - dt * dF2[i].rho_u);
        u2[i].E = 0.75 * u_conservation[i].E + 
                  0.25 * (u1[i].E - dt * dF2[i].E);
    }
    
    vector<Conserved> u2_ext;
    add_ghost_cells(u2, u2_ext);
    vector<Conserved> dF3(u_conservation.size());
    compute_dF(u2_ext, dF3, n, dx);
    
    for (int i = 0; i < u_conservation.size(); ++i)
    {
        u_conservation[i].rho = u_conservation[i].rho / 3.0 + 
            (2.0/3.0) * (u2[i].rho - dt * dF3[i].rho);
        u_conservation[i].rho_u = u_conservation[i].rho_u / 3.0 + 
            (2.0/3.0) * (u2[i].rho_u - dt * dF3[i].rho_u);
        u_conservation[i].E = u_conservation[i].E / 3.0 + 
            (2.0/3.0) * (u2[i].E - dt * dF3[i].E);
    }
}
/*void RungeKutta3FDS(vector<Conserved>& U, double dt, double dx, int diff_calc_type)
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

    CalcFluxFDS(U, dx, f);
    for (int i = 0; i < n; ++i)
    {
        U1[i].rho = U0[i].rho - dt * f[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * f[i].rho_u;
        U1[i].E = U0[i].E - dt * f[i].E;
    }
    ApplyBoundaryConditions(U1);

    CalcFluxFDS(U1, dx, f);
    for (int i = 0; i < n; ++i)
    {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * f[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * f[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * f[i].E);
    }
    ApplyBoundaryConditions(U2);

    CalcFluxFDS(U2, dx, f);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * f[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * f[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * f[i].E);
    }
    ApplyBoundaryConditions(U);
}*/
void RungeKutta3FVS(vector<Conserved>& U, double dt, double dx, int diff_calc_type)
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
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, diff_calc_type);
    for (int i = 0; i < n; ++i)
    {
        U1[i].rho = U0[i].rho - dt * diff[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * diff[i].rho_u;
        U1[i].E = U0[i].E - dt * diff[i].E;
    }
    ApplyBoundaryConditions(U1);
    for (int i = 0; i < n; ++i)
    {
        f_left[i] = {0.0, 0.0, 0.0};
        f_right[i] = {0.0, 0.0, 0.0};
        diff[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFVS(U1, f_left, f_right);
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, diff_calc_type);
    for (int i = 0; i < n; ++i)
    {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * diff[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * diff[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * diff[i].E);
    }
    ApplyBoundaryConditions(U2);
    for (int i = 0; i < n; ++i)
    {
        f_left[i] = {0.0, 0.0, 0.0};
        f_right[i] = {0.0, 0.0, 0.0};
        diff[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFVS(U2, f_left, f_right);
    CalcDiff(f_left, f_right, temp_left, temp_right, diff, dx, diff_calc_type);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * diff[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * diff[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * diff[i].E);
    }
    ApplyBoundaryConditions(U);
}
void RungeKutta3(vector<Conserved>& U, double dt, double dx, int flux_process_type, int diff_calc_type, int n)
{
    if (flux_process_type == 0)
    {
        // Use FDS flux processing
        RungeKutta3FDS(U, dt, n, dx);
    }
    else if (flux_process_type == 1)
    {
        // Use FVS flux processing
        RungeKutta3FVS(U, dt, dx, diff_calc_type);
    }
    else
    {
        // Handle other flux processing types if needed
        cout << "Unsupported flux process type: " << flux_process_type << endl;
    }
}