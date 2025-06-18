#include <iostream>
#include <vector>
#include "DataStructure.h"

using namespace std;

double VanLeerLimiter(double r)
{
    // Van Leer limiter
    return (r + abs(r)) / (1.0 + abs(r));
}
void CalcDiffTVD(const vector<Conserved>& U_left, const vector<Conserved>& U_right, vector<Conserved>& new_u_left, vector<Conserved>& new_u_right, vector<Conserved>& diff, double dx)
{
    // Use TVD to calculate differences with Van Leer limiter
    int n = U_left.size();
    vector<Conserved> diff_left(n), diff_right(n);

    for (int i = 1; i < n - 2; ++i)
    {
        double r_left_rho = (U_left[i].rho - U_left[i - 1].rho) / (U_left[i + 1].rho - U_left[i].rho + eps);
        double r_left_rho_u = (U_left[i].rho_u - U_left[i - 1].rho_u) / (U_left[i + 1].rho_u - U_left[i].rho_u + eps);
        double r_left_E = (U_left[i].E - U_left[i - 1].E) / (U_left[i + 1].E - U_left[i].E + eps);
        double r_right_rho = (U_right[i + 2].rho - U_right[i + 1].rho) / (U_right[i + 1].rho - U_right[i].rho + eps);
        double r_right_rho_u = (U_right[i + 2].rho_u - U_right[i + 1].rho_u) / (U_right[i + 1].rho_u - U_right[i].rho_u + eps);
        double r_right_E = (U_right[i + 2].E - U_right[i + 1].E) / (U_right[i + 1].E - U_right[i].E + eps);

        // Apply Van Leer limiter
        double limiter_rho_left = VanLeerLimiter(r_left_rho);
        double limiter_rho_u_left = VanLeerLimiter(r_left_rho_u);
        double limiter_E_left = VanLeerLimiter(r_left_E);
        double limiter_rho_right = VanLeerLimiter(r_right_rho);
        double limiter_rho_u_right = VanLeerLimiter(r_right_rho_u);
        double limiter_E_right = VanLeerLimiter(r_right_E);

        new_u_left[i].rho = U_left[i].rho + 0.5 * (limiter_rho_left * (U_left[i + 1].rho - U_left[i].rho));
        new_u_left[i].rho_u = U_left[i].rho_u + 0.5 * (limiter_rho_u_left * (U_left[i + 1].rho_u - U_left[i].rho_u));
        new_u_left[i].E = U_left[i].E + 0.5 * (limiter_E_left * (U_left[i + 1].E - U_left[i].E));
        new_u_right[i].rho = U_right[i + 1].rho - 0.5 * (limiter_rho_right * (U_right[i + 1].rho - U_right[i].rho));
        new_u_right[i].rho_u = U_right[i + 1].rho_u - 0.5 * (limiter_rho_u_right * (U_right[i + 1].rho_u - U_right[i].rho_u));
        new_u_right[i].E = U_right[i + 1].E - 0.5 * (limiter_E_right * (U_right[i + 1].E - U_right[i].E));
    }

    for (int i = 2; i < n - 2; ++i)
    {
        diff[i].rho = (new_u_left[i].rho - new_u_left[i - 1].rho + new_u_right[i].rho - new_u_right[i - 1].rho) / dx;
        diff[i].rho_u = (new_u_left[i].rho_u - new_u_left[i - 1].rho_u + new_u_right[i].rho_u - new_u_right[i - 1].rho_u) / dx;
        diff[i].E = (new_u_left[i].E - new_u_left[i - 1].E + new_u_right[i].E - new_u_right[i - 1].E) / dx;
    }
}
void CalcDiffWENO(const vector<Conserved>& U_left, const vector<Conserved>& U_right, vector<Conserved>& new_u_left, vector<Conserved>& new_u_right, vector<Conserved>& diff, double dx)
{
    int n = U_left.size();
    double local_eps = 1e-6; // eps value for WENO
    vector<double> C = {0.1, 0.6, 0.3}; // WENO coefficients

    for (int i = 2; i < n - 2; ++i)
    {
        // for left state
        double IS1_rho = 0.25 * pow((U_left[i - 2].rho - 4 * U_left[i - 1].rho + 3 * U_left[i].rho), 2) + (13.0 / 12) * pow((U_left[i - 2].rho - 2 * U_left[i - 1].rho + U_left[i].rho), 2);
        double IS1_rho_u = 0.25 * pow((U_left[i - 2].rho_u - 4 * U_left[i - 1].rho_u + 3 * U_left[i].rho_u), 2) + (13.0 / 12) * pow((U_left[i - 2].rho_u - 2 * U_left[i - 1].rho_u + U_left[i].rho_u), 2);
        double IS1_E = 0.25 * pow((U_left[i - 2].E - 4 * U_left[i - 1].E + 3 * U_left[i].E), 2) + (13.0 / 12) * pow((U_left[i - 2].E - 2 * U_left[i - 1].E + U_left[i].E), 2);
        double IS2_rho = 0.25 * pow((U_left[i - 1].rho - U_left[i + 1].rho), 2) + (13.0 / 12) * pow((U_left[i - 1].rho - 2 * U_left[i].rho + U_left[i + 1].rho), 2);
        double IS2_rho_u = 0.25 * pow((U_left[i - 1].rho_u - U_left[i + 1].rho_u), 2) + (13.0 / 12) * pow((U_left[i - 1].rho_u - 2 * U_left[i].rho_u + U_left[i + 1].rho_u), 2);
        double IS2_E = 0.25 * pow((U_left[i - 1].E - U_left[i + 1].E), 2) + (13.0 / 12) * pow((U_left[i - 1].E - 2 * U_left[i].E + U_left[i + 1].E), 2);
        double IS3_rho = 0.25 * pow((3 * U_left[i].rho - 4 * U_left[i + 1].rho + U_left[i + 2].rho), 2) + (13.0 / 12) * pow((U_left[i].rho - 2 * U_left[i + 1].rho + U_left[i + 2].rho), 2);
        double IS3_rho_u = 0.25 * pow((3 * U_left[i].rho_u - 4 * U_left[i + 1].rho_u + U_left[i + 2].rho_u), 2) + (13.0 / 12) * pow((U_left[i].rho_u - 2 * U_left[i + 1].rho_u + U_left[i + 2].rho_u), 2);
        double IS3_E = 0.25 * pow((3 * U_left[i].E - 4 * U_left[i + 1].E + U_left[i + 2].E), 2) + (13.0 / 12) * pow((U_left[i].E - 2 * U_left[i + 1].E + U_left[i + 2].E), 2);
        vector<double> IS_rho = {IS1_rho, IS2_rho, IS3_rho};
        vector<double> IS_rho_u = {IS1_rho_u, IS2_rho_u, IS3_rho_u};
        vector<double> IS_E = {IS1_E, IS2_E, IS3_E};
        vector<double> alpha_rho = {0.0, 0.0, 0.0};
        vector<double> alpha_rho_u = {0.0, 0.0, 0.0};
        vector<double> alpha_E = {0.0, 0.0, 0.0};
        vector<double> omega_rho = {0.0, 0.0, 0.0};
        vector<double> omega_rho_u = {0.0, 0.0, 0.0};
        vector<double> omega_E = {0.0, 0.0, 0.0};

        for (int i = 0; i < 3; ++i)
        {
            alpha_rho[i] = C[i] / pow(local_eps + IS_rho[i], 2);
            alpha_rho_u[i] = C[i] / pow(local_eps + IS_rho_u[i], 2);
            alpha_E[i] = C[i] / pow(local_eps + IS_E[i], 2);
        }
        for (int i = 0; i < 3; ++i)
        {
            omega_rho[i] = alpha_rho[i] / (alpha_rho[0] + alpha_rho[1] + alpha_rho[2]);
            omega_rho_u[i] = alpha_rho_u[i] / (alpha_rho_u[0] + alpha_rho_u[1] + alpha_rho_u[2]);
            omega_E[i] = alpha_E[i] / (alpha_E[0] + alpha_E[1] + alpha_E[2]);
        }

        vector<Conserved> f(3);
        f[0].rho = (1.0 / 3) * U_left[i - 2].rho - (7.0 / 6) * U_left[i - 1].rho + (11.0 / 6) * U_left[i].rho;
        f[0].rho_u = (1.0 / 3) * U_left[i - 2].rho_u - (7.0 / 6) * U_left[i - 1].rho_u + (11.0 / 6) * U_left[i].rho_u;
        f[0].E = (1.0 / 3) * U_left[i - 2].E - (7.0 / 6) * U_left[i - 1].E + (11.0 / 6) * U_left[i].E;
        f[1].rho = (-1.0 / 6) * U_left[i - 1].rho + (5.0 / 6) * U_left[i].rho + (1.0 / 3) * U_left[i + 1].rho;
        f[1].rho_u = (-1.0 / 6) * U_left[i - 1].rho_u + (5.0 / 6) * U_left[i].rho_u + (1.0 / 3) * U_left[i + 1].rho_u;
        f[1].E = (-1.0 / 6) * U_left[i - 1].E + (5.0 / 6) * U_left[i].E + (1.0 / 3) * U_left[i + 1].E;
        f[2].rho = (1.0 / 3) * U_left[i].rho + (5.0 / 6) * U_left[i + 1].rho - (1.0 / 6) * U_left[i + 2].rho;
        f[2].rho_u = (1.0 / 3) * U_left[i].rho_u + (5.0 / 6) * U_left[i + 1].rho_u - (1.0 / 6) * U_left[i + 2].rho_u;
        f[2].E = (1.0 / 3) * U_left[i].E + (5.0 / 6) * U_left[i + 1].E - (1.0 / 6) * U_left[i + 2].E;

        for (int j = 0; j < 3; ++j)
        {
            new_u_left[i].rho += omega_rho[j] * f[j].rho;
            new_u_left[i].rho_u += omega_rho_u[j] * f[j].rho_u;
            new_u_left[i].E += omega_E[j] * f[j].E;
        }

        // for right state
        IS1_rho = 0.25 * pow((U_right[i - 2].rho - 4 * U_right[i - 1].rho + 3 * U_right[i].rho), 2) + (13.0 / 12) * pow((U_right[i - 2].rho - 2 * U_right[i - 1].rho + U_right[i].rho), 2);
        IS1_rho_u = 0.25 * pow((U_right[i - 2].rho_u - 4 * U_right[i - 1].rho_u + 3 * U_right[i].rho_u), 2) + (13.0 / 12) * pow((U_right[i - 2].rho_u - 2 * U_right[i - 1].rho_u + U_right[i].rho_u), 2);
        IS1_E = 0.25 * pow((U_right[i - 2].E - 4 * U_right[i - 1].E + 3 * U_right[i].E), 2) + (13.0 / 12) * pow((U_right[i - 2].E - 2 * U_right[i - 1].E + U_right[i].E), 2);
        IS2_rho = 0.25 * pow((U_right[i - 1].rho - U_right[i + 1].rho), 2) + (13.0 / 12) * pow((U_right[i - 1].rho - 2 * U_right[i].rho + U_right[i + 1].rho), 2);
        IS2_rho_u = 0.25 * pow((U_right[i - 1].rho_u - U_right[i + 1].rho_u), 2) + (13.0 / 12) * pow((U_right[i - 1].rho_u - 2 * U_right[i].rho_u + U_right[i + 1].rho_u), 2);
        IS2_E = 0.25 * pow((U_right[i - 1].E - U_right[i + 1].E), 2) + (13.0 / 12) * pow((U_right[i - 1].E - 2 * U_right[i].E + U_right[i + 1].E), 2);
        IS3_rho = 0.25 * pow((3 * U_right[i].rho - 4 * U_right[i + 1].rho + U_right[i + 2].rho), 2) + (13.0 / 12) * pow((U_right[i].rho - 2 * U_right[i + 1].rho + U_right[i + 2].rho), 2);
        IS3_rho_u = 0.25 * pow((3 * U_right[i].rho_u - 4 * U_right[i + 1].rho_u + U_right[i + 2].rho_u), 2) + (13.0 / 12) * pow((U_right[i].rho_u - 2 * U_right[i + 1].rho_u + U_right[i + 2].rho_u), 2);
        IS3_E = 0.25 * pow((3 * U_right[i].E - 4 * U_right[i + 1].E + U_right[i + 2].E), 2) + (13.0 / 12) * pow((U_right[i].E - 2 * U_right[i + 1].E + U_right[i + 2].E), 2);
        IS_rho = {IS1_rho, IS2_rho, IS3_rho};
        IS_rho_u = {IS1_rho_u, IS2_rho_u, IS3_rho_u};
        IS_E = {IS1_E, IS2_E, IS3_E};
        alpha_rho = {0.0, 0.0, 0.0};
        alpha_rho_u = {0.0, 0.0, 0.0};
        alpha_E = {0.0, 0.0, 0.0};
        omega_rho = {0.0, 0.0, 0.0};
        omega_rho_u = {0.0, 0.0, 0.0};
        omega_E = {0.0, 0.0, 0.0};

        for (int i = 0; i < 3; ++i)
        {
            alpha_rho[i] = C[i] / pow(local_eps + IS_rho[i], 2);
            alpha_rho_u[i] = C[i] / pow(local_eps + IS_rho_u[i], 2);
            alpha_E[i] = C[i] / pow(local_eps + IS_E[i], 2);
        }
        for (int i = 0; i < 3; ++i)
        {
            omega_rho[i] = alpha_rho[i] / (alpha_rho[0] + alpha_rho[1] + alpha_rho[2]);
            omega_rho_u[i] = alpha_rho_u[i] / (alpha_rho_u[0] + alpha_rho_u[1] + alpha_rho_u[2]);
            omega_E[i] = alpha_E[i] / (alpha_E[0] + alpha_E[1] + alpha_E[2]);
        }
        f[0].rho = (1.0 / 3) * U_right[i + 2].rho - (7.0 / 6) * U_right[i + 1].rho + (11.0 / 6) * U_right[i].rho;
        f[0].rho_u = (1.0 / 3) * U_right[i + 2].rho_u - (7.0 / 6) * U_right[i + 1].rho_u + (11.0 / 6) * U_right[i].rho_u;
        f[0].E = (1.0 / 3) * U_right[i + 2].E - (7.0 / 6) * U_right[i + 1].E + (11.0 / 6) * U_right[i].E;
        f[1].rho = (-1.0 / 6) * U_right[i + 1].rho + (5.0 / 6) * U_right[i].rho + (1.0 / 3) * U_right[i - 1].rho;
        f[1].rho_u = (-1.0 / 6) * U_right[i + 1].rho_u + (5.0 / 6) * U_right[i].rho_u + (1.0 / 3) * U_right[i - 1].rho_u;
        f[1].E = (-1.0 / 6) * U_right[i + 1].E + (5.0 / 6) * U_right[i].E + (1.0 / 3) * U_right[i - 1].E;
        f[2].rho = (1.0 / 3) * U_right[i].rho + (5.0 / 6) * U_right[i - 1].rho - (1.0 / 6) * U_right[i - 2].rho;
        f[2].rho_u = (1.0 / 3) * U_right[i].rho_u + (5.0 / 6) * U_right[i - 1].rho_u - (1.0 / 6) * U_right[i - 2].rho_u;
        f[2].E = (1.0 / 3) * U_right[i].E + (5.0 / 6) * U_right[i - 1].E - (1.0 / 6) * U_right[i - 2].E;

        for (int j = 0; j < 3; ++j)
        {
            new_u_right[i].rho += omega_rho[j] * f[j].rho;
            new_u_right[i].rho_u += omega_rho_u[j] * f[j].rho_u;
            new_u_right[i].E += omega_E[j] * f[j].E;
        }
    }
    for (int i = 3; i < n - 3; ++i)
    {
        diff[i].rho = (new_u_left[i].rho - new_u_left[i - 1].rho + new_u_right[i + 1].rho - new_u_right[i].rho) / dx;
        diff[i].rho_u = (new_u_left[i].rho_u - new_u_left[i - 1].rho_u + new_u_right[i + 1].rho_u - new_u_right[i].rho_u) / dx;
        diff[i].E = (new_u_left[i].E - new_u_left[i - 1].E + new_u_right[i + 1].E - new_u_right[i].E) / dx;
    }
}
void CalcDiff(const vector<Conserved>& U_left, const vector<Conserved>& U_right, vector<Conserved>& new_u_left, vector<Conserved>& new_u_right, vector<Conserved>& diff, double dx, int diff_calc_type)
{
    if (diff_calc_type == 0)
    {
        // Use TVD method
        CalcDiffTVD(U_left, U_right, new_u_left, new_u_right, diff, dx);
    }
    else if (diff_calc_type == 1)
    {
        // Use WENO method
        CalcDiffWENO(U_left, U_right, new_u_left, new_u_right, diff, dx);
    }
    else
    {
        cerr << "Invalid diff_calc_type: " << diff_calc_type << endl;
        exit(EXIT_FAILURE);
    }
}