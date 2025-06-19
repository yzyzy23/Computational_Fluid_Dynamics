#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// Initialize constants
const double dx = 0.005;
const int N = static_cast<int>((1.0 - (-1.0)) / dx) + 1; // Number of grid points
const double end_time = 0.14; // End time for the simulation
const double gamma_val = 1.4;
const double R_val = 286.9;
const double CFL = 0.2;
const double x_min = -1.0;
const double x_max = 1.0;
const double eps = 1e-5; // Small value to avoid division by zero

struct State {double rho, u, p;};

struct Conserved
{double rho, rho_u, E;};

vector<Conserved> InitCondition(double dx)
{
    vector<Conserved> U(N);
    for (int i = 0; i < N; ++i)
    {
        double x = x_min + i * dx;
        if (x < 0)
        {
            U[i] = {1.0, 0.0, 1.0};
        }
        else
        {
            U[i] = {0.125, 0.0, 0.1};
        }
    }
    return U;
}

void ConservedToState(const vector<Conserved>& U, vector<double>& rho, vector<double>& u, vector<double>& p)
{
    for (int i = 0; i < N; ++i)
    {
        rho[i] = U[i].rho;
        u[i] = U[i].rho_u / U[i].rho;
        p[i] = (gamma_val - 1) * (U[i].E - 0.5 * U[i].rho * u[i] * u[i]);
    }
}

Conserved flux_roe(const Conserved& UL, const Conserved& UR)
{
    double rho_left = UL.rho, rho_right = UR.rho;
    double u_left = UL.rho_u / rho_left, u_right = UR.rho_u / rho_right;
    double p_left = (gamma_val - 1) * (UL.E - 0.5 * rho_left * u_left * u_left);
    double p_right = (gamma_val - 1) * (UR.E - 0.5 * rho_right * u_right * u_right);

    double H_left = (gamma_val / (gamma_val - 1)) * p_left / rho_left + 0.5 * u_left * u_left;
    double H_right = (gamma_val / (gamma_val - 1)) * p_right / rho_right + 0.5 * u_right * u_right;

    double sqrt_rho_left = sqrt(rho_left), sqrt_rho_right = sqrt(rho_right);
    double u_roe = (sqrt_rho_left * u_left + sqrt_rho_right * u_right) / (sqrt_rho_left + sqrt_rho_right);
    double H_roe = (sqrt_rho_left * H_left + sqrt_rho_right * H_right) / (sqrt_rho_left + sqrt_rho_right);
    double c_roe = sqrt((gamma_val - 1) * (H_roe - 0.5 * u_roe * u_roe));

    double lambda1 = u_roe - c_roe;
    double lambda2 = u_roe;
    double lambda3 = u_roe + c_roe;

    double eps = 1e-6;
    lambda1 = (abs(lambda1) < eps) ? (lambda1 * lambda1 + eps * eps) / (2 * eps) : lambda1;
    lambda3 = (abs(lambda3) < eps) ? (lambda3 * lambda3 + eps * eps) / (2 * eps) : lambda3;

    double drho = UR.rho - UL.rho;
    double drho_u = UR.rho_u - UL.rho_u;
    double dE = UR.E - UL.E;

    Conserved flux;
    flux.rho = 0.5 * (UL.rho * u_left + UR.rho * u_right) - 0.5 * (abs(lambda1) * drho);
    flux.rho_u = 0.5 * (UL.rho_u * u_left + UR.rho_u * u_right + p_left + p_right) - 0.5 * (abs(lambda2) * drho_u);
    flux.E = 0.5 * (u_left * (UL.E + p_left) + u_right * (UR.E + p_right)) - 0.5 * (abs(lambda3) * dE);

    return flux;
}

void CalcWENO(const vector<double>& p, const vector<double>& u, const vector<Conserved>& U,
    vector<vector<double>>& F_left, vector<vector<double>>& F_right, vector<vector<double>>& diff_left, vector<vector<double>>& diff_right)
{
    vector<vector<double>> lambda(3, vector<double>(N, 0.0));
    vector<vector<double>> lambda_left(3, vector<double>(N, 0.0));
    vector<vector<double>> lambda_right(3, vector<double>(N, 0.0));
    vector<double> C(N, 0.0);
    double eps = 1e-3;

    for (int i = 0; i < N; ++i)
    {
        C[i] = sqrt(abs(gamma_val * p[i] / U[i].rho));
        lambda[0][i] = u[i];
        lambda[1][i] = u[i] - C[i];
        lambda[2][i] = u[i] + C[i];
        for (int k = 0; k < 3; ++k)
        {
            lambda_left[k][i] = (abs(lambda[k][i]) + lambda[k][i]) / 2;
            lambda_right[k][i] = lambda[k][i] - lambda_left[k][i];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        F_left[0][i] = U[i].rho / (2 * gamma_val) * (2 * (gamma_val - 1) * lambda_left[0][i] + lambda_left[1][i] + lambda_left[2][i]);
        F_left[1][i] = U[i].rho / (2 * gamma_val) * (2 * (gamma_val - 1) * lambda_left[0][i] * u[i] + lambda_left[1][i] * (u[i] - C[i]) + lambda_left[2][i] * (u[i] + C[i]));
        F_left[2][i] = U[i].rho / (2 * gamma_val) * ((gamma_val - 1) * lambda_left[0][i] * u[i] * u[i] + 0.5 * lambda_left[1][i] * (u[i] - C[i]) * (u[i] - C[i]) + 0.5 * lambda_left[2][i] * (u[i] + C[i]) * (u[i] + C[i]) + (3 - gamma_val) / (2 * gamma_val - 2) * (lambda_left[1][i] + lambda_left[2][i]) * C[i] * C[i]);

        F_right[0][i] = U[i].rho / (2 * gamma_val) * (2 * (gamma_val - 1) * lambda_right[0][i] + lambda_right[1][i] + lambda_right[2][i]);
        F_right[1][i] = U[i].rho / (2 * gamma_val) * (2 * (gamma_val - 1) * lambda_right[0][i] * u[i] + lambda_right[1][i] * (u[i] - C[i]) + lambda_right[2][i] * (u[i] + C[i]));
        F_right[2][i] = U[i].rho / (2 * gamma_val) * ((gamma_val - 1) * lambda_right[0][i] * u[i] * u[i] + 0.5 * lambda_right[1][i] * (u[i] - C[i]) * (u[i] - C[i]) + 0.5 * lambda_right[2][i] * (u[i] + C[i]) * (u[i] + C[i]) + (3 - gamma_val) / (2 * gamma_val - 2) * (lambda_right[1][i] + lambda_right[2][i]) * C[i] * C[i]);
    }

    for (int i = 0; i < N - 4; ++i)
    {
        double IS1_left = 0.0, IS2_left = 0.0, IS3_left = 0.0;
        double IS1_right = 0.0, IS2_right = 0.0, IS3_right = 0.0;

        for (int k = 0; k < 3; ++k)
        {
            IS1_left += 1.0 / 4.0 * pow(F_left[k][i] - 4 * F_left[k][i + 1] + 3 * F_left[k][i + 2], 2) + 13.0 / 12.0 * pow(F_left[k][i] - 2 * F_left[k][i + 1] + F_left[k][i + 2], 2);
            IS2_left += 1.0 / 4.0 * pow(F_left[k][i + 1] - F_left[k][i + 3], 2) + 13.0 / 12.0 * pow(F_left[k][i + 1] - 2 * F_left[k][i + 2] + F_left[k][i + 3], 2);
            IS3_left += 1.0 / 4.0 * pow(3 * F_left[k][i + 2] - 4 * F_left[k][i + 3] + F_left[k][i + 4], 2) + 13.0 / 12.0 * pow(F_left[k][i + 2] - 2 * F_left[k][i + 3] + F_left[k][i + 4], 2);

            IS1_right += 1.0 / 4.0 * pow(F_right[k][i + 4] - 4 * F_right[k][i + 3] + 3 * F_right[k][i + 2], 2) + 13.0 / 12.0 * pow(F_right[k][i + 4] - 2 * F_right[k][i + 3] + F_right[k][i + 2], 2);
            IS2_right += 1.0 / 4.0 * pow(F_right[k][i + 3] - F_right[k][i + 1], 2) + 13.0 / 12.0 * pow(F_right[k][i + 3] - 2 * F_right[k][i + 2] + F_right[k][i + 1], 2);
            IS3_right += 1.0 / 4.0 * pow(3 * F_right[k][i + 2] - 4 * F_right[k][i + 1] + F_right[k][i], 2) + 13.0 / 12.0 * pow(F_right[k][i + 2] - 2 * F_right[k][i + 1] + F_right[k][i], 2);
        }

        double omega1_left = 0.1 / pow(IS1_left + eps, 2) / (0.1 / pow(IS1_left + eps, 2) + 0.6 / pow(IS2_left + eps, 2) + 0.3 / pow(IS3_left + eps, 2));
        double omega2_left = 0.6 / pow(IS2_left + eps, 2) / (0.1 / pow(IS1_left + eps, 2) + 0.6 / pow(IS2_left + eps, 2) + 0.3 / pow(IS3_left + eps, 2));
        double omega3_left = 1 - omega1_left - omega2_left;

        double omega1_right = 0.1 / pow(IS1_right + eps, 2) / (0.1 / pow(IS1_right + eps, 2) + 0.6 / pow(IS2_right + eps, 2) + 0.3 / pow(IS3_right + eps, 2));
        double omega2_right = 0.6 / pow(IS2_right + eps, 2) / (0.1 / pow(IS1_right + eps, 2) + 0.6 / pow(IS2_right + eps, 2) + 0.3 / pow(IS3_right + eps, 2));
        double omega3_right = 1 - omega1_right - omega2_right;

        for (int k = 0; k < 3; ++k)
        {
            diff_left[k][i + 2] = omega1_left * (2 * F_left[k][i] / 6 - 7 * F_left[k][i + 1] / 6 + 11 * F_left[k][i + 2] / 6) +
                    omega2_left * (-F_left[k][i + 1] / 6 + 5 * F_left[k][i + 2] / 6 + 2 * F_left[k][i + 3] / 6) +
                    omega3_left * (2 * F_left[k][i + 2] / 6 + 5 * F_left[k][i + 3] / 6 - F_left[k][i + 4] / 6);

            diff_right[k][i + 2] = omega1_right * (2 * F_right[k][i + 4] / 6 - 7 * F_right[k][i + 3] / 6 + 11 * F_right[k][i + 2] / 6) +
                    omega2_right * (-F_right[k][i + 3] / 6 + 5 * F_right[k][i + 2] / 6 + 2 * F_right[k][i + 1] / 6) +
                    omega3_right * (2 * F_right[k][i + 2] / 6 + 5 * F_right[k][i + 1] / 6 - F_right[k][i] / 6);
        }
    }
}
void RungeKutta2(vector<Conserved>& U, vector<vector<double>>& diff_left, vector<vector<double>>& diff_right, const vector<double>& p, const vector<double>& u, vector<vector<double>>& F_left, vector<vector<double>>& F_right, double dt)
{
    CalcWENO(p, u, U, F_left, F_right, diff_left, diff_right);
    for (int i = 2; i < N - 4; ++i)
    {
        U[i + 1].rho -= (dt / dx) * (diff_left[0][i + 1] - diff_left[0][i] + diff_right[0][i + 2] - diff_right[0][i + 1]);
        U[i + 1].rho_u -= (dt / dx) * (diff_left[1][i + 1] - diff_left[1][i] + diff_right[1][i + 2] - diff_right[1][i + 1]);
        U[i + 1].E -= (dt / dx) * (diff_left[2][i + 1] - diff_left[2][i] + diff_right[2][i + 2] - diff_right[2][i + 1]);
    }
}
double GVCLimiter(double r)
{
    // GVC limiter
    return abs(r) > 1.0 ? 1.0 : r;
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

        // Apply GVC limiter
        double limiter_rho_left = GVCLimiter(r_left_rho);
        double limiter_rho_u_left = GVCLimiter(r_left_rho_u);
        double limiter_E_left = GVCLimiter(r_left_E);
        double limiter_rho_right = GVCLimiter(r_right_rho);
        double limiter_rho_u_right = GVCLimiter(r_right_rho_u);
        double limiter_E_right = GVCLimiter(r_right_E);

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
void CalcFluxFVS(const vector<Conserved>& U, vector<Conserved>& F_left, vector<Conserved>& F_right)
{
    // Use Steger-Warming flux vector splitting
    int n = U.size();
    // Initialize eigenvalues and output
    vector<double> lambda(3), lambda_left(3), lambda_right(3);
    for (int i = 0; i < n; ++i)
    {
        F_left[i] = {0.0, 0.0, 0.0};
        F_right[i] = {0.0, 0.0, 0.0};
    }

    for (int i = 0; i < n; ++i)
    {
        State state = {U[i].rho, U[i].rho_u / U[i].rho, (gamma_val - 1.0) * (U[i].E / U[i].rho - 0.5 * (U[i].rho_u / U[i].rho) * (U[i].rho_u / U[i].rho))};
        double rho = state.rho;
        double u = state.u;
        double T = (gamma_val - 1.0) * (U[i].E / rho - 0.5 * u * u) / R_val;
        double p = rho * R_val * T;
        double c = sqrt(gamma_val * p / rho);
        double local_eps = 1e-3;

        lambda[0] = u;
        lambda[1] = u - c;
        lambda[2] = u + c;

        for (int j = 0; j < 3; ++j)
        {
            lambda_left[j] = (lambda[j] + sqrt(lambda[j] * lambda[j] + local_eps * local_eps)) / 2.0;
            lambda_right[j] = (lambda[j] - sqrt(lambda[j] * lambda[j] + local_eps * local_eps)) / 2.0;
        }

        double omega_left = ((3 - gamma_val) * (lambda_left[1] + lambda_left[2]) * c * c) / (2.0 * (gamma_val - 1.0));
        F_left[i].rho = (rho / (2 * gamma_val)) * ((2 * (gamma_val - 1) * lambda_left[0]) + lambda_left[1] + lambda_left[2]);
        F_left[i].rho_u = (rho / (2 * gamma_val)) * ((2 * (gamma_val - 1) * lambda_left[0] * u) + lambda_left[1] * (u - c) + lambda_left[2] * (u + c));
        F_left[i].E = (rho / (2 * gamma_val)) * ((gamma_val - 1) * lambda_left[0] * u * u + 0.5 * lambda_left[1] * (u - c) * (u - c) + 0.5 * lambda_left[2] * (u + c) * (u + c) + omega_left);

        double omega_right = ((3 - gamma_val) * (lambda_right[1] + lambda_right[2]) * c * c) / (2.0 * (gamma_val - 1.0));
        F_right[i].rho = (rho / (2 * gamma_val)) * ((2 * (gamma_val - 1) * lambda_right[0]) + lambda_right[1] + lambda_right[2]);
        F_right[i].rho_u = (rho / (2 * gamma_val)) * ((2 * (gamma_val - 1) * lambda_right[0] * u) + lambda_right[1] * (u - c) + lambda_right[2] * (u + c));
        F_right[i].E = (rho / (2 * gamma_val)) * ((gamma_val - 1) * lambda_right[0] * u * u + 0.5 * lambda_right[1] * (u - c) * (u - c) + 0.5 * lambda_right[2] * (u + c) * (u + c) + omega_right);
    }
}
void RungeKutta3(vector<Conserved>& U, double dt, double dx)
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
    CalcDiffTVD(f_left, f_right, temp_left, temp_right, diff, dx);
    for (int i = 0; i < n; ++i)
    {
        U1[i].rho = U0[i].rho - dt * diff[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * diff[i].rho_u;
        U1[i].E = U0[i].E - dt * diff[i].E;
    }

    for (int i = 0; i < n; ++i)
    {
        f_left[i] = {0.0, 0.0, 0.0};
        f_right[i] = {0.0, 0.0, 0.0};
        diff[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFVS(U1, f_left, f_right);
    CalcDiffTVD(f_left, f_right, temp_left, temp_right, diff, dx);
    for (int i = 0; i < n; ++i)
    {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * diff[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * diff[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * diff[i].E);
    }

    for (int i = 0; i < n; ++i)
    {
        f_left[i] = {0.0, 0.0, 0.0};
        f_right[i] = {0.0, 0.0, 0.0};
        diff[i] = {0.0, 0.0, 0.0};
    }

    CalcFluxFVS(U2, f_left, f_right);
    CalcDiffTVD(f_left, f_right, temp_left, temp_right, diff, dx);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * diff[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * diff[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * diff[i].E);
    }
}

int main()
{
    double dx = (x_max - x_min) / (N - 1);
    vector<Conserved> U = InitCondition(dx);
    vector<vector<double>> F_left(3, vector<double>(N, 0.0)), F_right(3, vector<double>(N, 0.0));
    vector<vector<double>> diff_left(3, vector<double>(N, 0.0)), diff_right(3, vector<double>(N, 0.0));

    vector<double> rho(N), u(N), p(N);
    double t = 0.0;

    while (t < end_time)
    {
        ConservedToState(U, rho, u, p);
        vector<double> a(N);
        for (int i = 0; i < N; ++i) {
            a[i] = sqrt(gamma_val * p[i] / rho[i]);
        }
        double max_speed = 0.0;
        for (int i = 0; i < N; ++i) 
        {
            max_speed = max(max_speed, abs(u[i] + a[i]));
        }
        double dt = CFL * dx / max_speed;
        dt = min(dt, end_time - t);

        RungeKutta2(U, diff_left, diff_right, p, u, F_left, F_right, dt);

        t += dt;
    }

    FILE* file = fopen("SodResult.txt", "w");
    if (file == NULL)
    {
        cout << "Error opening file: SodResult.txt!" << endl;
        return -1;
    }
    ConservedToState(U, rho, u, p);
    for (int i = 0; i < N; ++i)
    {

        fprintf(file, "%lf %lf %lf %lf\n", x_min + i * dx, rho[i], u[i], p[i]);
    }
    fclose(file);

    return 0;
}