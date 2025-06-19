#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const int N = 100;
const double gamma_val = 1.4;
const double end_time = 0.2;
const double CFL = 0.2;
const double xmin = -1.0, xmax = 1.0;

struct Conserved
{
    double rho;
    double rho_u;
    double E;
};

vector<Conserved> InitCondition(double dx)
{
    vector<Conserved> U(N);
    for (int i = 0; i < N; ++i)
    {
        double x = xmin + i * dx;
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

void MinModLimiter(const vector<double>& v, vector<double>& v_left, vector<double>& v_right)
{
    int n = v.size();
    v_left[0] = v[0];
    v_right[0] = v[0];
    v_left[n - 1] = v[n - 1];
    v_right[n - 1] = v[n - 1];

    for (int i = 1; i < n - 1; ++i)
    {
        double slope_l = v[i] - v[i - 1];
        double slope_r = v[i + 1] - v[i];
        double slope = (slope_l * slope_r <= 0) ? 0.0 : (slope_l > 0 ? 1 : -1) * min(abs(slope_l), abs(slope_r));
        v_left[i] = v[i] - 0.5 * slope;
        v_right[i] = v[i] + 0.5 * slope;
    }
}
void GVCLimiter(const vector<double>& v, vector<double>& v_left, vector<double>& v_right)
{
    int n = v.size();
    v_left[0] = v[0];
    v_right[0] = v[0];
    v_left[n - 1] = v[n - 1];
    v_right[n - 1] = v[n - 1];

    for (int i = 1; i < n - 1; ++i)
    {
        double slope_l = v[i] - v[i - 1];
        double slope_r = v[i + 1] - v[i];
        double r = slope_l / slope_r;
        double phi = abs(r) > 1.0 ? 1.0 : abs(r);
        v_left[i] = v[i] - 0.5 * phi * slope_r;
        v_right[i] = v[i] + 0.5 * phi * slope_r;
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

vector<Conserved> compute_flux_difference(const vector<Conserved>& U, double dx)
{
    vector<Conserved> U_left(N), U_right(N);
    vector<Conserved> F(N + 1), dF(N);

    for (int i = 0; i < N; ++i)
    {
        vector<double> v(N);
        for (int j = 0; j < N; ++j) v[j] = U[j].rho;
        vector<double> v_left(N), v_right(N);
        MinModLimiter(v, v_left, v_right);
        U_left[i].rho = v_left[i];
        U_right[i].rho = v_right[i];
    }
    for (int i = 0; i < N; ++i)
    {
        vector<double> v(N);
        for (int j = 0; j < N; ++j) v[j] = U[j].rho_u;
        vector<double> v_left(N), v_right(N);
        MinModLimiter(v, v_left, v_right);
        U_left[i].rho_u = v_left[i];
        U_right[i].rho_u = v_right[i];
    }
    for (int i = 0; i < N; ++i)
    {
        vector<double> v(N);
        for (int j = 0; j < N; ++j) v[j] = U[j].E;
        vector<double> v_left(N), v_right(N);
        MinModLimiter(v, v_left, v_right);
        U_left[i].E = v_left[i];
        U_right[i].E = v_right[i];
    }

    F[0] = flux_roe(U_left[0], U_right[0]);
    F[N] = flux_roe(U_left[N - 1], U_right[N - 1]);
    for (int i = 1; i < N; ++i)
    {
        F[i] = flux_roe(U_right[i - 1], U_left[i]);
    }

    for (int i = 0; i < N; ++i)
    {
        dF[i].rho = (F[i + 1].rho - F[i].rho) / dx;
        dF[i].rho_u = (F[i + 1].rho_u - F[i].rho_u) / dx;
        dF[i].E = (F[i + 1].E - F[i].E) / dx;
    }

    return dF;
}

vector<Conserved> RungeKutta3(const vector<Conserved>& U, double dt, double dx)
{
    vector<Conserved> U1(N), U2(N), U_new(N);
    vector<Conserved> U0 = U;

    vector<Conserved> dF1 = compute_flux_difference(U, dx);
    for (int i = 0; i < N; ++i)
    {
        U1[i].rho = U0[i].rho - dt * dF1[i].rho;
        U1[i].rho_u = U0[i].rho_u - dt * dF1[i].rho_u;
        U1[i].E = U0[i].E - dt * dF1[i].E;
    }

    vector<Conserved> dF2 = compute_flux_difference(U1, dx);
    for (int i = 0; i < N; ++i) {
        U2[i].rho = 0.75 * U0[i].rho + 0.25 * (U1[i].rho - dt * dF2[i].rho);
        U2[i].rho_u = 0.75 * U0[i].rho_u + 0.25 * (U1[i].rho_u - dt * dF2[i].rho_u);
        U2[i].E = 0.75 * U0[i].E + 0.25 * (U1[i].E - dt * dF2[i].E);
    }

    vector<Conserved> dF3 = compute_flux_difference(U2, dx);
    for (int i = 0; i < N; ++i)
    {
        U_new[i].rho = (1.0 / 3.0) * U0[i].rho + (2.0 / 3.0) * (U2[i].rho - dt * dF3[i].rho);
        U_new[i].rho_u = (1.0 / 3.0) * U0[i].rho_u + (2.0 / 3.0) * (U2[i].rho_u - dt * dF3[i].rho_u);
        U_new[i].E = (1.0 / 3.0) * U0[i].E + (2.0 / 3.0) * (U2[i].E - dt * dF3[i].E);
    }

    return U_new;
}

int main()
{
    double dx = (xmax - xmin) / (N - 1);
    vector<Conserved> U = InitCondition(dx);

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

        U = RungeKutta3(U, dt, dx);

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

        fprintf(file, "%lf %lf %lf %lf\n", xmin + i * dx, rho[i], u[i], p[i]);
    }
    fclose(file);

    return 0;
}