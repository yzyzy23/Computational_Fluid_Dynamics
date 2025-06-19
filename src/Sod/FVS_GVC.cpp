#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// Initialize constants
const double dx = 0.005;
const double end_time = 0.14; // End time for the simulation
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

double CalcSoundSpeed(const State& state)
{
    return sqrt(gamma_val * state.p / state.rho);
}

void InitCondition(vector<State>& states, vector<Conserved>& U) {
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
    for (int i = 0; i < n; ++i)
    {
        U[i] = StateToConserved(states[i]);
    }
}
double GVCLimiter(double r)
{
    // GVC limiter
    return abs(r) > 1.0 ? 1.0 : r;
}
void CalcDiffGVC(const vector<Conserved>& U_left, const vector<Conserved>& U_right, vector<Conserved>& new_u_left, vector<Conserved>& new_u_right, vector<Conserved>& diff, double dx)
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
        State state = ConservedToState(U[i]);
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
    CalcDiffGVC(f_left, f_right, temp_left, temp_right, diff, dx);
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
    CalcDiffGVC(f_left, f_right, temp_left, temp_right, diff, dx);
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
    CalcDiffGVC(f_left, f_right, temp_left, temp_right, diff, dx);
    for (int i = 0; i < n; ++i)
    {
        U[i].rho = (1.0/3.0) * U0[i].rho + (2.0/3.0) * (U2[i].rho - dt * diff[i].rho);
        U[i].rho_u = (1.0/3.0) * U0[i].rho_u + (2.0/3.0) * (U2[i].rho_u - dt * diff[i].rho_u);
        U[i].E = (1.0/3.0) * U0[i].E + (2.0/3.0) * (U2[i].E - dt * diff[i].E);
    }
}
int main() 
{
    int n = static_cast<int>(2.0 / dx);
    double dt = 0.0;

    vector<State> states(n);
    vector<Conserved> U(n);
    InitCondition(states, U);

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

        time += dt;
        RungeKutta3(U, dt, dx);
        
        for (int i = 0; i < n; ++i)
        {
            states[i] = ConservedToState(U[i]);
        }
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
        fprintf(file, "%lf %lf %lf %lf\n", x_min + i * dx, states[i].rho, states[i].u, states[i].p);
    }
    fclose(file);

    return 0;
}