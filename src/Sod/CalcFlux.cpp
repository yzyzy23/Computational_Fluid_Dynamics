#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "CalcDiff.cpp"

using namespace Eigen;
using namespace std;

void CalcFluxFDS(const vector<Conserved>& U, double dx, vector<Conserved>& F, vector<Conserved>& F_roe)
{
    int n = U.size();
    vector<Conserved> Conserve_left(n), Conserve_right(n), diff(n);
    for (int i = 0; i < n; ++i)
    {
        Conserve_left[i] = {0.0, 0.0, 0.0};
        Conserve_right[i] = {0.0, 0.0, 0.0};
    }

    CalcDiff(U, U, Conserve_left, Conserve_right, diff, dx, 1);

    // Calculate Roe averages
    for (int i = 3; i < n - 2; ++i)
    {
        Conserved U_left = Conserve_left[i];
        Conserved U_right = Conserve_right[i];
        State left = ConservedToState(U_left);
        State right = ConservedToState(U_right);
        Conserved F_left = {left.rho, left.rho * left.u * left.u + left.p, left.u * (U_left.E + left.p)};
        Conserved F_right = {right.rho, right.rho * right.u * right.u + right.p, right.u * (U_right.E + right.p)};
        double sqrt_rho_left = sqrt(left.rho);
        double sqrt_rho_right = sqrt(right.rho);
        double H_left = (U_left.E + left.p) / left.rho;
        double H_right = (U_right.E + right.p) / right.rho;
        double u_avg = (sqrt_rho_left * left.u + sqrt_rho_right * right.u) / (sqrt_rho_left + sqrt_rho_right);
        double H_avg = (sqrt_rho_left * H_left + sqrt_rho_right * H_right) / (sqrt_rho_left + sqrt_rho_right);
        double c_avg = sqrt((gamma_val - 1.0) * (H_avg - 0.5 * u_avg * u_avg));

        MatrixXd A = MatrixXd::Zero(3, 3);
        A(0, 0) = 0;
        A(0, 1) = 1.0;
        A(0, 2) = 0;
        A(1, 0) = (-1) * ((3 - gamma_val) / 2) * u_avg * u_avg;
        A(1, 1) = (3 - gamma_val) * u_avg;
        A(1, 2) = gamma_val - 1.0;
        A(2, 0) = (((gamma_val - 2) / 2) * u_avg * u_avg * u_avg) - ((u_avg * c_avg * c_avg) / (gamma_val - 1));
        A(2, 1) = ((c_avg * c_avg) / (gamma_val - 1)) + (((3 - gamma_val) / 2) * u_avg * u_avg);
        A(2, 2) = gamma_val * u_avg;

        SelfAdjointEigenSolver<MatrixXd> solver(A);
        MatrixXd V = solver.eigenvectors().real();
        MatrixXd G = solver.eigenvalues().real().asDiagonal();

        MatrixXd G_abs = MatrixXd::Zero(3, 3);
        // Entropy modification
        for (int i = 0; i < 3; ++i)
        {
            if (abs(G(i, i)) > eps)
            {
                G_abs(i, i) = abs(G(i, i));
            } 
            else
            {
                G_abs(i, i) = ((G(i, i) * G(i, i)) + (eps * eps)) / (2 * eps);
            }
        }


        MatrixXd V_inv = V.inverse();
        MatrixXd A_ave_abs = V_inv * G_abs * V;

        F[i].rho = 0.5 * (F_left.rho + F_right.rho) - 0.5 * (A_ave_abs(0, 0) * (U_right.rho - U_left.rho) + A_ave_abs(0, 1) * (U_right.rho_u - U_left.rho_u) + A_ave_abs(0, 2) * (U_right.E - U_left.E));
        F[i].rho_u = 0.5 * (F_left.rho_u + F_right.rho_u) - 0.5 * (A_ave_abs(1, 0) * (U_right.rho - U_left.rho) + A_ave_abs(1, 1) * (U_right.rho_u - U_left.rho_u) + A_ave_abs(1, 2) * (U_right.E - U_left.E));
        F[i].E = 0.5 * (F_left.E + F_right.E) - 0.5 * (A_ave_abs(2, 0) * (U_right.rho - U_left.rho) + A_ave_abs(2, 1) * (U_right.rho_u - U_left.rho_u) + A_ave_abs(2, 2) * (U_right.E - U_left.E));

        /*// Entropy modification
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
        F_roe[i].rho = 0.5 * (F_left.rho + F_right.rho) - 0.5 * F.rho;
        F_roe[i].rho_u = 0.5 * (F_left.rho_u + F_right.rho_u) - 0.5 * F.rho_u;
        F_roe[i].E = 0.5 * (F_left.E + F_right.E) - 0.5 * F.E;*/
    }
    for (int i = 4; i < n - 2; ++i)
    {
        F_roe[i].rho = (F[i].rho - F[i - 1].rho) / dx;
        F_roe[i].rho_u = (F[i].rho_u - F[i - 1].rho_u) / dx;
        F_roe[i].E = (F[i].E - F[i - 1].E) / dx;
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

        lambda[0] = u;
        lambda[1] = u - c;
        lambda[2] = u + c;

        for (int j = 0; j < 3; ++j)
        {
            lambda_left[j] = (lambda[j] + sqrt(lambda[j] * lambda[j] + eps * eps)) / 2.0;
            lambda_right[j] = (lambda[j] - sqrt(lambda[j] * lambda[j] + eps * eps)) / 2.0;
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
/*Conserved CalcFluxFDS(const State& left, const State& right)
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
    fluxes[0] = CalcFluxFDS(left, states[0]);
    State right = states[n-1];
    right.u = -right.u;
    fluxes[n] = CalcFluxFDS(states[n-1], right);

    // Calculate fluxes for interior states
    for (int i = 1; i < n; ++i)
    {
        fluxes[i] = CalcFluxFDS(states[i-1], states[i]);
    }
}*/
