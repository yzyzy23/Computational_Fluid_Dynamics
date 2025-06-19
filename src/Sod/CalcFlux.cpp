#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "CalcDiff.cpp"

using namespace Eigen;
using namespace std;

/*void CalcFluxFDS(const vector<Conserved>& U, double dx, vector<Conserved>& F, vector<Conserved>& F_roe, int diff_calc_type)
{
    int n = U.size();
    vector<Conserved> Conserve_left(n), Conserve_right(n), diff(n);
    for (int i = 0; i < n; ++i)
    {
        Conserve_left[i] = {0.0, 0.0, 0.0};
        Conserve_right[i] = {0.0, 0.0, 0.0};
    }

    CalcDiff(U, U, Conserve_left, Conserve_right, diff, dx, diff_calc_type);

    // Calculate Roe averages
    for (int i = 3; i < n - 2; ++i)
    {
        Conserved U_left = Conserve_left[i];
        Conserved U_right = Conserve_right[i];
        State left = ConservedToState(U_left);
        State right = ConservedToState(U_right);
        Conserved F_left = {left.rho * left.u, left.rho * left.u * left.u + left.p, left.u * (U_left.E + left.p)};
        Conserved F_right = {right.rho * right.u, right.rho * right.u * right.u + right.p, right.u * (U_right.E + right.p)};
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
        double local_eps = 0.1;
        for (int j = 0; j < 3; ++j)
        {
            if (abs(G(j, j)) > local_eps)
            {
                G_abs(j, j) = abs(G(j, j));
            } 
            else
            {
                //G_abs(j, j) = ((G(j, j) * G(j, j)) + (local_eps * local_eps)) / (2 * local_eps);
                G_abs(j, j) = local_eps;
            }
        }

        MatrixXd V_inv = V.inverse();
        MatrixXd A_avg_abs = V * G_abs * V_inv;

        F[i].rho = 0.5 * (F_left.rho + F_right.rho) - 0.5 * (A_avg_abs(0, 0) * (U_right.rho - U_left.rho) + A_avg_abs(0, 1) * (U_right.rho_u - U_left.rho_u) + A_avg_abs(0, 2) * (U_right.E - U_left.E));
        F[i].rho_u = 0.5 * (F_left.rho_u + F_right.rho_u) - 0.5 * (A_avg_abs(1, 0) * (U_right.rho - U_left.rho) + A_avg_abs(1, 1) * (U_right.rho_u - U_left.rho_u) + A_avg_abs(1, 2) * (U_right.E - U_left.E));
        F[i].E = 0.5 * (F_left.E + F_right.E) - 0.5 * (A_avg_abs(2, 0) * (U_right.rho - U_left.rho) + A_avg_abs(2, 1) * (U_right.rho_u - U_left.rho_u) + A_avg_abs(2, 2) * (U_right.E - U_left.E));

    }
    for (int i = 4; i < n - 2; ++i)
    {
        F_roe[i].rho = (F[i].rho - F[i - 1].rho) / dx;
        F_roe[i].rho_u = (F[i].rho_u - F[i - 1].rho_u) / dx;
        F_roe[i].E = (F[i].E - F[i - 1].E) / dx;
    }
}*/
/*void CalcFluxFDS(const vector<Conserved>& U, double dx, vector<Conserved>& dF) {
    int n = U.size();

    // 一阶插值，计算 U_l 和 U_r
    vector<Conserved> U_l(U.begin() + 1, U.end() - 2);
    vector<Conserved> U_r(U.begin() + 2, U.end() - 1);

    // 计算左、右状态的物理量
    vector<double> rho_l(U_l.size()), u_l(U_l.size()), e_l(U_l.size()), p_l(U_l.size()), H_l(U_l.size());
    vector<double> rho_r(U_r.size()), u_r(U_r.size()), e_r(U_r.size()), p_r(U_r.size()), H_r(U_r.size());

    for (int i = 0; i < U_l.size(); ++i) {
        rho_l[i] = U_l[i].rho;
        u_l[i] = U_l[i].rho_u / U_l[i].rho;
        e_l[i] = U_l[i].E / U_l[i].rho;
        p_l[i] = (gamma_val - 1) * (e_l[i] - 0.5 * u_l[i] * u_l[i]) * rho_l[i];
        H_l[i] = e_l[i] + p_l[i] / rho_l[i];

        rho_r[i] = U_r[i].rho;
        u_r[i] = U_r[i].rho_u / U_r[i].rho;
        e_r[i] = U_r[i].E / U_r[i].rho;
        p_r[i] = (gamma_val - 1) * (e_r[i] - 0.5 * u_r[i] * u_r[i]) * rho_r[i];
        H_r[i] = e_r[i] + p_r[i] / rho_r[i];
    }

    // 计算 Roe 平均值
    vector<double> rho_roe(U_l.size()), u_roe(U_l.size()), H_roe(U_l.size()), p_roe(U_l.size()), c_roe(U_l.size());
    for (int i = 0; i < U_l.size(); ++i) {
        rho_roe[i] = pow(0.5 * (sqrt(rho_l[i]) + sqrt(rho_r[i])), 2);
        u_roe[i] = (sqrt(rho_l[i]) * u_l[i] + sqrt(rho_r[i]) * u_r[i]) / (sqrt(rho_l[i]) + sqrt(rho_r[i]));
        H_roe[i] = (sqrt(rho_l[i]) * H_l[i] + sqrt(rho_r[i]) * H_r[i]) / (sqrt(rho_l[i]) + sqrt(rho_r[i]));
        p_roe[i] = (gamma_val - 1) / gamma_val * (rho_roe[i] * H_roe[i] - 0.5 * rho_roe[i] * u_roe[i] * u_roe[i]);
        c_roe[i] = sqrt((gamma_val - 1) * (H_roe[i] - 0.5 * u_roe[i] * u_roe[i]));
    }

    // 计算特征分解矩阵
    vector<vector<vector<double>>> lamda(U_l.size(), vector<vector<double>>(3, vector<double>(3, 0.0)));
    vector<vector<vector<double>>> S(U_l.size(), vector<vector<double>>(3, vector<double>(3, 0.0)));
    vector<vector<vector<double>>> S_inv(U_l.size(), vector<vector<double>>(3, vector<double>(3, 0.0)));
    vector<vector<vector<double>>> A_abs(U_l.size(), vector<vector<double>>(3, vector<double>(3, 0.0)));

    for (int i = 0; i < U_l.size(); ++i) {
        lamda[i][0][0] = u_roe[i];
        lamda[i][1][1] = u_roe[i] - c_roe[i];
        lamda[i][2][2] = u_roe[i] + c_roe[i];

        S[i][0][0] = 1;
        S[i][0][1] = 1;
        S[i][0][2] = 1;
        S[i][1][0] = u_roe[i];
        S[i][1][1] = u_roe[i] - c_roe[i];
        S[i][1][2] = u_roe[i] + c_roe[i];
        S[i][2][0] = 0.5 * u_roe[i] * u_roe[i];
        S[i][2][1] = 0.5 * u_roe[i] * u_roe[i] + gamma_val / (gamma_val - 1) * p_roe[i] / rho_roe[i] - u_roe[i] * c_roe[i];
        S[i][2][2] = 0.5 * u_roe[i] * u_roe[i] + gamma_val / (gamma_val - 1) * p_roe[i] / rho_roe[i] + u_roe[i] * c_roe[i];

        S_inv[i][0][0] = 1 - (gamma_val - 1) * u_roe[i] * u_roe[i] / (2 * c_roe[i] * c_roe[i]);
        S_inv[i][0][1] = (gamma_val - 1) * u_roe[i] / (c_roe[i] * c_roe[i]);
        S_inv[i][0][2] = -(gamma_val - 1) / (c_roe[i] * c_roe[i]);
        S_inv[i][1][0] = (gamma_val - 1) * u_roe[i] * u_roe[i] / (4 * c_roe[i] * c_roe[i]) + 0.5 * u_roe[i] / c_roe[i];
        S_inv[i][1][1] = -(gamma_val - 1) * u_roe[i] / (2 * c_roe[i] * c_roe[i]) - 1 / (2 * c_roe[i]);
        S_inv[i][1][2] = (gamma_val - 1) / (2 * c_roe[i] * c_roe[i]);
        S_inv[i][2][0] = (gamma_val - 1) * u_roe[i] * u_roe[i] / (4 * c_roe[i] * c_roe[i]) - 0.5 * u_roe[i] / c_roe[i];
        S_inv[i][2][1] = -(gamma_val - 1) * u_roe[i] / (2 * c_roe[i] * c_roe[i]) + 1 / (2 * c_roe[i]);
        S_inv[i][2][2] = (gamma_val - 1) / (2 * c_roe[i] * c_roe[i]);
    }

    for (int i = 0; i < U_l.size(); ++i) {
        vector<vector<double>> aa(3, vector<double>(3, 0.0));
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                aa[j][k] = S[i][j][k] * abs(lamda[i][j][k]);
            }
        }
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                A_abs[i][j][k] = 0.0;
                for (int l = 0; l < 3; ++l) {
                    A_abs[i][j][k] += aa[j][l] * S_inv[i][l][k];
                }
            }
        }
    }

    vector<Conserved> F(U_l.size());
    for (int i = 0; i < U_l.size(); ++i) {
        Conserved F_l = {
            rho_l[i] * u_l[i],
            rho_l[i] * u_l[i] * u_l[i] + p_l[i],
            u_l[i] * (rho_l[i] * e_l[i] + p_l[i])
        };

        Conserved F_r = {
            rho_r[i] * u_r[i],
            rho_r[i] * u_r[i] * u_r[i] + p_r[i],
            u_r[i] * (rho_r[i] * e_r[i] + p_r[i])
        };

        F[i].rho = 0.5 * (F_l.rho + F_r.rho);
        F[i].rho_u = 0.5 * (F_l.rho_u + F_r.rho_u);
        F[i].E = 0.5 * (F_l.E + F_r.E);

        for (int j = 0; j < 3; ++j) {
            F[i].rho -= 0.5 * A_abs[i][0][j] * (U_r[i].rho - U_l[i].rho);
            F[i].rho_u -= 0.5 * A_abs[i][1][j] * (U_r[i].rho_u - U_l[i].rho_u);
            F[i].E -= 0.5 * A_abs[i][2][j] * (U_r[i].E - U_l[i].E);
        }
    }

    for (int i = 0; i < F.size() - 1; ++i) {
        dF[i].rho = (F[i + 1].rho - F[i].rho) / dx;
        dF[i].rho_u = (F[i + 1].rho_u - F[i].rho_u) / dx;
        dF[i].E = (F[i + 1].E - F[i].E) / dx;
    }
}*/
void CalcFluxFDS(const vector<Conserved>& U, vector<Conserved>& fluxes)
{
    int ext_size = U.size();
    int num_faces = ext_size - 1;
    fluxes.resize(num_faces);

    for (int i = 0; i < num_faces; ++i)
    {
        Conserved U_l = U[i];
        Conserved U_r = U[i+1];

        State s_l = ConservedToState(U_l);
        State s_r = ConservedToState(U_r);

        // Calculate Roe averages
        double sqrt_rho_l = sqrt(s_l.rho);
        double sqrt_rho_r = sqrt(s_r.rho);
        double sum_sqrt = sqrt_rho_l + sqrt_rho_r;

        double rho_roe = 0.25 * pow(sum_sqrt, 2);
        double u_roe = (sqrt_rho_l * s_l.u + sqrt_rho_r * s_r.u) / sum_sqrt;

        double H_l = (s_l.p * gamma_val / (gamma_val - 1.0) + 
                0.5 * s_l.rho * s_l.u * s_l.u) / s_l.rho;
        double H_r = (s_r.p * gamma_val / (gamma_val - 1.0) + 
                0.5 * s_r.rho * s_r.u * s_r.u) / s_r.rho;

        double H_roe = (sqrt_rho_l * H_l + sqrt_rho_r * H_r) / sum_sqrt;

        double a_roe = sqrt((gamma_val - 1.0) * 
                    (H_roe - 0.5 * u_roe * u_roe));

        double lambda1 = abs(u_roe);
        double lambda2 = abs(u_roe - a_roe);
        double lambda3 = abs(u_roe + a_roe);

        Conserved F_l, F_r;
        F_l.rho = s_l.rho * s_l.u;
        F_l.rho_u = F_l.rho * s_l.u + s_l.p;
        F_l.E = s_l.u * (s_l.p * gamma_val / (gamma_val - 1.0) + 
                    0.5 * s_l.rho * s_l.u * s_l.u);

        F_r.rho = s_r.rho * s_r.u;
        F_r.rho_u = F_r.rho * s_r.u + s_r.p;
        F_r.E = s_r.u * (s_r.p * gamma_val / (gamma_val - 1.0) + 
                    0.5 * s_r.rho * s_r.u * s_r.u);

        // Calculate Roe fluxes
        double diff_rho = U_r.rho - U_l.rho;
        double diff_rho_u = U_r.rho_u - U_l.rho_u;
        double diff_E = U_r.E - U_l.E;

        double delta1 = (gamma_val - 1.0) / (a_roe * a_roe) * 
                ((H_roe - u_roe * u_roe) * diff_rho + 
                    u_roe * diff_rho_u - diff_E);
        double delta2 = (diff_rho + (u_roe + a_roe) * delta1 - 
                    (diff_rho_u + u_roe * diff_rho)) / (2.0 * a_roe);
        double delta3 = (diff_rho_u + u_roe * diff_rho) - delta2;

        Conserved absA_deltaU;
        absA_deltaU.rho = lambda1 * ((delta1 * a_roe * a_roe / (gamma_val - 1.0))) + 
                    lambda2 * delta2 + 
                    lambda3 * delta3;

        absA_deltaU.rho_u = lambda1 * (u_roe * (delta1 * a_roe * a_roe / (gamma_val - 1.0)) + 
                        a_roe * (delta3 - delta2)) + 
            lambda2 * ((u_roe - a_roe) * delta2) + 
            lambda3 * ((u_roe + a_roe) * delta3);

        absA_deltaU.E = lambda1 * (0.5 * u_roe * u_roe * (delta1 * a_roe * a_roe / (gamma_val - 1.0)) + 
                    u_roe * a_roe * (delta3 - delta2) / 2.0) + 
        lambda2 * (H_roe - u_roe * a_roe) * delta2 + 
        lambda3 * (H_roe + u_roe * a_roe) * delta3;

        fluxes[i].rho = 0.5 * (F_l.rho + F_r.rho) - 0.5 * absA_deltaU.rho;
        fluxes[i].rho_u = 0.5 * (F_l.rho_u + F_r.rho_u) - 0.5 * absA_deltaU.rho_u;
        fluxes[i].E = 0.5 * (F_l.E + F_r.E) - 0.5 * absA_deltaU.E;
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