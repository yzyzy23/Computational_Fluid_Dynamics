// =============================================================
// 程序名: LidDrivenCavity.cpp
// 功能描述: 求解典型的流体力学数值模拟问题： 顶盖驱动的方腔流动
// =============================================================
#include <iostream>
#include <cmath>
#include "SOR.cpp"
using namespace std;

const double EPS = 1e-5;

void InitVelocity(double** u, double** v, int xsize, int ysize)
{
    for (int i = 0; i < xsize; ++i)
    {
        u[i][0] = 0.0;
        u[i][ysize - 1] = sin(M_PI * i / xsize) * sin(M_PI * i / xsize);
        v[i][0] = 0.0;
        v[i][ysize - 1] = 0.0;
    }
    for (int i = 0; i < ysize; ++i)
    {
        u[0][i] = 0.0;
        u[xsize - 1][i] = 0.0;
        v[0][i] = 0.0;
        v[xsize - 1][i] = 0.0;
    }
}
void InitVorticity(double** u, double** v, double** omega, int xsize, int ysize, double h)
{
    for (int i = 0; i < xsize; ++i)
    {
        for (int j = 0; j < ysize; ++j)
        {
            double dv = v[i + 1][j] - v[i - 1][j];
            double du = u[i][j + 1] - u[i][j - 1];
            omega[i][j] = 0.5 * (dv - du) / h;
        }
    }
}
void ApplyBoundCondition(double** psi, double** u, double** omega, double** v, int xsize, int ysize, double h)
{
    for (int i = 0; i < xsize; ++i)
    {
        omega[i][0] = 2 * (psi[i][0] - psi[i][1]) / (h * h) + 2 * u[i][0] / h;
        omega[i][ysize - 1] = 2 * (psi[i][ysize - 1] - psi[i][ysize - 2]) / (h * h) - 2 * u[i][ysize - 1] / h;
    }
    for (int j = 0; j < ysize; ++j)
    {
        omega[0][j] = 2 * (psi[0][j] - psi[1][j]) / (h * h) - 2 * v[0][j] / h;
        omega[xsize - 1][j] = 2 * (psi[xsize - 1][j] - psi[xsize - 2][j]) / (h * h) + 2 * v[xsize - 1][j] / h;
    }
}
void UpdateVorticity(double** omega, double** psi, double h, double nu, int xsize, int ysize)
{
    double** temp = new double*[xsize];
    for (int i = 0; i < xsize; ++i)
    {
        temp[i] = new double[ysize];
    }
    for (int i = 1; i < xsize - 1; ++i)
    {
        for (int j = 1; j < ysize - 1; ++j)
        {
            double laplacian = (omega[i + 1][j] + omega[i - 1][j] + omega[i][j + 1] + omega[i][j - 1] - 4 * omega[i][j]) / (h * h);
            double convection = (psi[i + 1][j] - psi[i - 1][j]) * (omega[i][j + 1] - omega[i][j - 1]) / (4 * h * h) - (psi[i][j + 1] - psi[i][j - 1]) * (omega[i + 1][j] - omega[i - 1][j]) / (4 * h * h);
            temp[i][j] = nu * laplacian + convection + omega[i][j];
        }
    }
    for (int i = 1; i < xsize - 1; ++i)
    {
        for (int j = 1; j < ysize - 1; ++j)
        {
            omega[i][j] = temp[i][j];
        }
    }
    for (int i = 0; i < xsize; ++i)
    {
        delete[] temp[i];
    }
    delete[] temp;
}

int main()
{
    return 0;
}