// =============================================================
// 程序名: LidDrivenCavity.cpp
// 功能描述: 求解典型的流体力学数值模拟问题： 顶盖驱动的方腔流动
// =============================================================
#include <iostream>
#include <cmath>
using namespace std;

const double EPS = 1e-6;
const int MAX_ITER = 10000;

void InitVelocity(double** u, double** v, int xsize, int ysize)
{
    for (int i = 0; i < xsize; ++i)
    {
        u[i][0] = 0.0;
        u[i][ysize - 1] = sin(M_PI * double(i) / xsize) * sin(M_PI * double(i) / xsize);
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
    for (int i = 1; i < xsize - 1; ++i)
    {
        for (int j = 1; j < ysize - 1; ++j)
        {
            double dv = v[i + 1][j] - v[i - 1][j];
            double du = u[i][j + 1] - u[i][j - 1];
            omega[i][j] = 0.5 * (dv - du) / h;
        }
    }
}
void ApplyBoundCondition(double** psi, double** omega, double** u, double** v, int xsize, int ysize, double h)
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
void UpdateVorticity(double** omega, double** psi, double h, double nu, double dt, int xsize, int ysize)
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
            temp[i][j] = dt * (nu * laplacian + convection) + omega[i][j];
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
bool SOR(double** f, double** omega, int xsize, int ysize, double h, double relax_factor)
{
    double max_diff = 0.0;
    int iter = 0;
    do
    {
        max_diff = 0.0;
        for (int i = 1; i < xsize - 1; ++i)
        {
            for (int j = 1; j < ysize - 1; ++j)
            {
                double old_value = f[i][j];
                f[i][j] = (1 - relax_factor) * old_value + relax_factor * 0.25 * (f[i + 1][j] + f[i][j + 1] + f[i - 1][j] + f[i][j - 1] + h * h * omega[i][j]);
                double diff = abs(f[i][j] - old_value);
                if (diff > max_diff)
                {
                    max_diff = diff;
                }
            }
        }
        ++iter;
    } while (iter < MAX_ITER && max_diff > EPS);

    if (iter >= MAX_ITER || iter == 1)
    {
        return false;
    }
    return true;
}

int main()
{
    cout << "Please enter the grid spacing (h): ";
    double h;
    cin >> h;
    double nu = 0.001;
    double relax_factor = 1.5;
    double dt = 0.01;
    bool is_converge = false;
    int xsize = static_cast<int>(1.0 / h);
    int ysize = static_cast<int>(1.0 / h);
    double** u = new double*[xsize];
    double** v = new double*[xsize];
    double** omega = new double*[xsize];
    double** psi = new double*[xsize];
    for (int i = 0; i < xsize; ++i)
    {
        u[i] = new double[ysize];
        v[i] = new double[ysize];
        omega[i] = new double[ysize];
        psi[i] = new double[ysize];
        for (int j = 0; j < ysize; ++j)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            omega[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
    }

    InitVelocity(u, v, xsize, ysize);
    InitVorticity(u, v, omega, xsize, ysize, h);
    is_converge = SOR(psi, omega, xsize, ysize, h, relax_factor);

    while (is_converge)
    {
        ApplyBoundCondition(psi, omega, u, v, xsize, ysize, h);
        UpdateVorticity(omega, psi, h, nu, dt, xsize, ysize);
        is_converge = SOR(psi, omega, xsize, ysize, h, relax_factor);
    }

    FILE* fp1 = fopen("StreamFunction.txt", "w");
    if (fp1 == NULL)
    {
        cout << "Error opening file: StreamFunction.txt!" << endl;
        return -1;
    }
    for (int i = 0; i < xsize; ++i)
    {
        for (int j = 0; j < ysize; ++j)
        {
            fprintf(fp1, "%d %d %lf\n", i, j, psi[i][j]);
        }
    }
    fclose(fp1);
    FILE* fp2 = fopen("Vorticity.txt", "w");
    if (fp2 == NULL)
    {
        cout << "Error opening file: Vorticity.txt!" << endl;
        return -1;
    }
    for (int i = 0; i < xsize; ++i)
    {
        for (int j = 0; j < ysize; ++j)
        {
            fprintf(fp2, "%d %d %lf\n", i, j, omega[i][j]);
        }
    }
    fclose(fp2);

    for (int i = 0; i < xsize; ++i)
    {
        delete[] u[i];
        delete[] v[i];
        delete[] omega[i];
        delete[] psi[i];
    }
    delete[] u;
    delete[] v;
    delete[] omega;
    delete[] psi;
    cout << "Simulation completed. Results saved to LidDrivenCavity.txt" << endl;

    return 0;
}