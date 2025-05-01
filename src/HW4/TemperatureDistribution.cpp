// =============================================================
// 程序名: TemperatureDistribution.cpp
// 功能描述: 用SOR方法求解二维稳态热传导方程
//          ∂²u/∂x² + ∂²u/∂y² = 0
// 输入参数:
// h: 网格间距
// omega: 松弛因子
// 输出参数：
// 温度分布被存储在文件中，方便利用python进行可视化
// 在终端输出迭代次数
// 为了测试松弛因子和步长对收敛速率的影响，在输入时可以选择开启测试循环
// 这将对松弛因子在0.8-2.0，步长在0.1，0.2，0.25，0.5，1.0进行测试
// =============================================================
#include <iostream>
using namespace std;

const double EPS = 1e-5;

void GridInit(double** f, int xsize, int ysize)
{   
    for (int i = 1; i < xsize; ++i)
    {
        for (int j = 1; j < ysize; ++j)
        {
            f[i][j] = 0.0;
        }
    }
    for (int j = 0; j <= ysize; ++j)
    {
        f[0][j] = 293.15;
        f[xsize][j] = 293.15;
    }
    for (int i = 0; i <= xsize; ++i)
    {
        f[i][0] = 293.15;
        f[i][ysize] = 373.15;
    }
}
int SOR(double** f, int xsize, int ysize, double eps, int max_iter, double omega)
{
    double max_diff = 0.0;
    int iter = 0;
    do
    {
        max_diff = 0.0;
        for (int i = 1; i < xsize; ++i)
        {
            for (int j = 1; j < ysize; ++j)
            {
                double old_value = f[i][j];
                f[i][j] = (1 - omega) * old_value + omega * 0.25 * (f[i + 1][j] + f[i][j + 1] + f[i - 1][j] + f[i][j - 1]);
                double diff = abs(f[i][j] - old_value);
                if (diff > max_diff)
                {
                    max_diff = diff;
                }
            }
        }
        ++iter;
    } while (iter < max_iter && max_diff > eps);
    return iter;
}

int main()
{
    double h = 1.0, omega = 1.0;
    cout << "Start convergence rate test? (1 for yes, 0 for no): ";
    bool test;
    cin >> test;
    if (!test)
    {
        cout << "Enter h, the grid spacing: ";
        cin >> h;
        cout << "Enter omega, the relaxation parameter: ";
        cin >> omega;
    }
    
    // initialize T with grid size
    int N = 15, M = 12;
    int xsize = int(N / h);
    int ysize = int(M / h);
    int max_iter = 100000;
    double** T = new double*[xsize + 1];
    for (int i = 0; i <= xsize; ++i)
    {
        T[i] = new double[ysize + 1];
    }

    // execute SOR method
    if (!test)
    {
        GridInit(T, xsize, ysize);
        int iter = SOR(T, xsize, ysize, EPS, max_iter, omega);

        // output results to file
        FILE* fp = fopen("TemperatureDistribution.txt", "w");
        if (fp == NULL)
        {
            cout << "Error opening file: TemperatureDistribution.txt!" << endl;
            return 1;
        }
        for (int i = 0; i <= xsize; ++i)
        {
            for (int j = 0; j <= ysize; ++j)
            {
                fprintf(fp, "%lf ", T[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
        cout << "SOR method completed in " << iter << " iterations." << endl;
    }
    else
    {
        // reinitialize T for convergence rate test
        for (int i = 0; i <= xsize; ++i)
        {
            delete[] T[i];
        }
        delete[] T;

        FILE* fp = fopen("ConvergenceRate.txt", "w");
        if (fp == NULL)
        {
            cout << "Error opening file: ConvergenceRate.txt!" << endl;
            return 1;
        }
        double min_omega = 0.8, max_omega = 2, step = 0.01;
        double h_arr[4] = { 0.2, 0.25, 0.5, 1.0};
        for (int i = 0; i < 4; ++i)
        {
            h = h_arr[i];
            xsize = int(N / h);
            ysize = int(M / h);
            double** T = new double*[xsize + 1];
            for (int i = 0; i <= xsize; ++i)
            {
                T[i] = new double[ysize + 1];
            }
            fprintf(fp, "%lf\n", h);
            for (double omega = min_omega; omega <= max_omega; omega += step)
            {
                GridInit(T, xsize, ysize);
                int iter = SOR(T, xsize, ysize, EPS, max_iter, omega);
                fprintf(fp, "%lf %d\n", omega, iter);
            }
            for (int i = 0; i <= xsize; ++i)
            {
                delete[] T[i];
            }
            delete[] T;
            cout << "Convergence Rate test complete for h = " << h << endl;
        }
        fclose(fp);
    }

    // release memory for T
    if (!test)
    {
        for (int i = 0; i <= xsize; ++i)
        {
            delete[] T[i];
        }
        delete[] T;
    }
    

    return 0;
}