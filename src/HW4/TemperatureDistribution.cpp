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
    for (int i = 0; i <= xsize; ++i)
    {
        f[i][0] = 293.15;
        f[i][ysize] = 373.15;
    }
    for (int j = 0; j <= ysize; ++j)
    {
        f[0][j] = 293.15;
        f[xsize][j] = 293.15;
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
    double h, omega;
    cout << "Enter h, the grid spacing: ";
    cin >> h;
    cout << "Enter omega, the relaxation parameter: ";
    cin >> omega;

    // initialize T with grid size
    int N = 15, M = 12;
    int xsize = int(N / h);
    int ysize = int(M / h);
    double** T = new double*[xsize + 1];
    for (int i = 0; i <= xsize; ++i)
    {
        T[i] = new double[ysize + 1];
    }

    // execute SOR method
    GridInit(T, xsize, ysize);
    int max_iter = 100000;
    int iter = SOR(T, xsize, ysize, EPS, max_iter, omega);

    // output results to file
    FILE* fp = fopen("TemperatureDistribution.txt", "w");
    if (fp == NULL)
    {
        cout << "Error opening file!" << endl;
        return 1;
    }
    fprintf(fp, "%d %d\n", xsize, ysize);
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

    // release memory for T
    for (int i = 0; i <= xsize; ++i)
    {
        delete[] T[i];
    }
    delete[] T;

    return 0;
}