// =============================================================
// 程序名: HW3
// 功能描述: 用Lax格式、Lax-Wendroff格式和一阶迎风格式求解波动方程：
//          ∂u/∂t + ∂u/∂x = 0
//          u(x, 0) = sin(2πx) (0 < x < 3)
// 输入参数:
// N: 网格节点数
// c: CFL数
// T: 计算时间
// 输出参数：
// 
// =============================================================
#include <iostream>
#include <cmath>

using namespace std;

void Lax(double **u, int t, double c, int N)
{
    for (int i = 1; i < N - 1; i++)
    {
        u[i][t + 1] = 0.5 * (1 - c) * u[i + 1][t] + 0.5 * (1 + c) * u[i - 1][t];
    }
    u[0][t + 1] = 0.5 * (1 - c) * u[1][t] + 0.5 * (1 + c) * u[N - 1][t];
    u[N - 1][t + 1] = 0.5 * (1 - c) * u[0][t] + 0.5 * (1 + c) * u[N - 2][t];
    return;
}
void LaxWendroff(double **u, int t, double c, int N)
{
    for (int i = 1; i < N - 1; i++)
    {
        u[i][t + 1] = u[i][t] - 0.5 * c * (u[i + 1][t] - u[i - 1][t]) + 0.5 * c * c * (u[i + 1][t] - 2 * u[i][t] + u[i - 1][t]);
    }
    u[0][t + 1] = u[0][t] - 0.5 * c * (u[1][t] - u[N - 1][t]) + 0.5 * c * c * (u[1][t] - 2 * u[0][t] + u[N - 1][t]);
    u[N - 1][t + 1] = u[N - 1][t] - 0.5 * c * (u[0][t] - u[N - 2][t]) + 0.5 * c * c * (u[0][t] - 2 * u[N - 1][t] + u[N - 2][t]);
    return;
}
void Upwind(double **u, int t, double c, int N)
{
    for (int i = 1; i <= N - 1; i++)
    {
        u[i][t + 1] = u[i][t] - c * (u[i][t] - u[i - 1][t]);
    }
    u[0][t + 1] = u[0][t] - c * (u[0][t] - u[N - 1][t]);
    return;
}

int main() {
    int N = 0;
    double c = 0, T = 0;
    cin >> N >> c >> T;

    double dx = 3.0 / N;
    double dt = c * dx;
    int M = int(T / dt);

    // 动态分配二维数组
    double **u_lax = new double*[N];
    double **u_lax_wendroff = new double*[N];
    double **u_upwind = new double*[N];
    double **u_exact = new double*[N];
    for (int i = 0; i < N; i++) {
        u_lax[i] = new double[M+1];
        u_lax_wendroff[i] = new double[M+1];
        u_upwind[i] = new double[M+1];
        u_exact[i] = new double[M+1];
    }

    // 初始化
    for (int i = 0; i < N; i++) {
        u_lax[i][0] = sin(2 * M_PI * i * dx);
        u_lax_wendroff[i][0] = u_lax[i][0];
        u_upwind[i][0] = u_lax[i][0];
        u_exact[i][0] = u_lax[i][0];
    }

    // 计算数值解
    for (int t = 0; t < M; t++) {
        Lax(u_lax, t, c, N);
        LaxWendroff(u_lax_wendroff, t, c, N);
        Upwind(u_upwind, t, c, N);
    }

    // 输出结果
    cout << "Lax Scheme Result:" << endl;
    for (int i = 0; i < N; i++) {
        cout << u_lax[i][M] << " ";
    }
    cout << endl;
    cout << "Lax-Wendroff Scheme Result:" << endl;
    for (int i = 0; i < N; i++) {
        cout << u_lax_wendroff[i][M] << " ";
    }
    cout << endl;
    cout << "Upwind Scheme Result:" << endl;
    for (int i = 0; i < N; i++) {
        cout << u_upwind[i][M] << " ";
    }
    cout << endl;

    // 释放内存
    for (int i = 0; i < N; i++) {
        delete[] u_lax[i];
        delete[] u_lax_wendroff[i];
        delete[] u_upwind[i];
        delete[] u_exact[i];
    }
    delete[] u_lax;
    delete[] u_lax_wendroff;
    delete[] u_upwind;
    delete[] u_exact;

    return 0;
}