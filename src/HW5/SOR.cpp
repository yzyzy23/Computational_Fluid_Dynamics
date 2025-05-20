#include <iostream>
using namespace std;

const double EPS = 1e-5;
const int MAX_ITER = 10000;

int SOR(double** f, int xsize, int ysize, double eps, double omega)
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
    } while (iter < MAX_ITER && max_diff > eps);
    return iter;
}