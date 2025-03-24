// =============================================================
// 函数名: HW2
// 功能描述: 计算正弦函数f(x)=sin(x)在指定点处的数值导数
// 输入参数:
// h : 步长参数 (h > 0)，控制差分计算精度
// x : 计算点坐标 (x ∈ ℝ)
// 数学方法:
// 一阶导数差分格式:
// - 前向差分: (f(x + h) - f(x)) / h [O(h)精度]
//
// 二阶导数差分格式:
// - 前向差分: (f(x + 2h) - 2f(x + h) + f(x)) / h² [O(h)精度]
// 输出参数：
// 共四行。
// 第一行为前向差分和中心差分计算的一阶导数值（单、双精度均包括在内，下同），第二行为一阶导数值的精确值以及误差。
// 第三行为前向差分和中心差分计算的二阶导数值，第四行为二阶导数值的精确值以及误差。
// =============================================================
#include <iostream>
#include <cmath>

using namespace std;

// 以下用Template实现，方便改变数据类型
// 计算一阶导数和二阶导数的精确值
template <typename T>
T exact_1st(T x)
{
    return T(cos(x));
}
template <typename T>
T exact_2nd(T x)
{
    return T(-sin(x));
}
// 计算一阶导数和二阶导数的前向差分近似值
template <typename T>
T forward_1st(T x, T h)
{
    return T((sin(x + h) - sin(x)) / h);
}
template <typename T>
T forward_2nd(T x, T h)
{
    return T((sin(x + h * 2) - 2 * sin(x + h) + sin(x)) / (h * h));
}

int main()
{
    // 输入步长参数和计算点坐标
    double h_d, x_d;
    float h_f, x_f;
    cin >> h_d >> x_d;
    h_f = float(h_d);
    x_f = float(x_d);

    return 0;
}