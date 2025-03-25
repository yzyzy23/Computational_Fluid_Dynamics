// =============================================================
// 程序名: HW2
// 功能描述: 计算并比较一阶导数和二阶导数的前向差分、中心差分和精确值，以sin(x)为例
// 输入参数:
// h : 步长参数 (h > 0)，控制差分计算精度
// x : 计算点坐标 (x ∈ ℝ)
// 数学方法:
// 一阶导数差分格式:
// - 前向差分: (f(x + h) - f(x)) / h [O(h)精度]
// - 中心差分: (f(x + h) - f(x - h)) / 2h [O(h²)精度]
// 二阶导数差分格式:
// - 前向差分: (f(x + 2h) - 2f(x + h) + f(x)) / h² [O(h)精度]
// - 中心差分: (f(x + h) - 2f(x) + f(x - h)) / 2h² [O(h²)精度]
// 输出参数：
// 共十二行。
// 前六行为一阶导数的前向差分、中心差分和精确值，以及误差。
// 后六行为二阶导数的前向差分、中心差分和精确值，以及误差。
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
    return T((sin(x + h * 2) - 2 * sin(x + h) + sin(x)) / 2 * (h * h));
}
// 计算一阶导数和二阶导数的中心差分近似值
template <typename T>
T center_1st(T x, T h)
{
    return T((sin(x + h) - sin(x - h)) / (2 * h));
}
template <typename T>
T center_2nd(T x, T h)
{
    return T((sin(x + h) - 2 * sin(x) + sin(x - h)) / 2 * (h * h));
}

int main()
{
    // 输入步长参数和计算点坐标
    double h_d, x_d;
    float h_f, x_f;
    cin >> h_d >> x_d;
    h_f = float(h_d);
    x_f = float(x_d);

    // 计算一阶导数和二阶导数的精确值
    double exact_1st_d = exact_1st(x_d);
    double exact_2nd_d = exact_2nd(x_d);
    float exact_1st_f = exact_1st(x_f);
    float exact_2nd_f = exact_2nd(x_f);

    // 计算一阶导数和二阶导数的前向差分近似值
    double forward_1st_d = forward_1st(x_d, h_d);
    double forward_2nd_d = forward_2nd(x_d, h_d);
    float forward_1st_f = forward_1st(x_f, h_f);
    float forward_2nd_f = forward_2nd(x_f, h_f);

    // 计算一阶导数和二阶导数的中心差分近似值
    double center_1st_d = center_1st(x_d, h_d);
    double center_2nd_d = center_2nd(x_d, h_d);
    float center_1st_f = center_1st(x_f, h_f);
    float center_2nd_f = center_2nd(x_f, h_f);

    // 输出结果
    cout << "First Order Derivative:" << endl;
    cout << "Forward Difference: " << forward_1st_d << " " << forward_1st_f << endl;
    cout << "Center Difference: " << center_1st_d << " " << center_1st_f << endl;
    cout << "Exact Value: " << exact_1st_d << " " << exact_1st_f << endl;
    cout << "Error for Forward Difference: " << abs(exact_1st_d - forward_1st_d) << " " << abs(exact_1st_f - forward_1st_f) << endl;
    cout << "Error for Center Difference: " << abs(exact_1st_d - center_1st_d) << " " << abs(exact_1st_f - center_1st_f) << endl;
    cout << "Second Order Derivative:" << endl;
    cout << "Forward Difference: " << forward_2nd_d << " " << forward_2nd_f << endl;
    cout << "Center Difference: " << center_2nd_d << " " << center_2nd_f << endl;
    cout << "Exact Value: " << exact_2nd_d << " " << exact_2nd_f << endl;
    cout << "Error for Forward Difference: " << abs(exact_2nd_d - forward_2nd_d) << " " << abs(exact_2nd_f - forward_2nd_f) << endl;
    cout << "Error for Center Difference: " << abs(exact_2nd_d - center_2nd_d) << " " << abs(exact_2nd_f - center_2nd_f) << endl;

    return 0;
}