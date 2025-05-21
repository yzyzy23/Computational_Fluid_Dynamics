# 1.数理算法原理
顶盖驱动的方腔流动实际上就是要求解二维不可压缩流动：
$$
\begin{cases}
 \frac{\partial\omega}{\partial t}+(v \cdot\nabla)\omega=\nu\Delta\omega&(1.1)\\
 \Delta\psi=-\omega&(1.2)\\
 u=\frac{\partial\psi}{\partial y},~~~~~~v=-\frac{\partial\psi}{\partial x}\\
\end{cases}
$$
边界条件为$u(x,1)=sin^2(\pi x)$，其他边界为固壁。
下面对方程组分别进行离散。
对于式$(1.1)$，首先展开成以下形式：
$$\frac{\partial\omega}{\partial t}+\frac{\partial\psi}{\partial y}\frac{\partial\omega}{\partial x}-\frac{\partial\psi}{\partial x}\frac{\partial\omega}{\partial y}=\nu\Delta\omega$$
对空间导数做中心差分，并将Laplace算子离散为五点Laplace算子，时间导数做向前Euler差分：
$$\frac{\omega^{n+1}_{i,j}-\omega^{n}_{i,j}}{\Delta t}+\frac{\psi^{n}_{i,j+1}-\psi^{n}_{i,j-1}}{2h}\frac{\omega^{n}_{i+1,j}-\omega^{n}_{i-1,j}}{2h}-\frac{\psi^{n}_{i+1,j}-\psi^{n}_{i-1,j}}{2h}\frac{\omega^{n}_{i,j+1}-\omega^{n}_{i,j-1}}{2h}=\nu\frac{\omega^{n}_{i,j+1}+\omega^{n}_{i,j-1}+\omega^{n}_{i+1,j}+\omega^{n}_{i-1,j}-4\omega^{n}_{i,j}}{h^2}$$
从中可以导出内点的$\omega$迭代式：
$$\omega^{n+1}_{i,j}=\Delta t*(-(\frac{\psi^{n}_{i,j+1}-\psi^{n}_{i,j-1}}{2h})(\frac{\omega^{n}_{i+1,j}-\omega^{n}_{i-1,j}}{2h})+(\frac{\psi^{n}_{i+1,j}-\psi^{n}_{i-1,j}}{2h})(\frac{\omega^{n}_{i,j+1}-\omega^{n}_{i,j-1}}{2h})+\nu\frac{\omega^{n}_{i,j+1}+\omega^{n}_{i,j-1}+\omega^{n}_{i+1,j}+\omega^{n}_{i-1,j}-4\omega^{n}_{i,j}}{h^2}+\omega^{n}_{i,j})$$
对于式$(1.2)$，等价为解Poisson方程，可直接调用HW4中的SOR法进行求解。
对于边界上的涡量，一般形式的Thom公式为$\omega_0=-\frac{2(\psi_1-\psi_0+vh)}{h^2}$
带入得到边界条件：
$$\omega_{0,j}=\frac{2(\psi_{0,j}-\psi_{1,j})}{h^2}-\frac{2}{h}v_{0,j}$$
$$\omega_{n-1,j}=\frac{2(\psi_{n-1,j}-\psi_{n-2,j})}{h^2}+\frac{2}{h}v_{n-1,j}$$
$$\omega_{i,0}=\frac{2(\psi_{i,0}-\psi_{i,1})}{h^2}+\frac{2}{h}u_{i,0}$$
$$\omega_{i,m-1}=\frac{2(\psi_{i,m-1}-\psi_{i,m-2})}{h^2}-\frac{2}{h}u_{i,m-1}$$
最后，这里给出的初始条件是关于速度的，可以根据涡量的定义式$\omega=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}$将给定的速度条件转化为涡量的初值。
于是得到算法的主要步骤如下：
1.从涡量定义式将给出的速度关系转化为初始涡量
2.调用SOR得到初始流函数
3.用Thom公式得到边界涡量
4.迭代内点涡量
5.调用SOR得到当前流函数
6.重复3 4 5
# 2.代码生成和调试
在代码中实现如下函数：
1.InitVelocity：用给定的速度进行初始化
2.InitVorticity：用速度对涡量进行初始化
3.SOR：通过超松弛迭代得到各点的流函数值，若得到迭代前后流函数值的最大改变量小于阈值，则认为收敛，特别的，若迭代次数为1次，即迭代1次就收敛，则我们认为流动达到稳定状态，可以结束整个程序。最后返回bool类型变量，用于判断是否结束整个程序。
4.UpdateVorticity：用流函数值更新内点的涡量
5.ApplyBoundCondition：用Thom公式更新边界的涡量
最终核心部分如下：

    InitVelocity(u, v, xsize, ysize);
    InitVorticity(u, v, omega, xsize, ysize, h);
    is_converge = SOR(psi, omega, xsize, ysize, h, relax_factor);

    while (is_converge)
    {
        ApplyBoundCondition(psi, omega, u, v, xsize, ysize, h);
        UpdateVorticity(omega, psi, h, nu, dt, xsize, ysize);
        is_converge = SOR(psi, omega, xsize, ysize, h, relax_factor);
    }

程序输入h，即指定格点的长度，得到结果为两个文件StreamFunction.txt和Vorticity.txt，文件每行三个数i j data，分别为离散化之后的格点坐标和对应的函数值。
接下来用Visualize.py对两个文件进行可视化，要注意python中画图的坐标系对二维数组来说是列优先而非行优先，因此在读入数据时需要反转坐标。