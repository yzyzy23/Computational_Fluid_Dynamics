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
对于式$(1.1)$，首先展开成以下形式：
$$\frac{\partial\omega}{\partial t}+\frac{\partial\psi}{\partial y}\frac{\partial\omega}{\partial x}-\frac{\partial\psi}{\partial x}\frac{\partial\omega}{\partial y}=\nu\Delta\omega$$
对空间导数做中心差分，并将Laplace算子离散为五点Laplace算子，时间导数做向前Euler差分：
$$\omega^{n+1}_{i,j}-\omega^{n}_{i,j}+\frac{\psi^{n}_{i,j+1}-\psi^{n}_{i,j-1}}{2h}\frac{\omega^{n}_{i+1,j}-\omega^{n}_{i-1,j}}{2h}-\frac{\psi^{n}_{i+1,j}-\psi^{n}_{i-1,j}}{2h}\frac{\omega^{n}_{i,j+1}-\omega^{n}_{i,j-1}}{2h}=\nu\frac{\omega^{n}_{i,j+1}+\omega^{n}_{i,j-1}+\omega^{n}_{i+1,j}+\omega^{n}_{i-1,j}-4\omega^{n}_{i,j}}{h^2}$$
从中可以导出内点的$\omega$迭代式：
$$\omega^{n+1}_{i,j}=-(\frac{\psi^{n}_{i,j+1}-\psi^{n}_{i,j-1}}{2h})(\frac{\omega^{n}_{i+1,j}-\omega^{n}_{i-1,j}}{2h})+(\frac{\psi^{n}_{i+1,j}-\psi^{n}_{i-1,j}}{2h})(\frac{\omega^{n}_{i,j+1}-\omega^{n}_{i,j-1}}{2h})+\nu\frac{\omega^{n}_{i,j+1}+\omega^{n}_{i,j-1}+\omega^{n}_{i+1,j}+\omega^{n}_{i-1,j}-4\omega^{n}_{i,j}}{h^2}+\omega^{n}_{i,j}$$
对于式$(1.2)$，等价为解Poisson方程
对于边界上的涡量，一般形式的Thom公式为$\omega_0=-\frac{2(\psi_1-\psi_0+vh)}{h^2}$
带入得到边界条件：
$$\omega_{0,j}=\frac{2(\psi_{0,j}-\psi_{1,j})}{h^2}-\frac{2}{h}v_{0,j}$$
$$\omega_{n-1,j}=\frac{2(\psi_{n-1,j}-\psi_{n-2,j})}{h^2}+\frac{2}{h}v_{n-1,j}$$
$$\omega_{i,0}=\frac{2(\psi_{i,0}-\psi_{i,1})}{h^2}+\frac{2}{h}u_{i,0}$$
$$\omega_{i,m-1}=\frac{2(\psi_{i,m-1}-\psi_{i,m-2})}{h^2}-\frac{2}{h}u_{i,m-1}$$
对于边界条件，可以根据涡量的定义式$\omega=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}$从给定的速度条件中导出。