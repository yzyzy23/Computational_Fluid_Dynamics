# 1.数理原理推导
对于一维的Sod激波管问题，实际上相当于求解一维的欧拉方程
$$\frac{\partial U}{\partial t}+\frac{\partial f}{\partial x}=0$$
初始条件为
$$\rho,u,p=1,0,1~~(x<0)$$
$$\rho,u,p=0.125,0,0.1~~(x>0)$$
## 1.精确解介绍
首先我们介绍求解Riemann精确解的方法及其结果。
从流体力学课上讲过的知识我们可以知道，激波管中可能出现三种波：激波、接触间断和膨胀波，其分别对应不同的边界条件。接下来，我们可以写出质量通量、动量通量和能量通量守恒的条件，并取控制体微元进行求解。这里查阅文档可知Sod激波管的实际演化为右边为激波，左边为膨胀波，因此我们可以得到方程如下：
左边膨胀波满足等熵条件：
$$\frac{p_{left}}{\rho_{left}^{\gamma}}=\frac{p}{\rho^{\gamma}}$$
$$u_{left}+\frac{2c_{left}}{\gamma -1}=u+\frac{2c}{\gamma -1}$$
右边激波满足Rankine-Hugoniot关系式：
$$\rho(u-Z)=\rho_{right}(u_{right}-Z)$$
$$\rho u(u-Z)+p=\rho_{right}u_{right}(u_{right}-Z)+p_{right}$$
$$E(u-Z)+up=E_{right}(u_{right}-Z)+u_{right}p_{right}$$
以上共5个方程5个未知数，因此可以用Newton迭代法解出。
在实际代码实现方面，我们选择用网上开源的Python程序进行求解，并保证代码结果有较高的可信度。
得到精确解的图形如下：
![精确解](../photo/Sod/ExactSolve.png)
## 2.近似解求法
接下来，我们介绍求近似解的方法。
要求解一维欧拉方程，现有的技术一般分为流通矢量分裂(Flux
 Vector Splitting, FVS) 与流通差分分裂 (Flux Vector Splitting, FVS) 两种方法
### 1.流通矢量分裂（FVS）
方程的原始守恒变量U和流通矢量f(U)：
$$U=(\rho,\rho u,E)^T$$
记为$\omega_i$
则$$f(U)=(\omega_2,\frac{\omega_2^2}{\omega_1}+p,\frac{\omega_2}{\omega_1}(\omega_3+p))^T$$
再利用完全气体状态方程可以将f(U)做进一步的变形，于是可以得到Jacobi矩阵A：
$$A=\frac{\partial f}{\partial U}$$
$$A(U)=\begin{matrix}
    0 & 1 & 0 \\
    \frac{\gamma-3}{2}u^2 & (3-\gamma)u & \gamma-1\\
    -\frac{\gamma pu^2}{(\gamma-1)\rho}+\frac{\gamma-2}{2}u^3 & \frac{\gamma p}{(\gamma-1)\rho}-(\gamma-\frac{3}{2})u^2 & \gamma u
\end{matrix}$$
在引入了声速$c^2=\frac{\gamma p}{\rho}$的前提下，矩阵A的特征值和特征向量可解，其特征值为u-c,u,u+c
于是，要求的$\frac{\partial f}{\partial x}=\frac{\partial AU}{\partial x}=A\frac{\partial U}{\partial x}$，可以用特征值来求解。
下面我们对f(U)进行分裂：
$$\lambda = \lambda^++\lambda^-$$
则前向、后向的通量可求。
在本项目中，我们采用 Steger-Warming (S-W) 分裂法，其要求：
$$\lambda=\frac{\lambda+|\lambda|}{2}$$
于是我们可以将f化为$\lambda_i$的函数，也就得到了分裂后的流通矢量。
### 2.流通差分分裂（FDS）
流通差分分裂主要包括以下步骤：
1.用各种差分格式计算在格点i+1/2处的左右守恒通量
2.用近似解法计算每个格点i+1/2处的F
3.用差分法得到$\frac{\partial f}{\partial x}$
本项目中我们用Roe格式实现上述步骤。
上面我们已经推导了将f(U)转化成AU，进而用特征值进行计算来降低计算量的方法。这里我们对上面的步骤进行改造，对于非线性的Jacobi矩阵A，如果我们能够将其近似为线性的，那么也就化为求解方程$\frac{\partial U}{\partial t}+a\frac{\partial U}{\partial x}=0$，而这是我们所熟悉的方程。
注意到这里的a实际上是导数，因此我们可以用区间平均变化率来进行近似，且由Largrange中值定理，区间内必有这样的点（记为$u_{roe}$）使得其变化率满足要求。
于是近似的Riemann解为：
$$F=\frac{F(U_L)+F(U_R)}{2}-\frac{1}{2}S^{-1}|\Gamma|S(U_R-U_L)$$
这里$S^{-1}\Gamma S$为A的特征值分解。
由Roe公式，我们可以把A(U)写成A($\bar{U}$)，其中$\bar{U}$满足：
$$\bar{\rho}=(\frac{\sqrt{\rho_l}+\sqrt{\rho_r}}{2})^2$$
$$\bar{u}=\frac{\sqrt{\rho_l}u_l+\sqrt{\rho_r}u_r}{\sqrt{\rho_l}+\sqrt{\rho_r}}$$
$$\bar{H}=\frac{\sqrt{\rho_l}H_l+\sqrt{\rho_r}H_r}{\sqrt{\rho_l}+\sqrt{\rho_r}}$$
于是，我们就可以计算F
综上，Roe法的计算步骤如下：
1.用各种差分格式计算在格点i+1/2处的左右守恒通量U
2.用Roe平均计算$\bar{U}$
3.对得到的A进行特征值分解（参考我们在FVS中得到的结论）
4.计算每个格点i+1/2处的F
5.用差分法得到$\frac{\partial f}{\partial x}$
这里需要注意的是，Roe格式本质上是一种迎风格式，而迎风格式本身数值耗散可能会使得膨胀波等被耗散掉，而在特征值为0左右，耗散会变成0，因此我们需要在0处认为增加耗散，也就是熵修正。
这里进行的熵修正为：
$$\lambda=|\lambda|,~~~~~~|\lambda|>\epsilon$$
$$\lambda=\epsilon,~~~~~~|\lambda|<\epsilon$$
### 3.激波捕捉格式
在上面的步骤分析中，我们发现剩下的关键步骤是用差分格式进行激波捕捉，常用的方法有TVD格式，WENO格式和GVC格式。下面我们分别进行介绍。

1.TVD格式

本格式是对现有的差分进行改造，使得其总变差不增
也就是将格式构造为如下形式：
$$u^{n+1}_i=u^n_i+A(u^n_{i+1}-u^n_i)-B(u^n_i-u^n_{i-1})$$
一种TVD的构造格式是
$$u_{j+mid}=u_{j+mid}^l=u_j+\phi(r_j^l)\frac{u_{j+1}-u_j}{2}$$
$$u_{j-mid}=u_{j-mid}^r=u_j-\phi(r_j^r)\frac{u_{j}-u_{j-1}}{2}$$
这里$r_j=\frac{u_{j+1}-u_j}{u_j-u_{j-1}}$，$\phi(r_j^l)$为限制器，我们这里选择Van Leer限制器，其算法为：
$$\phi(r)=\frac{|r|+r}{|r|+1}$$
该限制器具有良好的光滑性。

2.WENO格式

WENO格式（加权本质无振荡方法）的基本思想是通过对多个插值区域进行加权平均，以提高数值通量的精度。具体要求通过为每个插值区域分配权重，确保实现更高的精度。
本项目中我们采用五点WENO格式，其计算公式如下（为了方便起见这里只列出正通量）：
$$IS_0 = \frac{13}{12} (f_{i-2} - 2f_{i-1} + f_i)^2 + \frac{1}{4} (f_{i-2} - 4f_{i-1} + 3f_i)^2$$
$$IS_1 = \frac{13}{12} (f_{i-1} - 2f_i + f_{i+1})^2 + \frac{1}{4} (f_{i-1} - f_{i+1})^2$$
$$IS_2 = \frac{13}{12} (f_i - 2f_{i+1} + f_{i+2})^2 + \frac{1}{4} (3f_i - 4f_{i+1} + f_{i+2})^2$$
$$\alpha_k=\frac{C_k}{(\epsilon+IS_k)^2}$$
$$\omega_k=\frac{\alpha_k}{\alpha_0+\alpha_1+\alpha_2}$$
$$f^{(0)}_{j+\frac{1}{2}} = \frac{1}{3}f_{j-2} - \frac{7}{6}f_{j-1} + \frac{11}{6}f_j$$
$$f^{(1)}_{j+\frac{1}{2}} = -\frac{1}{6}f_{j-1} + \frac{5}{6}f_j + \frac{1}{3}f_{j+1}$$
$$f^{(2)}_{j+\frac{1}{2}} = \frac{1}{3}f_j + \frac{5}{6}f_{j+1} - \frac{1}{6}f_{j+2}$$
$$f_{j+\frac{1}{2}}=\omega_0f^{(0)}_{j+\frac{1}{2}}+\omega_1f^{(1)}_{j+\frac{1}{2}}+\omega_2f^{(2)}_{j+\frac{1}{2}}$$
其中C0=0.1，C1=0.6,C2=0.3。于是我们可以构建差分：
$$\frac{\partial f_j}{\partial x}=\frac{f_{j+\frac{1}{2}}-f_{j-\frac{1}{2}}}{\Delta x}$$

3.GVC格式

GVC格式是为了解决激波附近的数值解因频散误差而出现非物理振荡而提出的，其原理为在在激波前使用慢格式（抑制超前扰动），激波后使用快格式（加速滞后扰动），从而达到抑制振荡的目的。
这里我们选用二阶迎风格式作为快格式，二阶中心差分格式作为慢格式。

### 4.时间推进
在有了以上的步骤后，我们最终的工作是求解形如：
$$\frac{\partial U}{\partial t}=Q(U)$$
的方程。以上方程已经被我们化简为常微分形式，因此我们可以采用计算常微分方程的Runge-Kutta法来进行计算。
本项目中我们选择3阶的Runge-Kutta方法，其公式如下：
$$U_1=U^n+\Delta tQ(U^n)$$
$$U_2=\frac{3}{4}U^n+\frac{1}{4}(U_1+\Delta tQ(U_1))$$
$$U^{n+1}=\frac{1}{3}U^n+\frac{2}{3}(U_2+\Delta tQ(U_2))$$

# 2.代码生成和调试
本项目形式比较分散，代码数量多而繁杂，因此我们采取如下代码生成方法：
在主程序中仅保留初始化和结果输出等内容，将上面的核心内容分散到三个程序CalcDiff.cpp CalcFlux.cpp TimeIter.cpp中，这样代码的可读性和可维护性均得到较大提升。
最后，使用python程序SodVisualize.py对结果进行可视化。