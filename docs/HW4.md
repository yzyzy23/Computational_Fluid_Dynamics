# 1.数理算法原理
流体的热传导方程：
$$\frac{\partial T}{\partial t}-\alpha\nabla^2T=0$$
最终达到稳态时$\frac{\partial T}{\partial t}=0$，因此在二维条件下原方程化为：
$$\frac{\partial^2T}{\partial^2x}+\frac{\partial^2T}{\partial^2y}=0$$
其中温度分布$T=T(x,y)$（温度单位：K，长度单位：cm）
边界条件为：
$$T(x,0)=T(0,y)=T(15,y)=293.15$$
$$T(x,12)=373.15$$
$$0\leq x\leq 15, 0\leq y\leq 12$$
要用迭代法求解以上方程，可以先用五点差分格式进行离散：
$$\frac{T(i+1,j)-2T(i,j)+T(i-1,j)}{\Delta x^2}+\frac{T(i,j+1)-2T(i,j)+T(i,j-1)}{\Delta y^2}=0$$
整理得：
$$T(i,j)=\frac{\Delta y^2(T(i+1,j)+T(i-1,j))+\Delta x^2(T(i,j+1)+T(i,j-1))}{2(\Delta x^2+\Delta y^2)}$$
假设采用均匀网格，即$\Delta x=\Delta y$,则上式写为：
$$T(i,j)=\frac{T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)}{4}$$
Gauss-Seidel迭代法为:
$$T^{k+1}(i,j)=\frac{T^k(i+1,j)+T^{k+1}(i-1,j)+T^k(i,j+1)+T^{k+1}(i,j-1)}{4}$$
松弛法要求：$T^{k+1}(i,j)=T^k(i,j)+\omega\Delta T$，这里$\omega$是松弛因子
下面计算$\Delta T$：
$$\Delta T=T^{k+1}(i,j)-T^k(i,j)=\frac{T^k(i+1,j)+T^{k+1}(i-1,j)+T^k(i,j+1)+T^{k+1}(i,j-1)}{4}-T^k(i,j)$$
带入原方程并整理得到：
$$T^{k+1}(i,j)=(1-\omega)T^k(i,j)+\frac{\omega}{4}[T^{k+1}(i-1,j)+T^{k+1}(i,j-1)+T^k(i+1,j)+T^k(i,j+1)]$$
这就是松弛法的迭代公式，终止条件为$max|T^{k+1}-T^k|<\epsilon$，下面取$\epsilon=10^{-5}$求解本问题。