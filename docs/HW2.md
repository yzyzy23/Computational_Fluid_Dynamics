# 1.数理算法原理
f(x)的泰勒展开：
$$f(x+\Delta x)=f(x)+\Delta xf'(x)+\Delta x^2\frac{f''(x)}{2!}+\ldots+\Delta x^n\frac{f^{(n)}(x)}{n!}~~~~~~(1)$$

将(1)式保留到三阶项，可得：
$$f(x+\Delta x)=f(x)+\Delta xf'(x)+\Delta x^2\frac{f''(x)}{2!}+\Delta x^3\frac{f'''(x)}{3!}+O(\Delta x^4)~~~~~~(2)$$
$$f(x-\Delta x)=f(x)-\Delta xf'(x)+\Delta x^2\frac{f''(x)}{2!}-\Delta x^3\frac{f'''(x)}{3!}+O(\Delta x^4)$$
两式相减得到一阶导数的中心差分格式：
$$f'(x)=\frac{f(x+\Delta x)-f(x-\Delta x)}{2\Delta x}+\Delta x^2\frac{f'''(x)}{3!}+O(\Delta x^3)$$
同理对(1)式保留到四阶项，两式相加得到二阶导数的中心差分格式:
$$f''(x)=\frac{f(x+\Delta x)-2f(x)+f(x-\Delta x)}{2\Delta x^2}+\Delta x^2\frac{f''''(x)}{4!}+O(\Delta x^3)$$
从推导中可以看到二者均为二阶精度

\
将(1)式保留到二阶项并移项得到一阶导数的前向差分格式：
$$f'(x)=\frac{f(x+\Delta x)-f(x)}{\Delta x}+\Delta x\frac{f''(x)}{2!}+O(\Delta x^2)$$
将(1)式做替换并保留到三阶项可得到：
$$f(x+2\Delta x)=f(x)+2\Delta xf'(x)+4\Delta x^2\frac{f''(x)}{2!}+8\Delta x^3\frac{f'''(x)}{3!}+O(\Delta x^4)~~~~~~(3)$$
(3) - 2(2)得到二阶导数的前向差分格式：
$$f''(x)=\frac{f(x+2\Delta x)-2f(x+\Delta x)+f(x)}{\Delta x^2}+\Delta xf'''(x)+O(x^2)$$
从推导中可以看到二者均为一阶精度