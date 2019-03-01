// ODEHeun_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-26
版本   : 0.0.0
------------------------------------------------------
    常微分方程的Heun解法
理论：
    对于常微分方程
     dy
    ---- = f(x, y)
     dx
    y(x0) = y0, x0 <= x

    Heun法为
    1. p_(k+1) = yk+hf(xk,yk)   //欧拉法
    2. y_(k+1) = yk+h(f(xk,yk)+f(x_(k+1),p_(k+1))/2   //梯形法
       k = 0,1,2,3,...

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 9.3
------------------------------------------------------
输入   :
    fun     被积分函数
    x0, y0  初值
    h       步长
    n       迭代次数
输出   :
    sol     解矩阵，nx2
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// ODEHeun 常微分方程的Heun解法
func ODEHeun(fun func(float64, float64) float64, x0, y0, h float64, n int) (goNum.Matrix, bool) {
	/*
		常微分方程的Heun解法
		输入   :
		    fun     被积分函数
		    x0, y0  初值
		    h       步长
		    n       迭代次数
		输出   :
		    sol     解矩阵，nx2
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if n < 0 {
		panic("Error in goNum.ODEHeun: n is not a positive value")
	}

	sol := goNum.ZeroMatrix(n+1, 2)
	p := goNum.ZeroMatrix(n+1, 2)
	var err bool = false

	//初值
	sol.SetMatrix(0, 0, x0)
	sol.SetMatrix(0, 1, y0)

	for i := 1; i < n+1; i++ {
		p.SetMatrix(i, 0, sol.GetFromMatrix(i-1, 0)+h)   //xi=x_(i-1)+h
		sol.SetMatrix(i, 0, sol.GetFromMatrix(i-1, 0)+h) //xi=x_(i-1)+h

		soltemp := fun(sol.GetFromMatrix(i-1, 0), sol.GetFromMatrix(i-1, 1))
		p.SetMatrix(i, 1, sol.GetFromMatrix(i-1, 1)+h*soltemp)

		soltemp = h * (soltemp + fun(sol.GetFromMatrix(i, 0), p.GetFromMatrix(i, 1))) / 2.0
		sol.SetMatrix(i, 1, sol.GetFromMatrix(i-1, 1)+soltemp)
	}

	err = true
	return sol, err
}

func fun53(x, y float64) float64 {
	return (x - y) / 2.0
}

func BenchmarkODEHeun(b *testing.B) {
	for i := 0; i < b.N; i++ {
		goNum.ODEHeun(fun53, 0.0, 1.0, 1e-3, 3e3) //解析值y(3.0)=1.669390
	}
}
