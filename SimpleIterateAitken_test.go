// SimpleIterateAitken_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-01
版本   : 0.0.0
------------------------------------------------------
    简单迭代求解类x=g(x)方程的解 xn+1=g(xn)
理论：
    1. g(x)在区间[a, b]可导；
    2. 当xE[a, b]，g(x)E[a, b]；
    3. 对于任意xE[a, b]，|g‘(x)| <= L < 1
    线性收敛

    则求解所得的根xn与真实根xr的的误差：
                L^n
    |xn-xr| <= ----- |x1-x0|
                1-L
Aitken boost method
    xn1 = g(xn)
    xn2 = g(xn1)
            xn2*xn - xn1^2
    xn1 = ------------------ (n = 0, 1, 2,...)
           xn2 - 2*xn1 + xn
------------------------------------------------------
输入   :
    fn      g(x)函数，定义为等式右侧部分，左侧为x
    a, b    求解区间
    c       求解初值
    N       步数上限
    tol     误差上限
输出   :
    sol     解值，指针
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"
)

func SimpleIterateAitken(fn func(float64) float64, a, b, c float64, N int, tol float64) (*float64, bool) {
	/*
		简单迭代求解类x=g(x)方程的解 xn+1=g(xn)
		输入   :
		    fn      g(x)函数，定义为等式右侧部分，左侧为x
		    a, b    求解区间
		    c       求解初值
		    N       步数上限
		    tol     误差上限
		输出   :
		    sol     解值，指针
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var sol, n2 float64
	var err bool = false

	// 判断端点和初值是否为所求之解
	switch {
	case math.Abs(fn(a)-a) < tol:
		sol = a
		err = true
		return &sol, err
	case math.Abs(fn(b)-b) < tol:
		sol = b
		err = true
		return &sol, err
	case math.Abs(fn(c)-c) < tol:
		sol = c
		err = true
		return &sol, err
	}

	//求解
	sol = fn(c)
	n2 = fn(sol)
	sol = (n2*c - sol*sol) / (n2 - 2.0*sol + c)
	for i := 0; i < N; i++ {
		if (math.Abs(sol - c)) < tol {
			err = true
			return &sol, err
		}
		c = sol
		sol = fn(c)
		n2 = fn(sol)
		sol = (n2*c - sol*sol) / (n2 - 2.0*sol + c)
	}
	return &sol, err
}

func BenchmarkSimpleIterateAitken(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SimpleIterate(func(x float64) float64 { return math.Pow(math.Pow(x, 2.0)+1.0, 1.0/3.0) }, 1.4, 1.5, 1.5, 1000, 1e-6)
	}
}
