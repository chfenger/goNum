// NewtonIterate_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-01
版本   : 0.0.0
------------------------------------------------------
    牛顿迭代求解非线性方程 f(x)=0 在区间[a, b]内的根
理论：
    （局部收敛定律）
    1. f(x)在区间[a, b]具有二阶连续导数；
    2. 当xE[a, b]，f'(x) != 0；
    （非局部收敛定律）
    1. 当xE[a, b]，f'(x)、f''(x)连续且不变号
    2. 选取初值x0E[a, b]，使f(x0)*f''(x0) > 0
    平方收敛
------------------------------------------------------
输入   :
    fn      f(x)函数，定义为等式左侧部分，右侧为0
    fn1     f'(x)函数
    a, b    求解区间
    c       求解初值
    N       步数上限
    tol     误差上限
输出   :
    sol     解值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"
)

func NewtonIterate(fn, fn1 func(float64) float64, a, b, c float64, N int, tol float64) (float64, bool) {
	/*
		牛顿迭代求解非线性方程 f(x)=0 在区间[a, b]内的根
		输入   :
		    fn      f(x)函数，定义为等式左侧部分，右侧为0
		    fn1     f'(x)函数
		    a, b    求解区间
		    c       求解初值
		    N       步数上限
		    tol     误差上限
		输出   :
		    sol     解值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var sol float64
	var err bool = false

	// 判断端点和初值是否为所求之解
	switch {
	case math.Abs(fn(a)) < tol:
		sol = a
		err = true
		return sol, err
	case math.Abs(fn(b)) < tol:
		sol = b
		err = true
		return sol, err
	case math.Abs(fn(c)) < tol:
		sol = c
		err = true
		return sol, err
	}

	//求解
	sol = c - fn(c)/fn1(c)
	for i := 0; i < N; i++ {
		if math.Abs(sol-c) < tol {
			err = true
			return sol, err
		}
		c = sol
		sol = c - fn(c)/fn1(c)
	}
	return sol, err
}

func BenchmarkNewtonIterate(b *testing.B) {
	for i := 0; i < b.N; i++ {
		NewtonIterate(func(x float64) float64 { return math.Pow(x, 3.0) - math.Pow(x, 2.0) - 1.0 }, func(x float64) float64 { return 3.0*math.Pow(x, 2.0) - 2.0*x }, 1.4, 1.5, 1.5, 1000, 1e-6)
	}
}
