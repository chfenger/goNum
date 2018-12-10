// Secant2P_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-02
版本   : 0.0.0
------------------------------------------------------
    双点弦截法求解方程 f(x)=0 在区间[a, b]内的根
理论：
    1. 当xE[a, b]，f''(x)连续，f'(x) != 0

           xn0*f(xn1) - xn1*f(xn0)
    xn2 = -------------------------
               f(xn1) - f(xn0)

    超线性收敛，收敛阶(1+5^0.5)/2
------------------------------------------------------
输入   :
    fn      f(x)函数，定义为等式左侧部分，右侧为0
    a, b    求解区间
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

func Secant2P(fn func(float64) float64, a, b float64, N int, tol float64) (*float64, bool) {
	/*
		双点弦截法求解方程 f(x)=0 在区间[a, b]内的根
		输入   :
		    fn      f(x)函数，定义为等式左侧部分，右侧为0
		    a, b    求解区间
		    N       步数上限
		    tol     误差上限
		输出   :
		    sol     解值，指针
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var sol float64
	var err bool = false

	//判断a b的次序
	if (b < a) || (fn(a)*fn(b) > 0) {
		return &sol, err
	}
	// 求解
	sol = (a*fn(b) - b*fn(a)) / (fn(b) - fn(a))
	for i := 0; i < N; i++ {
		//判断是否解得
		if (fn(a)*fn(sol) > 0) && (math.Abs(sol-a) < tol) {
			err = true
			return &sol, err
		} else if (fn(a)*fn(sol) < 0) && (math.Abs(sol-b) < tol) {
			err = true
			return &sol, err
		}
		//下一步
		switch {
		case fn(a)*fn(sol) > 0:
			a = sol
		default:
			b = sol
		}
		sol = (a*fn(b) - b*fn(a)) / (fn(b) - fn(a))
	}
	return &sol, err
}

func BenchmarkSecant2P(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Secant2P(func(x float64) float64 { return x*x*x - 2.0*x - 5.0 }, 2.0, 3.0, 1000, 1e-6)
	}
}
