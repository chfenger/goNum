// Secant1P
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-02
版本   : 0.0.0
------------------------------------------------------
    单点弦截法求解方程 f(x)=0 在区间[a, b]内的根
理论：
    1. 当xE[a, b]，f'(x)、f''(x)连续且不变号
    2. 选取初值x0E[a, b]，使f(x0)*f''(x0) > 0，x0选取其中
       一个，则x1选另外一个
    线性收敛
------------------------------------------------------
输入   :
    fn      f(x)函数，定义为等式左侧部分，右侧为0
    fn2     f''(x)函数
    a, b    求解区间
    N       步数上限
    tol     误差上限
输出   :
    sol     解值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// Secant1P 单点弦截法求解方程 f(x)=0 在区间[a, b]内的根
func Secant1P(fn, fn2 func(float64) float64, a, b float64,
	N int, tol float64) (float64, bool) {
	/*
		单点弦截法求解方程 f(x)=0 在区间[a, b]内的根
		输入   :
		    fn      f(x)函数，定义为等式左侧部分，右侧为0
		    fn2     f''(x)函数
		    a, b    求解区间
		    N       步数上限
		    tol     误差上限
		输出   :
		    sol     解值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var sol float64
	var err bool = false

	//判断a b的次序以及选取初值
	if b < a {
		return sol, err
	}
	switch {
	// a
	case fn(a)*fn2(a) > 0:
		for i := 0; i < N; i++ {
			sol = a - (b-a)*fn(a)/(fn(b)-fn(a))
			// 求解成功
			if math.Abs(sol-b) < tol {
				err = true
				return sol, err
			}
			b = sol
		}
		return sol, err
	// b
	case fn(b)*fn2(b) > 0:
		for i := 0; i < N; i++ {
			sol = b - (a-b)*fn(b)/(fn(a)-fn(b))
			// 求解成功
			if math.Abs(sol-a) < tol {
				err = true
				return sol, err
			}
			a = sol
		}
		return sol, err
	}
	return sol, err
}
