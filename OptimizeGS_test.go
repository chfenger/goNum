// OptimizeGS_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-24
版本   : 0.0.0
------------------------------------------------------
    黄金分割法(Golden Section)求单峰单自变量极小值
理论：
    对于在区间[a, b]内有定义的凹函数f(x)，取黄金分割点：
    c = a+(1-r)(b-a)
    d = b-(1-r)(b-a)
    其中r为黄金分割比例(Sqrt(5)-1)/2

    如果f(c) <= f(d)，则将d赋予b，继续迭代；
    如果f(c) > f(d)，则将c赋予a，继续迭代。
    迭代终止条件为Abs(f(a)-f(b)) < tol，取小值（c或d）

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 8.1.1.1
------------------------------------------------------
输入   :
    fun     函数
    a, b    区间范围
    tol     控制误差
    N       最大迭代步数
输出   :
    sol     解
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

func OptimizeGS(fun func(float64) float64, a, b, tol float64, N int) (float64, bool) {
	/*
		黄金分割法(Golden Section)求单峰单自变量极小值
		输入   :
		    fun     函数
		    a, b    区间范围
		    tol     控制误差
		    N       最大迭代步数
		输出   :
		    sol     解
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断a和b的关系
	if math.Abs(fun(a)-fun(b)) < tol {
		if fun(a) < fun(b) {
			return a, true
		} else {
			return b, true
		}
	}

	var sol float64
	var err bool = false
	r1 := 1.0 - (math.Sqrt(5.0)-1.0)/2.0 //1-r

	for i := 0; i < N; i++ {
		ba := b - a //b-a
		c := a + r1*ba
		d := b - r1*ba
		//区间压缩
		if fun(c) > fun(d) {
			a = c
		} else { //fun(c)<=fun(d)
			b = d
		}
		//误差判断
		if math.Abs(fun(a)-fun(b)) < tol {
			err = true
			if fun(c) < fun(d) {
				sol = c
			} else {
				sol = d
			}
			return sol, err
		}
	}

	return sol, err
}

func fun50(x float64) float64 {
	return x*x - math.Sin(x)
}

func BenchmarkOptimizeGS(b *testing.B) {
	for i := 0; i < b.N; i++ {
		goNum.OptimizeGS(fun50, 0.0, 1.0, 1e-10, 1e3)
	}
}
