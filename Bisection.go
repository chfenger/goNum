// Bisection
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-01
版本   : 0.0.0
------------------------------------------------------
    此程序设计使用二分法来求解连续、单自变量、单调函数（区间
内）指定有限区间上的解
    线性收敛
------------------------------------------------------
输入   :
    fn      函数，定义为等式左侧部分，右侧为零
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

func Bisection(fn func(float64) float64, a, b float64, N int, tol float64) (float64, bool) {
	/*
		用二分法来求解连续、单自变量、单调函数（区间内）指定有限区间上的解
		输入   :
		    fn      函数，定义为等式左侧部分，右侧为零
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

	//判断在[a,b]区间是否有解
	if (fn(a) > 0 && fn(b) > 0) || (fn(a) < 0 && fn(b) < 0) {
		return sol, err
	}

	//求解
	for i := 0; i < N; i++ {
		sol = (a + b) / 2
		//解出
		if math.Abs(fn(sol)) < tol {
			err = true
			return sol, err
		}
		//未解出，重置区间边界
		switch {
		case fn(sol) < 0 && fn(a) < 0:
			a = sol
		case fn(sol) > 0 && fn(a) < 0:
			b = sol
		case fn(sol) < 0 && fn(b) < 0:
			b = sol
		case fn(sol) > 0 && fn(b) < 0:
			a = sol
		default:
			return sol, err
		}
	}
	return sol, err
}
