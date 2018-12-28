// ODEDiff_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-26
版本   : 0.0.0
------------------------------------------------------
    差分方法求解常微分方程
理论：
    对于常微分方程：x''(t) = p(t)x'(t)+q(t)(t)+r(t)
    机器边值x(a) = x0, x(b) = xN
    使用中心差分公式可得
     x_(j+1)-2xj+x_(j-1)       x_(j+1)-x_(j-1)
    --------------------- = pj----------------- + qj*xj+rj
            h^2                      2h
    即
      -h                                h
    (---pj-1)x_(j-1) + (2+h^2*qj)xj + (---pj-1)x_(j+1) = -h^2*rj
      2                                 2
    下标j表示*(tj), tj=a+j*h，(区间[a, b]等分为N等份)
    整理成N-1阶线性方程组:
    |2+h^2*q1  h*p1/2-1                                                      |
    |-h*p2/2-1 2+h^2*q2  h*p2/2-1                                            |
    |          -h*p3/2-1 2+h^2*q3 h*p3/2-1                                   |*
    |                    ......                                              |
    |                             -h*p_(N-2)/2-1 2+h^2*q_(N-2)  h*p_(N-2)/2-1|
    |                                            -h*p_(N-1)/2-1 2+h^2*q_(N-1)|

    [x1 x2 x3 ... x_(N-1)]' =

    [-h^2*r1+e0 -h^2*r2 -h^2*r3 ... -h^2*r_(N-2) -h^2*r_(N-1)+eN]'
    其中：
    e0 = (h*p1/2+1)x0, eN = (h*p_(N-1)/2+1)xN

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 9.9
------------------------------------------------------
输入   :
    funp, funq, funr 被积分函数系数
    x0      初值,2x2, 按列a, b
    Nn      积分步数
输出   :
    sol     解矩阵, 2x(Nn+1)
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

func ODEDiff(funp, funq, funr func(float64) float64, x0 goNum.Matrix, Nn int) (goNum.Matrix, bool) {
	/*
		差分方法求解常微分方程
		输入   :
		    funp, funq, funr 被积分函数系数
		    x0      初值,2x2, 按列a, b
		    Nn      积分步数
		输出   :
		    sol     解矩阵, 2x(Nn+1)
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断x0维数
	if (x0.Rows != 2) || (x0.Columns != 2) {
		panic("Error in goNum.ODEDiff: Initial values error")
	}
	//判断Nn
	if Nn < 1 {
		panic("Error in goNum.ODEDiff: Nn must greater than zero")
	} else if Nn == 1 {
		return x0, true
	}

	sol := goNum.ZeroMatrix(2, Nn+1)
	var err bool = false
	Aa := goNum.ZeroMatrix(Nn-1, Nn-1)
	Bb := goNum.ZeroMatrix(Nn-1, 1)
	h := (x0.GetFromMatrix(0, 1) - x0.GetFromMatrix(0, 0)) / float64(Nn)

	//ti
	for i := 0; i < Nn+1; i++ {
		sol.SetMatrix(0, i, x0.GetFromMatrix(0, 0)+h*float64(i))
	}
	//x0, xN
	sol.SetMatrix(1, 0, x0.GetFromMatrix(1, 0))
	sol.SetMatrix(1, Nn, x0.GetFromMatrix(1, 1))

	//第一行
	Aa.SetMatrix(0, 0, 2.0+h*h*funq(sol.GetFromMatrix(0, 1)))
	Aa.SetMatrix(0, 1, h*funp(sol.GetFromMatrix(0, 1))/2.0-1.0)
	e0 := (h*funp(sol.GetFromMatrix(0, 1))/2.0 + 1.0) * sol.GetFromMatrix(1, 0)
	Bb.SetMatrix(0, 0, -1.0*h*h*funr(sol.GetFromMatrix(0, 1))+e0)
	for i := 1; i < Nn-2; i++ {
		Aa.SetMatrix(i, i-1, -1.0*h*funp(sol.GetFromMatrix(0, i+1))/2.0-1.0)
		Aa.SetMatrix(i, i, 2.0+h*h*funq(sol.GetFromMatrix(0, i+1)))
		Aa.SetMatrix(i, i+1, h*funp(sol.GetFromMatrix(0, i+1))/2.0-1.0)
		Bb.SetMatrix(i, 0, -1.0*h*h*funr(sol.GetFromMatrix(0, i+1)))
	}
	//最后行
	Aa.SetMatrix(Nn-2, Nn-2-1, -1.0*h*funp(sol.GetFromMatrix(0, Nn-1))/2.0-1.0)
	Aa.SetMatrix(Nn-2, Nn-2, 2.0+h*h*funq(sol.GetFromMatrix(0, Nn-1)))
	eN := (-1.0*h*funp(sol.GetFromMatrix(0, Nn-1))/2.0 + 1.0) * sol.GetFromMatrix(1, Nn)
	Bb.SetMatrix(Nn-2, 0, -1.0*h*h*funr(sol.GetFromMatrix(0, Nn-1))+eN)

	//求解线性方程组LEs_Chasing
	xTemp, errTemp := goNum.LEs_Chasing(Aa, Bb)
	if errTemp != true {
		panic("Error in goNum.ODEDiff: Solve error")
	}

	//xTemp赋予sol
	for i := 1; i < Nn; i++ {
		sol.SetMatrix(1, i, xTemp.GetFromMatrix(i-1, 0))
	}

	err = true
	return sol, err
}

func fun54p(t float64) float64 {
	return 2.0 * t / (1.0 + t*t)
}

func fun54q(t float64) float64 {
	return -2.0 / (1.0 + t*t)
}

func fun54r(t float64) float64 {
	return 1.0
}

func BenchmarkODEDiff(b *testing.B) {
	x54 := goNum.NewMatrix(2, 2, []float64{0.0, 4.0, 1.25, -0.95})
	for i := 0; i < b.N; i++ {
		goNum.ODEDiff(fun54p, fun54q, fun54r, x54, 160)
	}
}
