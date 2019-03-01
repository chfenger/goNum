// FittingBezier_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-23
版本   : 0.0.0
------------------------------------------------------
    Bezier曲线拟合控制点
理论：
    给定控制点集(xi, yi), i=0,1,...,N
    则Bezier曲线可以表示为：
    |        N
    |x(t) = Sum xi*B_(i,N)(t)
    |       i=0
    |
    |        N
    |y(t) = Sum yi*B_(i,N)(t)
    |       i=0
    其中，
    B_(i,N)(t)为Bernstein多项式：
                  N-i
    B_(i,N)(t) = C   *t^i*(1-t)^(N-i)
                  N
    0 <= t <= 1

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 5.5
------------------------------------------------------
输入   :
    XY      数据对，nx2，x-y
输出   :
    sol     解，(N+1)x2，x(t)-y(t)
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

//BernsteinPoly Bernstein Polynomial
func BernsteinPoly(i, N int) goNum.Matrix {
	cni := goNum.Cnm(N, i)
	sol := goNum.ZeroMatrix(N+1, 1)
	soltemp := goNum.ZeroMatrix(N+1, 1)

	soltemp.Data[0] = 1.0
	soltemp.Data[1] = -1.0 //1-t
	//(1-t)^(N-i)
	if N-i > 1 {
		for j := 2; j < N-i+1; j++ {
			for k := j; k > 0; k-- {
				soltemp.Data[k] = soltemp.Data[k] - soltemp.Data[k-1]
			}
		}
	}
	//(1-t)^(N-i) * t^i
	for j := N; j >= i; j-- {
		sol.Data[j] = float64(cni) * soltemp.Data[j-i]
	}

	return sol
}

// FittingBezier Bezier曲线拟合控制点
func FittingBezier(XY goNum.Matrix) (goNum.Matrix, bool) {
	/*
		Bezier曲线拟合控制点
		输入   :
		    XY      数据对，nx2，x-y
		输出   :
		    sol     解，(N+1)x2，x(t)-y(t)
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断维数
	if XY.Columns < 2 {
		panic("Error in goNum.FittingBezier: At least 2 columns of XY needed")
	}
	n := XY.Rows - 1 //N-1
	sol := goNum.ZeroMatrix(n+1, 2)
	var err bool = false

	//计算
	for i := 0; i < n+1; i++ { //n+1项BernsteinPoly
		soltemp := goNum.BernsteinPoly(i, n)
		xi := XY.GetFromMatrix(i, 0)
		yi := XY.GetFromMatrix(i, 1)
		for j := 0; j < n+1; j++ { //n次BernsteinPoly
			sol.SetMatrix(j, 0, sol.GetFromMatrix(j, 0)+xi*soltemp.Data[j])
			sol.SetMatrix(j, 1, sol.GetFromMatrix(j, 1)+yi*soltemp.Data[j])
		}
	}

	err = true
	return sol, err
}

func BenchmarkFittingBezier(b *testing.B) {
	xy49 := goNum.NewMatrix(4, 2, []float64{
		2.0, 2.0,
		1.0, 1.5,
		3.5, 0.0,
		4.0, 1.0})
	for i := 0; i < b.N; i++ {
		goNum.FittingBezier(xy49)
	}
}
