// FittingTriPoly
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-23
版本   : 0.0.0
------------------------------------------------------
    基于傅立叶（Fourier）级数的三角多项式拟合
理论：
    若f(x)周期为2pi，则存在M（2M<N）阶傅立叶（Fourier）级数
    使得N+1个数据对（xi等距分布）的拟合表示为：
            a0     M
    TM(x) = --- + Sum (aj*cos(jx)+bj*sin(jx))
             2    j=1
    其中
          2  N
    aj = ---Sum yk*cos(j*xk), j=0,1,2,...,M
          N k=1
          2  N
    bj = ---Sum yk*sin(j*xk), j=1,2,...,M
          N k=1

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 5.4.1
------------------------------------------------------
输入   :
    XY      数据对，nx2，x-y
    M       傅立叶级数，< N/2
输出   :
    sol     解，(M+1)x2
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func FittingTriPoly(XY Matrix, M int) (Matrix, bool) {
	/*
		基于傅立叶（Fourier）级数的三角多项式拟合
		输入   :
		    XY      数据对，nx2，x-y
		    M       傅立叶级数，< N/2
		输出   :
		    sol     解，(M+1)x2
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断维数
	if XY.Columns < 2 {
		panic("Error in goNum.FittingTriPoly: At least 2 columns of XY needed")
	}
	N := XY.Rows
	//判断M
	if float64(M) >= float64(N)/2.0 {
		panic("Error in goNum.FittingTriPoly: M is wrong")
	}

	sol := ZeroMatrix(M+1, 2) //b0=0.0
	var err bool = false

	//a0
	var a0 float64
	for k := 1; k < N; k++ {
		// a0 += XY.GetFromMatrix(k, 1) * math.Cos(0.0*XY.GetFromMatrix(k, 0))
		a0 += XY.GetFromMatrix(k, 1)
	}
	sol.SetMatrix(0, 0, 2.0*a0/float64(N))

	//aj, bj
	for j := 1; j < M+1; j++ {
		var aj, bj float64
		for k := 1; k < N; k++ {
			aj += XY.GetFromMatrix(k, 1) * math.Cos(float64(j)*XY.GetFromMatrix(k, 0))
			bj += XY.GetFromMatrix(k, 1) * math.Sin(float64(j)*XY.GetFromMatrix(k, 0))
		}
		sol.SetMatrix(j, 0, 2.0*aj/float64(N))
		sol.SetMatrix(j, 1, 2.0*bj/float64(N))
	}

	err = true
	return sol, err
}
