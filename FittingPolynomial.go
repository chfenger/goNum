// FittingPolynomial
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-11
版本   : 0.0.0
------------------------------------------------------
    多项式拟合
理论：
    对于单自变量单因变量的N个数据对
    假设其一个低于N-1次的多项式为：
    y(x) = a0 + a1x + a2x^2 + ... + amx^m (m < N-1)
    建立矛盾方程组Ax=b，即

     N
    Sum ai*xj^i = bj   (i=0, 1, 2, ..., m)
    j=1

    求解ai (i=0, 1, 2, ..., m)代入多项式即得拟合函数

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 136-138.
------------------------------------------------------
输入   :
    xy      单自变量单因变量的N个数据对，Nx2
    m       多项式次数，m < N-1
输出   :
    sol     解向量，从0到m对应a0到am
    RMS     均方误差
    MaxErr  最大误差
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
扩展   ：
    可以修改为适应log、exp、sin等拟合方法
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// FittingPolynomial 多项式拟合
func FittingPolynomial(xy Matrix, m int) (Matrix, float64, float64, bool) {
	/*
		多项式拟合
		输入   :
		    xy      单自变量单因变量的N个数据对，Nx2
		    m       多项式次数，m < N-1
		输出   :
		    sol     解向量，从0到m对应a0到am
		    RMS     均方误差
		    MaxErr  最大误差
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断m是否小于N-1
	if (m > xy.Rows-2) || (m < 0) {
		panic("Error in goNum.FittingPolynomial: Order m is wrong number")
	}

	N := xy.Rows

	//构建矛盾方程组系数矩阵A, b=xy.ColumnOfMatrix(1)
	A := ZeroMatrix(N, m+1)
	for i := 0; i < N; i++ {
		A.SetMatrix(i, 0, 1.0)
		for j := 1; j < m+1; j++ {
			temp := xy.GetFromMatrix(i, 0)
			A.SetMatrix(i, j, math.Pow(temp, float64(j)))
		}
	}
	//求解矛盾方程组
	sol, err := InconsistentLSQ(A, Slices1ToMatrix(xy.ColumnOfMatrix(1)))
	//判断结果
	if err != true {
		panic("Error in goNum.FittingPolynomial: Solve error")
	}

	errSub := make([]float64, N)
	var RMS float64
	for i := 0; i < N; i++ {
		fit := sol.Data[0]
		for j := 1; j < m+1; j++ {
			fit += sol.Data[j] * math.Pow(xy.GetFromMatrix(i, 0), float64(j))
		}
		errSub[i] = fit - xy.GetFromMatrix(i, 1)
		RMS += errSub[i] * errSub[i]
	}
	RMS = math.Sqrt(RMS)
	MaxErr, _, _ := MaxAbs(errSub)

	return sol, RMS, math.Abs(MaxErr), err
}
