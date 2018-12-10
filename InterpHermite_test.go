// InterpHermite_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-7
版本   : 0.0.0
------------------------------------------------------
    计算x点不高于2n+1次Hermite插值结果，拟合n+1个函数值数据
        点和对应的n+1个一阶导数点
    满阶插值，即阶数不高于2n+1
理论：
                n
    H2n+1(x) = Sum (alphaj(x)*yj+betaj(x)*mj)
               j=0
    yj, mj分别为函数值和一阶导数值
                            n         1
    alphaj(x) = (1-2(x-xj)*Sum     -------)lj^2(x)
                          k=0,k!=j  xj-xk
    betaj(x) = (x-xj)lj^2(x)
              (x-x0)(x-x1)...(x-xn)
    lj(x) = --------------------------, (被减数不含xj项)
             (xj-x0)(xj-x1)...(xj-xn)
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 111-113.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；第三列y'i
    xq      插值点, xq!=xi
输出   :
    sol     xq点插值结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

//求解lj(xq)
func ljxq_InterpHermite(A goNum.Matrix, xq float64, j int) float64 {
	sol := 1.0
	xj := A.GetFromMatrix(j, 0)
	for i := 0; i < A.Rows; i++ {
		xi := A.GetFromMatrix(i, 0)
		if i != j {
			sol = sol * (xq - xi) / (xj - xi)
		}
	}
	return sol
}

//求解alplhaj(xq)
func alphajxq_InterpHermite(A goNum.Matrix, xq float64, j int) float64 {
	var temp0 float64
	xj := A.GetFromMatrix(j, 0)
	for k := 0; k < A.Rows; k++ {
		if k != j {
			temp0 += 1.0 / (xj - A.GetFromMatrix(k, 0))
		}
	}
	temp1 := ljxq_InterpHermite(A, xq, j)
	return (1.0 - 2.0*(xq-xj)*temp0) * temp1 * temp1
}

//求解betaj(xq)
func betajxq_InterpHermite(A goNum.Matrix, xq float64, j int) float64 {
	xj := A.GetFromMatrix(j, 0)
	temp0 := ljxq_InterpHermite(A, xq, j)
	return (xq - xj) * temp0 * temp0
}

func InterpHermite(A goNum.Matrix, xq float64) (float64, bool) {
	/*
		计算x点不高于2n+1次Hermite插值结果，拟合n+1个函数值数据点和对应的n+1个一阶导数点
		输入   :
		    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；第三列y'i
		    xq      插值点, xq!=xi
		输出   :
		    sol     xq点插值结果
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断A列数是否为3
	if A.Columns != 3 {
		panic("Error in goNum.InterpHermite: give me xi, yi and y'i")
	}

	var sol float64
	var err bool = false

	//开始计算
	for j := 0; j < A.Rows; j++ {
		sol += alphajxq_InterpHermite(A, xq, j) * A.GetFromMatrix(j, 1)
		sol += betajxq_InterpHermite(A, xq, j) * A.GetFromMatrix(j, 2)
	}

	err = true
	return sol, err
}

func BenchmarkInterpHermite(b *testing.B) {
	//2.4x^3+1.5x^2+0.3x-1.63
	A26 := goNum.NewMatrix(4, 3, []float64{-10.0, -2254.63, 690.3,
		-4.0, -132.43, 103.5,
		4.0, 177.17, 127.5,
		10.0, 2551.37, 750.3})
	for i := 0; i < b.N; i++ {
		InterpHermite(A26, 2.4) //40.9076
	}
}
