// LEs_Chasing_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-8
版本   : 0.0.0
------------------------------------------------------
    追赶法求解严格对角占优的三对角矩阵
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 59-61.
------------------------------------------------------
输入   :
    A       系数矩阵, nxn
    BA      常数值向量, nx1
输出   :
    sol     解向量, nx1
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// LEs_Chasing 追赶法求解严格对角占优的三对角矩阵
func LEs_Chasing(A, BA goNum.Matrix) (goNum.Matrix, bool) {
	/*
		追赶法求解严格对角占优的三对角矩阵
		输入   :
		    A       系数矩阵, nxn
		    BA      常数值向量, nx1
		输出   :
		    sol     解向量, nx1
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断A是否方阵
	if A.Rows != A.Columns {
		panic("Error in goNum.LEs_Chasing: A is not a square matrix")
	}
	//判断BA是否与A行数相等
	if A.Rows != BA.Rows {
		panic("Error in goNum.LEs_Chasing: Rows of A and BA are not equal")
	}

	var err bool = false
	n := A.Rows
	ai := goNum.ZeroMatrix(n, 1) //第一位无效
	bi := goNum.ZeroMatrix(n, 1)
	ci := goNum.ZeroMatrix(n-1, 1)
	gamma := goNum.ZeroMatrix(n, 1)   //gammai
	beta := goNum.ZeroMatrix(n, 1)    //beta, 第一位无效
	delta := goNum.ZeroMatrix(n-1, 1) //deltai
	y := goNum.ZeroMatrix(n, 1)       //yi
	sol := goNum.ZeroMatrix(n, 1)     //xi

	//ai, bi, ci
	bi.Data[0] = A.GetFromMatrix(0, 0)
	ci.Data[0] = A.GetFromMatrix(0, 1)
	for i := 1; i < n-1; i++ {
		ai.Data[i] = A.GetFromMatrix(i, i-1)
		bi.Data[i] = A.GetFromMatrix(i, i)
		ci.Data[i] = A.GetFromMatrix(i, i+1)
	}
	ai.Data[n-1] = A.GetFromMatrix(n-1, n-2)
	bi.Data[n-1] = A.GetFromMatrix(n-1, n-1)

	//解gamma, beta和delta
	gamma.Data[0] = bi.Data[0]
	delta.Data[0] = ci.Data[0] / gamma.Data[0]
	for i := 1; i < n-1; i++ {
		beta.Data[i] = ai.Data[i]
		gamma.Data[i] = bi.Data[i] - beta.Data[i]*delta.Data[i-1]
		delta.Data[i] = ci.Data[i] / gamma.Data[i]
	}
	beta.Data[n-1] = ai.Data[n-1]
	gamma.Data[n-1] = bi.Data[n-1] - beta.Data[n-1]*delta.Data[n-2]

	//解yi
	y.Data[0] = BA.Data[0] / gamma.Data[0]
	for i := 1; i < BA.Rows; i++ {
		y.Data[i] = (BA.Data[i] - beta.Data[i]*y.Data[i-1]) / gamma.Data[i]
	}

	//解xi
	sol.Data[n-1] = y.Data[n-1]
	for i := n - 2; i >= 0; i-- {
		sol.Data[i] = y.Data[i] - delta.Data[i]*sol.Data[i+1]
	}

	err = true
	return sol, err
}

func BenchmarkLEs_Chasing(b *testing.B) {
	A29 := goNum.NewMatrix(3, 3, []float64{4.0, -1.0, 0.0,
		-1.0, 4.0, -1.0,
		0.0, -1.0, 4.0})
	BA29 := goNum.NewMatrix(3, 1, []float64{1.0, 4.0, -3.0})
	for i := 0; i < b.N; i++ {
		LEs_Chasing(A29, BA29) //{0.5, 1, -0.5}
	}
}
