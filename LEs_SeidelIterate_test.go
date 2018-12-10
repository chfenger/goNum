// LEs_SeidelIterate_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-22
版本   : 0.0.0
------------------------------------------------------
    解n阶线性方程组的Seidel迭代法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 68-72.
    收敛的条件：（B为变化后的系数矩阵）
       1. 矩阵B的谱半径小于1，或者
       2. 矩阵B的1范数小于1，或者
       3. 矩阵B的无穷范数小于1，或者
       4. 系数矩阵A严格对角占优
------------------------------------------------------
输入   :
    A       系数矩阵
    b       常数值向量
    tol     最大容许误差
    n       最大迭代步数
输出   :
    sol     解向量
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

func LEs_SeidelIterate(A, b, x0 goNum.Matrix, tol float64, n int) ([]float64, bool) {
	/*
		解n阶线性方程组的Seidel迭代法
		输入   :
		    A       系数矩阵
		    b       常数值向量
		    tol     最大容许误差
		    n       最大迭代步数
		输出   :
		    sol     解向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	B := goNum.ZeroMatrix(A.Rows, A.Columns)
	g := goNum.ZeroMatrix(A.Rows, 1)
	x1 := goNum.ZeroMatrix(A.Rows, 1)
	xtemp := goNum.ZeroMatrix(A.Rows, 1)
	sol := goNum.ZeroMatrix(A.Rows, 1)
	var err bool = false

	//方程组迭代化变换，求得矩阵B
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Columns; j++ {
			if j != i {
				B.SetMatrix(i, j, -1.0*A.GetFromMatrix(i, j)/A.GetFromMatrix(i, i))
			}
		}
		g.Data[i] = b.Data[i] / A.GetFromMatrix(i, i)
	}

	//判断B，是否收敛
	temp0, _ := goNum.Norm1(B)
	temp1, _ := goNum.NormInf(B)
	if (temp0 >= 1) || (temp1 >= 1) {
		return sol.Data, err
	}

	//求解
	for i := 0; i < n; i++ {
		for i0 := 0; i0 < B.Rows; i0++ {
			dotP := goNum.DotPruduct(goNum.NewMatrix(1, B.Columns, B.RowOfMatrix(i0)), x0)
			x1.Data[i0] = dotP.Data[0] + g.Data[i0]
			xtemp.Data[i0] = x1.Data[i0]
		}
		sol = goNum.SubMatrix(x1, x0)
		max, _, _ := goNum.Max(sol.Data)
		if math.Abs(max) < tol {
			sol = x1
			err = true
			return sol.Data, err
		}

		for i0 := 0; i0 < x0.Rows; i0++ {
			x0.Data[i0] = x1.Data[i0]
		}
	}

	return make([]float64, A.Rows), err
}

func BenchmarkLEs_SeidelIterate(b *testing.B) {
	A17 := goNum.NewMatrix(3, 3, []float64{10.0, -2.0, -1.0, -2.0, 10.0, -1.0, -1.0, -2.0, 5.0})
	b17 := goNum.NewMatrix(3, 1, []float64{3.0, 15.0, 10})
	x17 := goNum.ZeroMatrix(3, 1)

	for i := 0; i < b.N; i++ {
		LEs_SeidelIterate(A17, b17, x17, 1e-6, 1e6)
	}
}
