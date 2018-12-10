// LEs_JocobiIterate
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-22
版本   : 0.0.0
------------------------------------------------------
    解n阶线性方程组的Jocobi迭代法（简单迭代法）
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 61-68.
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

package goNum

import (
	"math"
)

func LEs_JocobiIterate(A, b, x0 Matrix, tol float64, n int) ([]float64, bool) {
	/*
		解n阶线性方程组的Jocobi迭代法（简单迭代法）
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
	B := ZeroMatrix(A.Rows, A.Columns)
	g := ZeroMatrix(A.Rows, 1)
	x1 := ZeroMatrix(A.Rows, 1)
	sol := ZeroMatrix(A.Rows, 1)
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
	temp0, _ := Norm1(B)
	temp1, _ := NormInf(B)
	if (temp0 >= 1) || (temp1 >= 1) {
		return sol.Data, err
	}

	//求解
	for i := 0; i < n; i++ {
		x1 = AddMatrix(DotPruduct(B, x0), g)
		sol = SubMatrix(x1, x0)
		max, _, _ := Max(sol.Data)
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
