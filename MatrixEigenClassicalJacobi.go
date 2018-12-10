// MatrixEigenClassicalJacobi
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-30
版本   : 0.0.0
------------------------------------------------------
    求解n阶对称矩阵A的全部特征值及其特征向量，经典雅可比法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 84-89.
------------------------------------------------------
输入   :
    A       系数矩阵
    tol     最大容许误差
    n       最大迭代步数
输出   :
    Bbar    主特征值矩阵（n阶对角矩阵）
    Rbar    主特征值所对应的特征向量（n维矩阵，第i列即对应于
            第i个特征值的特征向量）
   (err)    解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

//判断矩阵是否对称
func isSymMatrix_MatrixEigenClassicalJacobi(A Matrix) bool {
	//是否方阵
	if A.Columns != A.Rows {
		return false
	}

	//是否对称
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Columns; j++ {
			if j != i {
				if A.GetFromMatrix(i, j) != A.GetFromMatrix(j, i) {
					return false
				}
			}
		}
	}

	return true
}

//取最大非对角元素
func maxElementElse_MatrixEigenClassicalJacobi(A Matrix) (int, int) {
	var p, q int
	var max float64
	for i := 0; i < A.Rows-1; i++ {
		for j := i + 1; j < A.Rows; j++ {
			c := A.GetFromMatrix(i, j)
			if math.Abs(max) < math.Abs(c) {
				max = c
				p = i
				q = j
			}
		}
	}
	return p, q
}

//计算cos(theta)及sin(theta)
func cosSinTheta_MatrixEigenClassicalJacobi(A Matrix, p, q int) (float64, float64) {
	apq := A.GetFromMatrix(p, q)
	app := A.GetFromMatrix(p, p)
	aqq := A.GetFromMatrix(q, q)

	//情况1，app-aqq=0
	if app == aqq {
		c := math.Sqrt(2.0) / 2.0
		switch {
		case apq < 0:
			return c, -1.0 * c
		case apq > 0:
			return c, c
		default:
			panic("MatrixEigenClassicalJacobi/cosSinTheta: A(p, q) = 0")
		}
	}

	//其他情况
	c := (app - aqq) / (2.0 * apq)
	d := 2.0 * apq / (app - aqq)
	//如果|apq| << |app - aqq|，用d
	if 1000.0*math.Abs(apq) < math.Abs(app-aqq) {
		tant := d / (1.0 + math.Sqrt(1.0+d*d))
		cost := 1.0 / math.Sqrt(1.0+tant*tant)
		sint := tant * cost
		return cost, sint
	}
	// 用c
	switch {
	case c > 0:
		tant := 1.0 / (c + math.Sqrt(1.0+c*c))
		cost := 1.0 / math.Sqrt(1.0+tant*tant)
		sint := tant * cost
		return cost, sint
	case c < 0:
		tant := -1.0 / (math.Abs(c) + math.Sqrt(1.0+c*c))
		cost := 1.0 / math.Sqrt(1.0+tant*tant)
		sint := tant * cost
		return cost, sint
	default:
		panic("MatrixEigenClassicalJacobi/cosSinTheta: A(p, q) = 0, c = 0")
	}
	return 0.0, 0.0
}

//非主对角元素平方之和
func sum2Else_MatrixEigenClassicalJacobi(A Matrix) float64 {
	var sum2 float64
	for i := 0; i < A.Rows-1; i++ {
		for j := i + 1; j < A.Columns; j++ {
			if i != j {
				sum2 += A.GetFromMatrix(i, j) * A.GetFromMatrix(i, j)
			}
		}
	}
	return 2.0 * sum2
}

func MatrixEigenClassicalJacobi(A Matrix, tol float64, n int) (Matrix, Matrix, bool) {
	/*
		求解n阶对称矩阵A的全部特征值及其特征向量，经典雅可比法
		输入   :
		    A       系数矩阵
		    tol     最大容许误差
		    n       最大迭代步数
		输出   :
		    Bbar    主特征值矩阵（n阶对角矩阵）
		    Rbar    主特征值所对应的特征向量（n维矩阵，第i列即对应于
		            第i个特征值的特征向量）
		   (err)    解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	//判断A是否对称矩阵
	if !isSymMatrix_MatrixEigenClassicalJacobi(A) {
		return ZeroMatrix(A.Rows, A.Columns), ZeroMatrix(A.Rows, A.Columns), false
	}

	//第一步
	//Rbar最终为特征向量矩阵
	Rbar := IdentityE(A.Rows)
	//复制A矩阵为B，B为迭代过程中逐渐改变的矩阵，最终将成为特征值矩阵
	B := ZeroMatrix(A.Rows, A.Columns)
	Bbar := ZeroMatrix(A.Rows, A.Columns)
	for i := 0; i < len(A.Data); i++ {
		B.Data[i] = A.Data[i]
	}

	//迭代步
	for i := 0; i < n; i++ {
		//第二步，最大元素所在行列号
		p, q := maxElementElse_MatrixEigenClassicalJacobi(B)
		//第三步，计算cos(theta)及sin(theta)
		cost, sint := cosSinTheta_MatrixEigenClassicalJacobi(B, p, q)
		//第四步，计算R及Rbar，迭代B矩阵
		//R为迭代矩阵
		R := IdentityE(A.Rows)
		R.SetMatrix(p, p, cost)
		R.SetMatrix(p, q, sint)
		R.SetMatrix(q, p, -1.0*sint)
		R.SetMatrix(q, q, cost)
		Bbar = DotPruduct(DotPruduct(R, B), R.Transpose()) //A1 = RARt
		//Rbar = Rbar*Rt
		Rbar = DotPruduct(Rbar, R.Transpose())
		//第五步，判断误差
		if sum2Else_MatrixEigenClassicalJacobi(Bbar) <= tol {
			return Bbar, Rbar, true
		}
		//A = A1
		for i := 0; i < len(Bbar.Data); i++ {
			B.Data[i] = Bbar.Data[i]
		}
	}

	return ZeroMatrix(A.Rows, A.Columns), ZeroMatrix(A.Rows, A.Columns), false
}
