// MatrixEigenJacobiPass_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-30
版本   : 0.0.0
------------------------------------------------------
    求解n阶对称矩阵A的全部特征值及其特征向量，雅可比过关法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 90.
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

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// MatrixEigenJacobiPass 求解n阶对称矩阵A的全部特征值及其特征向量，雅可比过关法
func MatrixEigenJacobiPass(A goNum.Matrix, tol float64, n int) (goNum.Matrix, goNum.Matrix, bool) {
	/*
		求解n阶对称矩阵A的全部特征值及其特征向量，雅可比过关法
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
		return goNum.ZeroMatrix(A.Rows, A.Columns), goNum.ZeroMatrix(A.Rows, A.Columns), false
	}
	//1.
	//Rbar最终为特征向量矩阵
	Rbar := goNum.IdentityE(A.Rows)
	//复制A矩阵为B，B为迭代过程中逐渐改变的矩阵，最终将成为特征值矩阵
	B := goNum.ZeroMatrix(A.Rows, A.Columns)
	Bbar := goNum.ZeroMatrix(A.Rows, A.Columns)
	for i := 0; i < len(A.Data); i++ {
		B.Data[i] = A.Data[i]
	}
	//2. 计算非对角元素平方和
	v0 := sum2Else_MatrixEigenClassicalJacobi(B)
	//3. 设置阀值v1
	v1 := v0 / float64(B.Rows)

	//迭代步
	for i := 0; i < n; i++ {
		for i0 := 0; i0 < B.Rows-1; i0++ {
			for j := i + 1; j < B.Columns; j++ {
				//逐个扫描，判断是否大于阀值
				if B.GetFromMatrix(i, j) > v1 {
					//Jocobi正交相似变换（古典Jocobi法）
					//计算cos(theta)及sin(theta)
					cost, sint := cosSinTheta_MatrixEigenClassicalJacobi(B, i, j)
					//R为迭代矩阵
					R := goNum.IdentityE(A.Rows)
					R.SetMatrix(i, i, cost)
					R.SetMatrix(i, j, sint)
					R.SetMatrix(j, i, -1.0*sint)
					R.SetMatrix(j, j, cost)
					Bbar = goNum.DotPruduct(goNum.DotPruduct(R, B), R.Transpose()) //A1 = RARt
					//Rbar = Rbar*Rt
					Rbar = goNum.DotPruduct(Rbar, R.Transpose())
				}
			}
		}
		//计算并判断v1是否满足误差需求，否则迭代
		v0 = sum2Else_MatrixEigenClassicalJacobi(Bbar)
		if v0 < tol {
			return Bbar, Rbar, true
		}
		v1 = v0 / float64(B.Rows)
		//A = A1
		for i := 0; i < len(Bbar.Data); i++ {
			B.Data[i] = Bbar.Data[i]
		}
	}

	return goNum.ZeroMatrix(A.Rows, A.Columns), goNum.ZeroMatrix(A.Rows, A.Columns), false
}

func BenchmarkMatrixEigenJacobiPass(b *testing.B) {
	A21 := goNum.NewMatrix(3, 3, []float64{2.0, -1.0, 0.0,
		-1.0, 2.0, -1.0,
		0.0, -1.0, 2.0})
	for i := 0; i < b.N; i++ {
		MatrixEigenJacobiPass(A21, 1e-6, 1e6)
	}
}
