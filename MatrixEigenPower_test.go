// MatrixEigenPower_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-23
版本   : 0.0.0
------------------------------------------------------
    求解n阶矩阵A的主特征值（按模最大）及其特征向量
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 78-81.
------------------------------------------------------
输入   :
    A       系数矩阵
    u       n维初始向量
    tol     最大容许误差
    n       最大迭代步数
输出   :
    sol     主特征值
    v       主特征值所对应的特征向量
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

// MatrixEigenPower 求解n阶矩阵A的主特征值（按模最大）及其特征向量
func MatrixEigenPower(A, u0 goNum.Matrix, tol float64, n int) (float64, []float64, bool) {
	/*
		求解n阶矩阵A的主特征值（按模最大）及其特征向量
		输入   :
		    A       系数矩阵
		    u       n维初始向量
		    tol     最大容许误差
		    n       最大迭代步数
		输出   :
		    sol     主特征值
		    v       主特征值所对应的特征向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断输入正确与否
	if A.Rows != u0.Rows {
		panic("goNum.MatrixEigenPower: A and u are not matched")
	}

	u1 := goNum.ZeroMatrix(u0.Rows, u0.Columns)
	var l0, l1 float64
	v1 := make([]float64, u0.Rows)
	var err bool = false
	var j int

	u1 = goNum.DotPruduct(A, u0)
	for i0 := 0; i0 < u0.Rows; i0++ {
		if (math.Abs(u0.Data[i0]) > 1e-3) && (math.Abs(u1.Data[i0]) > 1e-3) {
			j = i0
			l0 = u1.Data[i0] / u0.Data[i0]
		}
		u0.Data[i0] = u1.Data[i0]
	}

	for i := 0; i < n; i++ {
		u1 = goNum.DotPruduct(A, u0)
		l1 = u1.Data[j] / u0.Data[j]

		//计算最大值，并进行规范化处理
		for i0 := 0; i0 < u0.Rows; i0++ {
			v1[i0] = math.Abs(u1.Data[i0])
		}
		_, j0, _ := goNum.Max(v1)
		max := u1.Data[j0]
		if max > 1e6 {
			for i0 := 0; i0 < u0.Rows; i0++ {
				u1.Data[i0] = u1.Data[i0] / max
			}
		}

		//判断算出否，并计算对应的特征向量
		if math.Abs(l1-l0) < tol {
			for i0 := 0; i0 < u0.Rows; i0++ {
				u1.Data[i0] = u1.Data[i0] / max
			}
			err = true
			return l1, u1.Data, err
		}

		//准备下次迭代
		l0 = l1
		for i0 := 0; i0 < u0.Rows; i0++ {
			u0.Data[i0] = u1.Data[i0]
		}
	}

	return 0.0, make([]float64, u0.Rows), err
}

func BenchmarkMatrixEigenPower(b *testing.B) {
	A19 := goNum.NewMatrix(3, 3, []float64{6.0, -12.0, 6.0,
		-21.0, -3.0, 24.0,
		-12.0, -12.0, 51.0})
	u19 := goNum.NewMatrix(3, 1, []float64{1.0, 1.0, 1.0})

	for i := 0; i < b.N; i++ {
		MatrixEigenPower(A19, u19, 1e-3, 1e3)
	}
}
