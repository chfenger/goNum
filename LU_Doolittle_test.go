// LU_Doolittle_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-21
版本   : 0.0.0
------------------------------------------------------
    求矩阵Doolittlede LU分解
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 53-56.
------------------------------------------------------
输入   :
    A       矩阵
输出   :
    L, U    下三角矩阵和上三角矩阵
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

func LU_Doolittle(A goNum.Matrix) (goNum.Matrix, goNum.Matrix, bool) {
	/*
		求矩阵Doolittlede LU分解
		输入   :
		    A       矩阵
		输出   :
		    L, U    下三角矩阵和上三角矩阵
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var err bool = false

	if A.Rows != A.Columns {
		panic("goNum.LU_Doolittle: A is not a square matrix")
	}

	L := goNum.ZeroMatrix(A.Rows, A.Columns)
	U := goNum.ZeroMatrix(A.Rows, A.Columns)

	for j := 0; j < A.Rows; j++ {
		U.SetMatrix(0, j, A.GetFromMatrix(0, j))
	}
	for i := 1; i < A.Rows; i++ {
		L.SetMatrix(i, 0, A.GetFromMatrix(i, 0)/U.GetFromMatrix(0, 0))
	}

	for k := 1; k < A.Rows; k++ {

		for j := k; j < A.Rows; j++ {
			var sum float64
			for m := 0; m < k; m++ {
				sum += L.GetFromMatrix(k, m) * U.GetFromMatrix(m, j)
			}
			U.SetMatrix(k, j, A.GetFromMatrix(k, j)-sum)
		}

		for i := k + 1; i < A.Rows; i++ {
			var sum float64
			for m := 0; m < k; m++ {
				sum += L.GetFromMatrix(i, m) * U.GetFromMatrix(m, k)
			}
			L.SetMatrix(i, k, (A.GetFromMatrix(i, k)-sum)/U.GetFromMatrix(k, k))
		}
	}

	//L.PrintMatrix()
	//U.PrintMatrix()

	err = true
	return L, U, err
}

func BenchmarkLU_Doolittle(b *testing.B) {
	A := goNum.ZeroMatrix(3, 3)
	A.SetMatrix(0, 0, 2.0)
	A.SetMatrix(0, 1, 1.0)
	A.SetMatrix(0, 2, 2.0)
	A.SetMatrix(1, 0, 4.0)
	A.SetMatrix(1, 1, 5.0)
	A.SetMatrix(1, 2, 4.0)
	A.SetMatrix(2, 0, 6.0)
	A.SetMatrix(2, 1, -3.0)
	A.SetMatrix(2, 2, 5.0)

	for i := 0; i < b.N; i++ {
		LU_Doolittle(A)
	}
}
