// LLT_Decompose
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-8
版本   : 0.0.0
------------------------------------------------------
    求对称正定矩阵的平方根分解法
理论：
    A = LL'
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 57-58.
------------------------------------------------------
输入   :
    A       矩阵,对称正定
输出   :
    L       下三角矩阵, 上三角矩阵为其转置
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func LLT_Decompose(A Matrix) (Matrix, bool) {
	/*
		求对称正定矩阵的平方根分解法
		输入   :
		    A       矩阵,对称正定
		输出   :
		    L
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断对称
	if A.Rows != A.Columns {
		panic("Error in goNum.LLT_Decompose: A is not symmetry")
	}
	n := A.Rows
	L := ZeroMatrix(n, n)
	var err bool = false

	//计算开始
	//第一列
	L.SetMatrix(0, 0, math.Sqrt(A.GetFromMatrix(0, 0)))
	l11 := L.GetFromMatrix(0, 0)
	for j := 1; j < n; j++ {
		L.SetMatrix(j, 0, A.GetFromMatrix(0, j)/l11)
	}
	//其它列
	for k := 1; k < n; k++ {
		//主对角元lkk
		var temp0 float64
		for m := 0; m < k; m++ {
			temp0 += L.GetFromMatrix(k, m) * L.GetFromMatrix(k, m)
		}
		temp0 = A.GetFromMatrix(k, k) - temp0
		L.SetMatrix(k, k, math.Sqrt(temp0))
		//k列其它元
		for j := k + 1; j < n; j++ {
			var temp1 float64
			for m := 0; m < k; m++ {
				temp1 += L.GetFromMatrix(k, m) * L.GetFromMatrix(j, m)
			}
			temp1 = (A.GetFromMatrix(k, j) - temp1) / L.GetFromMatrix(k, k)
			L.SetMatrix(j, k, temp1)
		}
	}

	err = true
	return L, err
}
