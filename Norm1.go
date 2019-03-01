// Norm1
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-21
版本   : 0.0.0
------------------------------------------------------
    求矩阵1范数
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 65.
    ||A||1 = Maxj(Sumi(|aij|))
------------------------------------------------------
输入   :
    A       矩阵
输出   :
    sol     范数值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// Norm1 求矩阵1范数
func Norm1(A Matrix) (float64, bool) {
	/*
		求矩阵1范数
		输入   :
		    A       矩阵
		输出   :
		    sol     范数值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	var sol float64
	var err bool = false

	//求取列绝对值的和
	col := make([]float64, A.Columns)
	for j := 0; j < A.Columns; j++ {
		thisColumn := A.ColumnOfMatrix(j)
		for i := 0; i < len(thisColumn); i++ {
			col[j] += math.Abs(thisColumn[i])
		}
	}

	//求取最大值
	sol, _, _ = Max(col)

	err = true
	return sol, err
}
