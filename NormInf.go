// NormInf
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-21
版本   : 0.0.0
------------------------------------------------------
    求矩阵无穷范数
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 65.
    ||A||Inf = Maxi(Sumj(|aij|))
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

import "math"

func NormInf(A Matrix) (float64, bool) {
	/*
		求矩阵无穷范数
		输入   :
		    A       矩阵
		输出   :
		    sol     范数值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	var sol float64
	var err bool = false

	//求取行绝对值的和
	row := make([]float64, A.Rows)
	for i := 0; i < A.Rows; i++ {
		thisRow := A.RowOfMatrix(i)
		for j := 0; j < len(thisRow); j++ {
			row[i] += math.Abs(thisRow[j])
		}
	}

	//求取最大值
	sol, _, _ = Max(row)

	err = true
	return sol, err
}
