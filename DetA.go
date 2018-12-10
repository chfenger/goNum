// DetA
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-20
版本   : 0.0.0
------------------------------------------------------
    求矩阵行列式值的列主元消去法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 52.
------------------------------------------------------
输入   :
    a       矩阵
输出   :
    sol     解值，值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func DetA(a [][]float64) (float64, bool) {
	/*
		求矩阵行列式值的列主元消去法
		输入   :
		    a       矩阵
		输出   :
		    sol     解值，值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	var sol float64 = 1.0
	var err bool = false
	var count0 int
	n := len(a)
	temp0 := make([]float64, n)

	// 判断是否方阵
	if len(a) != len(a[0]) {
		return sol, err
	}

	//主元消去
	for i := 0; i < n; i++ {
		//求第i列的主元素并调整行顺序
		acol := make([]float64, n-i)
		for icol := i; icol < n; icol++ {
			acol[icol-i] = a[icol][i]
		}
		_, ii, _ := MaxAbs(acol)
		if ii+i != i {
			count0++
			temp0 = a[ii+i]
			a[ii+i] = a[i]
			a[i] = temp0
		}

		//列消去
		for j := i + 1; j < n; j++ {
			mul := a[j][i] / a[i][i]
			for k := i; k < n; k++ {
				a[j][k] = a[j][k] - a[i][k]*mul
			}
		}
		sol = math.Pow(-1.0, float64(count0)) * sol * a[i][i]
	}

	err = true
	return sol, err
}
