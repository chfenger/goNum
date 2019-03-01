// InverseA
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-20
版本   : 0.0.0
------------------------------------------------------
    求矩阵逆的列主元消去法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 51.
------------------------------------------------------
输入   :
    a       矩阵
输出   :
    sol     解值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// InverseA 求矩阵逆的列主元消去法
func InverseA(a [][]float64) ([][]float64, bool) {
	/*
		求矩阵逆的列主元消去法
		输入   :
		    a       矩阵
		输出   :
		    sol     解值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	var err bool = false
	n := len(a)
	temp0, _ := E_Mat(n)
	b := temp0
	sol := b
	temp1 := make([]float64, n)

	//判断是否方阵
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
			temp1 = a[ii+i]
			a[ii+i] = a[i]
			a[i] = temp1
			temp1 = b[ii+i]
			b[ii+i] = b[i]
			b[i] = temp1
		}

		//列消去
		//本行主元置一
		mul := a[i][i]
		for j := 0; j < n; j++ {
			a[i][j] = a[i][j] / mul
			b[i][j] = b[i][j] / mul
		}
		//其它列置零
		for j := 0; j < n; j++ {
			if j != i {
				mul = a[j][i] / a[i][i]
				for k := 0; k < n; k++ {
					a[j][k] = a[j][k] - a[i][k]*mul
					b[j][k] = b[j][k] - b[i][k]*mul
				}
			}
		}
	}

	sol = b
	err = true
	return sol, err
}
