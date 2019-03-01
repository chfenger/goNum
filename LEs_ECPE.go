// LEs_ECPE
// linear equations - elemination of column principle
//    element
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-19
版本   : 0.0.0
------------------------------------------------------
    线性代数方程组的列主元消去法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 47-49.

    乘除运算的次数 n^3/3+n^2-n/3
------------------------------------------------------
输入   :
    a       a x = b线性代数方程组的系数矩阵
    b       a x = b线性代数方程组的右侧常数列向量
输出   :
    sol     解值
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// LEs_ECPE 线性代数方程组的列主元消去法
func LEs_ECPE(a [][]float64, b []float64) ([]float64, bool) {
	/*
		线性代数方程组的列主元消去法
		输入   :
		    a       a x = b线性代数方程组的系数矩阵
		    b       a x = b线性代数方程组的右侧常数列向量
		输出   :
		    sol     解值
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//方程个数为n
	var err bool = false
	atemp := a
	btemp := b
	n := len(btemp)
	sol := make([]float64, n)
	temp0 := make([]float64, n)
	var temp1 float64

	// 输入判断
	if len(atemp) != n {
		return sol, err
	}

	//求解
	//消去，求得上三角矩阵
	for i := 0; i < n-1; i++ {
		//求第i列的主元素并调整顺序
		acol := make([]float64, n-i)
		for icol := i; icol < n; icol++ {
			acol[icol-i] = atemp[icol][i]
		}
		_, ii, _ := MaxAbs(acol)
		if ii+i != i {
			temp0 = atemp[ii+i]
			atemp[ii+i] = atemp[i]
			atemp[i] = temp0
			temp1 = btemp[ii+i]
			btemp[ii+i] = btemp[i]
			btemp[i] = temp1
		}

		//列消去
		for j := i + 1; j < n; j++ {
			mul := atemp[j][i] / atemp[i][i]
			for k := i; k < n; k++ {
				atemp[j][k] = atemp[j][k] - atemp[i][k]*mul
			}
			btemp[j] = btemp[j] - btemp[i]*mul
		}
	}

	//回代
	sol[n-1] = btemp[n-1] / atemp[n-1][n-1]
	for i := n - 2; i >= 0; i-- {
		temp1 = 0.0
		for j := i + 1; j < n; j++ {
			temp1 = temp1 + atemp[i][j]*sol[j]
		}
		sol[i] = (btemp[i] - temp1) / atemp[i][i]
	}

	//返回结果
	err = true
	return sol, err
}
