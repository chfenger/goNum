// SelectSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-05
版本   : 0.0.0
------------------------------------------------------
    选择排序法
理论：
    时间复杂度: O(n^2)
    最好情况  : O(n^2)
    最坏情况  : O(n^2)
    空间复杂度: O(1)
    稳定性    : 不稳定
------------------------------------------------------
输入   :
    in      输入矩阵, 1xn
输出   :
    sol     排序结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// SelectSort 选择排序法
func SelectSort(in goNum.Matrix) (goNum.Matrix, bool) {
	/*
	      选择排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if in.Rows != 1 {
		panic("Error in goNum.SelectSort: Input Matrix error")
	}
	if in.Columns < 1 {
		panic("Error in goNum.SelectSort: Empty input Matrix")
	} else if in.Columns == 1 {
		return in, true
	}

	n := in.Columns
	sol := goNum.ZeroMatrix(1, n)
	var err bool = false

	//初始化sol
	for i := 0; i < n; i++ {
		sol.Data[i] = in.Data[i]
	}
	//排序开始
	for i := 0; i < n-1; i++ {
		mini := i
		for j := i + 1; j < n; j++ {
			if sol.Data[mini] > sol.Data[j] {
				mini = j
			}
		}
		sol.Data[i], sol.Data[mini] = sol.Data[mini], sol.Data[i]
	}

	err = true
	return sol, err
}

func BenchmarkSelectSort(b *testing.B) {
	x59 := goNum.NewMatrix(1, 5, []float64{5.0, 3.2, 1.8, 2.4, 0.1})
	for i := 0; i < b.N; i++ {
		goNum.SelectSort(x59)
	}
}
