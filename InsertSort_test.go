// InsertSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-05
版本   : 0.0.0
------------------------------------------------------
    插入排序法
理论：
    时间复杂度: O(n^2)
    最好情况  : O(n)
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

// InsertSort 插入排序法
func InsertSort(in goNum.Matrix) (goNum.Matrix, bool) {
	/*
	      插入排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if in.Rows != 1 {
		panic("Error in goNum.InsertSort: Input Matrix error")
	}
	if in.Columns < 1 {
		panic("Error in goNum.InsertSort: Empty input Matrix")
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
	for i := 1; i < n; i++ {
		index := i - 1
		min := sol.Data[i]
		for (index >= 0) && (sol.Data[index] > min) {
			sol.Data[index+1] = sol.Data[index]
			index--
		}
		sol.Data[index+1] = min
	}

	err = true
	return sol, err
}

func BenchmarkInsertSort(b *testing.B) {
	x60 := goNum.NewMatrix(1, 5, []float64{5.0, 3.2, 1.8, 2.4, 0.1})
	for i := 0; i < b.N; i++ {
		goNum.InsertSort(x60)
	}
}
