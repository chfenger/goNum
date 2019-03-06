// QuickSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    快速排序法
理论：
    时间复杂度: O(nlog2(n))
    最好情况  : O(nlog2(n))
    最坏情况  : O(n^2)
    空间复杂度: O(nlog2(n))
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

// quickSort_sort
// i0 --- first
// i2 --- last
func quickSort_sort(sol *goNum.Matrix, i0, i2 int) {
	if i0 >= i2 {
		return
	}
	i := i0
	j := i2
	ref := (*sol).Data[i] //第一个元素作为分区元素

	for i != j {
		for i < j && (*sol).Data[j] > ref {
			j -= 1
		}
		(*sol).Data[i] = (*sol).Data[j]

		for i < j && (*sol).Data[i] < ref {
			i += 1
		}
		(*sol).Data[j] = (*sol).Data[i]
	}
	(*sol).Data[i] = ref
	quickSort_sort(sol, i0, i-1)
	quickSort_sort(sol, i+1, i2)
}

// QuickSort 快速排序法
func QuickSort(in goNum.Matrix) (goNum.Matrix, bool) {
	/*
	      快速排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if in.Rows != 1 {
		panic("Error in goNum.QuickSort: Input Matrix error")
	}
	if in.Columns < 1 {
		panic("Error in goNum.QuickSort: Empty input Matrix")
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
	quickSort_sort(&sol, 0, n-1)

	err = true
	return sol, err
}

func BenchmarkQuickSort(b *testing.B) {
	x63 := goNum.NewMatrix(1, 5, []float64{5.0, 3.2, 1.8, 2.4, 0.1})
	for i := 0; i < b.N; i++ {
		goNum.QuickSort(x63)
	}
}
