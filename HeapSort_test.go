// HeapSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    堆排序法
理论：
    时间复杂度: O(nlog2(n))
    最好情况  : O(nlog2(n))
    最坏情况  : O(nlog2(n))
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

// heapSort_maxHeap 建立最大顶堆
func heapSort_maxHeap(sol *goNum.Matrix, n *int) {
	for i := *n / 2; i >= 0; i-- {
		heapSort_heapify(sol, i, n)
	}
}

// heapSort_heapify 堆调整
func heapSort_heapify(sol *goNum.Matrix, i int, n *int) {
	i0 := 2*i + 1
	i2 := 2*i + 2
	max := i

	if i0 < *n && (*sol).Data[i0] > (*sol).Data[i] {
		max = i0
	}
	if i2 < *n && (*sol).Data[i2] > (*sol).Data[i] {
		max = i2
	}

	if max != i {
		(*sol).Data[i], (*sol).Data[max] = (*sol).Data[max], (*sol).Data[i]
		heapSort_heapify(sol, max, n)
	}
}

// HeapSort 堆排序法
func HeapSort(in goNum.Matrix) (goNum.Matrix, bool) {
	/*
	      堆排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if in.Rows != 1 {
		panic("Error in goNum.HeapSort: Input Matrix error")
	}
	if in.Columns < 1 {
		panic("Error in goNum.HeapSort: Empty input Matrix")
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
	heapSort_maxHeap(&sol, &n)
	for i := in.Columns - 1; i > 0; i-- {
		sol.Data[0], sol.Data[i] = sol.Data[i], sol.Data[0]
		n--
		heapSort_heapify(&sol, 0, &n)
	}

	err = true
	return sol, err
}

func BenchmarkHeapSort(b *testing.B) {
	x64 := goNum.NewMatrix(1, 5, []float64{5.0, 3.2, 1.8, 2.4, 0.1})
	for i := 0; i < b.N; i++ {
		goNum.HeapSort(x64)
	}
}
