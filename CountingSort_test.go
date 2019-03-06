// CountingSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    计数排序法
理论：
    时间复杂度: O(n+k)
    最好情况  : O(n+k)
    最坏情况  : O(n+k)
    空间复杂度: O(n+k)
    稳定性    : 稳定
------------------------------------------------------
输入   :
    in      输入切片, 1xn
输出   :
    sol     排序结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
注意：
   仅对整数排序有效
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// IntMax 整数切片中最大数
func IntMax(in []int) int {
	max := in[0]
	for i := 1; i < len(in); i++ {
		if in[i] > max {
			max = in[i]
		}
	}
	return max
}

// CountingSort 计数排序法
func CountingSort(in []int) ([]int, bool) {
	/*
	      计数排序法
	   输入   :
	       in      输入切片, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if len(in) < 1 {
		panic("Error in goNum.CountingSort: Empty input Matrix")
	} else if len(in) == 1 {
		return in, true
	}

	n := len(in)
	sol := make([]int, n)
	max := IntMax(in)
	temp := make([]int, max+1)
	ind := 0
	var err bool = false

	//初始化sol
	for i := 0; i < n; i++ {
		sol[i] = in[i]
	}
	//排序开始

	for i := 0; i < n; i++ {
		// if !temp[sol[i]] {
		// 	temp[sol[i]] = 0
		// }
		temp[sol[i]]++
	}
	for i := 0; i < max+1; i++ {
		for temp[i] > 0 {
			sol[ind] = i
			ind++
			temp[i]--
		}
	}

	err = true
	return sol, err
}

func BenchmarkCountingSort(b *testing.B) {
	x65 := []int{20, 50, 12, 80, 5, 36}
	for i := 0; i < b.N; i++ {
		goNum.CountingSort(x65)
	}
}
