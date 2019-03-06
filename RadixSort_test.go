// RadixSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    基数排序法
理论：
    时间复杂度: O(n*k)
    最好情况  : O(n*k)
    最坏情况  : O(n*k)
    空间复杂度: O(n+k)
    稳定性    : 稳定
------------------------------------------------------
输入   :
    in      输入矩阵, 1xn
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

// RadixSort 基数排序法
func RadixSort(in []int) ([]int, bool) {
	/*
	      基数排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if len(in) < 1 {
		panic("Error in goNum.BucketSort: Empty input Matrix")
	} else if len(in) == 1 {
		return in, true
	}

	n := len(in)
	var err bool = false
	var maxDigit, mod, div int = 0, 10, 1

	soltemp := make([]int, n)
	sol := make([]int, n)
	for i := 0; i < n; i++ {
		soltemp[i] = in[i]
	}
	temp := make([][]int, 10)

	//排序开始
	//最大数的位数
	max := IntMax(soltemp)
	for max != 0 {
		max = max / 10
		maxDigit++
	}

	for i := 0; i < maxDigit; i++ {
		for j := 0; j < n; j++ {
			var num int = (soltemp[j] % mod) / div
			temp[num] = append(temp[num], soltemp[j])
		}
		ind := 0
		for j := 0; j < len(temp); j++ {
			for k := 0; k < len(temp[j]); k++ {
				sol[ind] = temp[j][k]
				ind++
				temp[j] = []int{} //必须置空
			}
		}
		mod = mod * 10
		div = div * 10
	}

	err = true
	return sol, err
}

func BenchmarkRadixSort(b *testing.B) {
	x67 := []int{20, 50, 12, 80, 5, 36}
	for i := 0; i < b.N; i++ {
		goNum.RadixSort(x67)
	}
}
