// BucketSort_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    桶排序法
理论：
    时间复杂度: O(n+k)
    最好情况  : O(n)
    最坏情况  : O(n^2)
    空间复杂度: O(n+k)
    稳定性    : 稳定
------------------------------------------------------
输入   :
    in      输入矩阵, 1xn
    bucketSize 桶中元素数
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
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

// IntMin 整数切片中最小数
func IntMin(in []int) int {
	min := in[0]
	for i := 1; i < len(in); i++ {
		if min > in[i] {
			min = in[i]
		}
	}
	return min
}

// bucketSort_sort
func bucketSort_sort(temp0 []int, bucketSize int) []int {
	if (temp0 == nil) || (len(temp0) < 2) {
		return temp0
	}
	temp2 := make([]int, 0)
	min := IntMin(temp0)
	max := IntMax(temp0)
	bucketCount := int(math.Floor(float64((max-min)/bucketSize))) + 1
	bucket := make([][]int, bucketCount) //第一维为桶数量，第二维为桶容量|| (bucketSize == 0)

	//排序开始
	//利用映射函数将数据分配到各个桶中
	for i := 0; i < len(temp0); i++ {
		indi := int(math.Floor(float64((temp0[i] - min) / bucketSize)))
		bucket[indi] = append(bucket[indi], temp0[i])
	}
	//桶中排序
	for i := 0; i < bucketCount; i++ {
		if bucketCount == 1 {
			bucketSize--
		}
		temp1 := bucketSort_sort(bucket[i], bucketSize)
		for j := 0; j < len(temp1); j++ {
			temp2 = append(temp2, temp1[j])
		}
	}
	return temp2
}

// BucketSort 桶排序法
func BucketSort(in []int, bucketSize int) ([]int, bool) {
	/*
	      桶排序法
	   输入   :
	       in      输入矩阵, 1xn
	       bucketSize 桶中元素数
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
	soltemp := make([]int, n)
	for i := 0; i < n; i++ {
		soltemp[i] = in[i]
	}

	//排序开始
	sol := bucketSort_sort(soltemp, bucketSize)

	err = true
	return sol, err
}

func BenchmarkBucketSort(b *testing.B) {
	x66 := []int{20, 50, 12, 80, 5, 36}
	for i := 0; i < b.N; i++ {
		goNum.BucketSort(x66，2)
	}
}
