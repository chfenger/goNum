// MergeSort
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-03-06
版本   : 0.0.0
------------------------------------------------------
    归并排序法
理论：
    时间复杂度: O(nlog2(n))
    最好情况  : O(nlog2(n))
    最坏情况  : O(nlog2(n))
    空间复杂度: O(n)
    稳定性    : 稳定
------------------------------------------------------
输入   :
    in      输入矩阵, 1xn
输出   :
    sol     排序结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// mergeSort_merge
// i0 --- first
// i1 --- mid
// i2 --- last
func mergeSort_merge(sol *Matrix, i0, i1, i2 int) {
	temp := ZeroMatrix(1, (*sol).Columns)
	i := i0
	j := i1 + 1
	k := 0
	for i <= i1 && j <= i2 {
		if (*sol).Data[i] < (*sol).Data[j] {
			temp.Data[k] = (*sol).Data[i]
			k++
			i++
		} else {
			temp.Data[k] = (*sol).Data[j]
			k++
			j++
		}
	}

	for i <= i1 {
		temp.Data[k] = (*sol).Data[i]
		k++
		i++
	}

	for j <= i2 {
		temp.Data[k] = (*sol).Data[j]
		k++
		j++
	}

	for i = 0; i < k; i++ {
		(*sol).Data[i0+i] = temp.Data[i]
	}
}

// mergeSort_sort
// i0 --- first
// i2 --- last
func mergeSort_sort(sol *Matrix, i0, i2 int) {
	if i0 < i2 {
		var i1 int = (i0 + i2) / 2
		mergeSort_sort(sol, i0, i1)
		mergeSort_sort(sol, i1+1, i2)
		mergeSort_merge(sol, i0, i1, i2)
	}
}

// MergeSort 归并排序法
func MergeSort(in Matrix) (Matrix, bool) {
	/*
	      归并排序法
	   输入   :
	       in      输入矩阵, 1xn
	   输出   :
	       sol     排序结果
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断初值维数
	if in.Rows != 1 {
		panic("Error in goNum.MergeSort: Input Matrix error")
	}
	if in.Columns < 1 {
		panic("Error in goNum.MergeSort: Empty input Matrix")
	} else if in.Columns == 1 {
		return in, true
	}

	n := in.Columns
	sol := ZeroMatrix(1, n)
	var err bool = false

	//初始化sol
	for i := 0; i < n; i++ {
		sol.Data[i] = in.Data[i]
	}
	//排序开始
	mergeSort_sort(&sol, 0, n-1)

	err = true
	return sol, err
}
