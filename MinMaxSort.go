// MinMaxSort
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-19
版本   : 0.0.0
------------------------------------------------------
    向量从小到大的排序
------------------------------------------------------
输入   :
    a       a 被排序向量
输出   :
    sol     解值，指针
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

func MinMaxSort(a []float64) (*[]float64, bool) {
	/*
		向量从小到大的排序
		输入   :
		    a       a 被排序向量
		输出   :
		    sol     解值，指针
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var err bool = false
	var temp float64
	var n int = len(a)
	sol := a

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if sol[i] > sol[j] {
				temp = sol[j]
				sol[j] = sol[i]
				sol[i] = temp
			}
		}
	}
	err = true
	return &sol, err
}
