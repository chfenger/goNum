// Max
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-19
版本   : 0.0.0
------------------------------------------------------
    向量第一个最大值及其位置
------------------------------------------------------
输入   :
    a       a 被处理向量
输出   :
    sol     解值
    ii      第一个最大值位置
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

func Max(a []float64) (float64, int, bool) {
	/*
		向量第一个最大值及其位置
		输入   :
		    a       a 被处理向量
		输出   :
		    sol     解值
		    ii      第一个最大值位置
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var sol float64
	var ii int
	var err bool = false

	n := len(a)
	ii = 0
	sol = a[ii]
	for i := 1; i < n; i++ {
		if sol < a[i] {
			ii = i
			sol = a[i]
		}
	}

	err = true
	return sol, ii, err
}
