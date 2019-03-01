// SearchByStep_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-10-31
版本   : 0.0.0
------------------------------------------------------
    此程序设计使用搜索法来求解连续、单自变量函数指定有限区间
上的解
------------------------------------------------------
输入   :
    fn      函数，定义为等式左侧部分，右侧为零
    a, b    求解区间，一般要求a<b，但不严格
    N       步数，区间细分粒度
    tol     误差上限
输出   :
    sol     解值
    err     解出标志：false-未全部解出；true-全部解出
------------------------------------------------------
*/
package goNum_test

import (
	"math"
	"testing"
)

// SearchByStep 搜索法来求解连续、单自变量函数指定有限区间上的解
func SearchByStep(fn func(float64) float64, a, b float64,
	N int, tol float64) ([]float64, bool) {
	/*
		搜索法来求解连续、单自变量函数指定有限区间上的解
		输入   :
		    fn      函数，定义为等式左侧部分，右侧为零
		    a, b    求解区间，一般要求a<b，但不严格
		    N       步数，区间细分粒度
		    tol     误差上限
		输出   :
		    sol     解值
		    err     解出标志：false-未全部解出；true-全部解出
	*/
	//初始化
	ab0 := make([]float64, 0, 1000)
	ab1 := make([]float64, 0, 1000)
	sol := make([]float64, 0, 1000)
	err := false

	j := 0                    //解的数量
	h := (b - a) / float64(N) //搜索步长，应小于最近两解的距离
	//确定单解区间，并存入对应数组
	for i := 1; i < N+1; i++ {
		if (fn(a+float64(i)*h) > 0 && fn(a+float64(i-1)*h) < 0) || (fn(a+float64(i)*h) < 0 && fn(a+float64(i-1)*h) > 0) {
			ab0 = append(ab0, a+float64(i-1)*h)
			ab1 = append(ab1, a+float64(i)*h)
			sol = append(sol, (ab0[j]+ab1[j])/2.0)
			j++
		}
	}

	//单解区间内循环细化，直至精度满足要求

	for i := 0; i < j; i++ {
		Nn := 0     //死循环约束
		solved := 0 //解得标志

		for {
			Nn += 1
			//循环超过一定数
			if Nn > 1000 {
				err = false
				return sol, err
			}

			h = (ab1[i] - ab0[i]) / float64(N)

			for ii := 1; ii < N+1; ii++ {
				if (fn(ab0[i]+float64(ii)*h) > 0 && fn(ab0[i]+float64(ii-1)*h) < 0) || (fn(ab0[i]+float64(ii)*h) < 0 && fn(ab0[i]+float64(ii-1)*h) > 0) {
					ab0[i] = ab0[i] + float64(ii-1)*h
					ab1[i] = ab0[i] + float64(ii)*h
					//是否满足精度要求
					if math.Abs(fn((ab0[i]+ab1[i])/2.0)) < tol {
						sol[i] = (ab0[i] + ab1[i]) / 2.0
						solved = 1
					}
					break //退出此区间的搜索循环
				}
			}
			//如果解除此区间的解，则退出死循环
			if solved == 1 {
				break
			}
		}
	}
	//返回
	err = true
	return sol, err
}

func BenchmarkSearchByStep(b *testing.B) {
	for i := 0; i < b.N; i++ {
		SearchByStep(func(x float64) float64 { return math.Cos(x) }, -1.0, 5.0, 100, 1e-3)
	}
}
