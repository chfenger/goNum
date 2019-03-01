// IntegralCompositeNewtonCotesHalf_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-12
版本   : 0.0.0
------------------------------------------------------
    1-8级逐次分半复化Newton-Cotes求积分公式
理论：
    对于积分
    b           n
    |f(x)dx ~= Sum Ak*f(xk)
    a          k=0

               (n)
    Ak = (b-a)C
               k

    (n)   (-1)^(n-k)  n
   C   = ------------ |t(t-1)(t-2)...(t-(k-1))(t-(k+1))...(t-n)dt
    k     k!(n-k)!n   0

    特别的，n=1为复化梯形公式；
           n=2为复化Simpson（辛浦生）公式；
           n=4为复化Cotes（科特斯）公式

    将区间[a, b]等分为Nn个子区间，每个子区间上使用Newton-Cotes求积分公式
    子区间逐次分半以满足精度要求

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 159-160.
------------------------------------------------------
输入   :
    fun     被积分函数
    a, b    积分范围
    tol     精度要求
    n       Newton-Cotes公式级数
    Nmax    最大循环次数
输出   :
    sol     解
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
注意   ：
    由于误差得不到有效控制，稳定性无法保证，故而并不是n值越
    大越好，实际应用中很少使用n值较大的Newton-Cotes公式
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

// IntegralCompositeNewtonCotesHalf 1-8级逐次分半复化Newton-Cotes求积分公式
func IntegralCompositeNewtonCotesHalf(fun func(float64) float64, a, b, tol float64, n, Nmax int) (float64, bool) {
	/*
		1-8级逐次分半复化Newton-Cotes求积分公式
		输入   :
		    fun     被积分函数
		    a, b    积分范围
		    tol     精度要求
		    n       Newton-Cotes公式级数
		    Nmax    最大循环次数
		输出   :
		    sol     解
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if (n < 1) || (n > 8) {
		panic("Error in goNum.IntegralNewtonCotes: n is not correct")
	}
	//判断a, b
	if a == b {
		return 0.0, true
	}

	var sol0, sol1 float64
	var err bool = false
	var Nn int = 1

	//循环
	for k := 0; k < Nmax; k++ {
		sol1 = 0.0
		//子区间长度
		Hh := (b - a) / float64(Nn)
		//调用IntegralNewtonCotes循环累加
		for i := 1; i < Nn+1; i++ {
			soltemp, errtemp := goNum.IntegralNewtonCotes(fun, a+Hh*float64(i-1), a+Hh*float64(i), n)
			if errtemp != true {
				panic("Error in goNum.IntegralNewtonCotes: Error in calling IntegralNewtonCotes")
			}
			sol1 += soltemp
		}
		if math.Abs(sol1-sol0) < tol {
			return sol1, true
		}
		Nn = 2 * Nn
		sol0 = sol1
	}

	return 0.0, err
}

func fun36(x float64) float64 {
	return 4.0 / (1 + x*x)
}

func BenchmarkIntegralCompositeNewtonCotesHalf(b *testing.B) {
	for i := 0; i < b.N; i++ {
		goNum.IntegralCompositeNewtonCotesHalf(fun36, 0.0, 1.0, 1e-6, 1, 1e3)
	}
}
