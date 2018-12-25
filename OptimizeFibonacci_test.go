// OptimizeFibonacci_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-24
版本   : 0.0.0
------------------------------------------------------
    Fibonacci搜索法求单峰单自变量极小值
理论：
    对于在区间[a, b]内有定义的凹函数f(x)，取点：
    ck = ak+(1-r)(bk-ak)
    d = ak+rk(bk-ak)
    其中r为Fibonacci数列值之比F_(n-k-1)/F_(n-k)

    迭代次数n应使得Fn > (b0-a0)/tol

    如果f(c) <= f(d)，则将d赋予b，c赋予d，继续迭代；
    如果f(c) > f(d)，则将c赋予a，d赋予c，继续迭代。
    迭代终止条件为Abs(f(a)-f(b)) < tol，取区间中值

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 8.1.1.2，并改进
------------------------------------------------------
输入   :
    fun     函数
    a, b    区间范围
    tol     控制误差
输出   :
    sol     解
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

func OptimizeFibonacci(fun func(float64) float64, a, b, tol float64) (float64, bool) {
	/*
		Fibonacci搜索法求单峰单自变量极小值
		输入   :
		    fun     函数
		    a, b    区间范围
		    tol     控制误差
		输出   :
		    sol     解
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断a和b的关系
	if math.Abs(fun(a)-fun(b)) < tol {
		if fun(a) < fun(b) {
			return a, true
		} else {
			return b, true
		}
	}

	var sol float64
	var err bool = false
	var n, cdFlag int = 0, 0 //cdFlag---下一步计算c（cdFlag=0）还是d（cdFlag=1）

	//计算n
	bat := (fun(b) - fun(a)) / tol
	for i := 0; i < 1e6; i++ {
		if float64(goNum.Fibonacci(i)) > bat {
			n = i
			break
		}
	}

	//计算
	//第一步计算两次，c、d
	fnn := float64(goNum.Fibonacci(n-1)) / float64(goNum.Fibonacci(n))
	ba := b - a
	c := a + (1.0-fnn)*ba
	d := a + fnn*ba
	fc := fun(c)
	fd := fun(d)
	if fc <= fd {
		b = d
		d = c
		fd = fc
		cdFlag = 0
	} else {
		a = c
		c = d
		fc = fd
		cdFlag = 1
	}
	//0 < k < n-3
	for k := 1; k < n-3; k++ {
		fnn = float64(goNum.Fibonacci(n-k-1)) / float64(goNum.Fibonacci(n-k))
		ba = b - a
		if cdFlag == 0 { //计算c
			c = a + (1.0-fnn)*ba
			fc = fun(c)
		} else { //计算d
			d = a + fnn*ba
			fd = fun(d)
		}
		//下一步
		if fc <= fd {
			b = d
			d = c
			fd = fc
			cdFlag = 0
		} else {
			a = c
			c = d
			fc = fd
			cdFlag = 1
		}
	}
	//k=n-3, F2/F3 = 1/2, 不放入循环是为减少if判断的损耗
	fnn = 0.5 - 0.01 //加区别常数0.01
	ba = b - a
	if cdFlag == 0 { //计算c
		c = a + (1.0-fnn)*ba
		fc = fun(c)
	} else { //计算d
		d = a + fnn*ba
		fd = fun(d)
	}
	if fc <= fd {
		b = d
	} else {
		a = c
	}
	sol = (b + a) / 2.0

	err = true
	return sol, err
}

func fun50f(x float64) float64 {
	return x*x - math.Sin(x)
}

func BenchmarkOptimizeFibonacci(b *testing.B) {
	for i := 0; i < b.N; i++ {
		goNum.OptimizeFibonacci(fun50f, 0.0, 1.0, 1e-10)
	}
}
