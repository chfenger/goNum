// Muller_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-20
版本   : 0.0.0
------------------------------------------------------
    Muller法求解非线性方程f(x)=0的解
理论：

    参考 John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 2.5.2.
------------------------------------------------------
输入   :
    fun     求解函数
    x0      初值自变量，三个不同点，3x1
    tol     控制误差
    n       最大迭代步数
输出   :
    sol     解
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

// Muller Muller法求解非线性方程f(x)=0的解
func Muller(fun func(float64) float64, x0 goNum.Matrix, tol float64, n int) (float64, bool) {
	/*
	   Muller法求解非线性方程f(x)=0的解
	   输入   :
	       fun     求解函数
	       x0      初值自变量，三个不同点，3x1
	       tol     控制误差
	       n       最大迭代步数
	   输出   :
	       sol     解
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断tol
	if tol <= 0.0 {
		panic("Error in goNum.Muller: tol less than or euqals to zero")
	}

	var sol float64
	var err bool = false

	//x0赋给p0并计算对应的y0
	p0 := goNum.ZeroMatrix(x0.Rows, x0.Columns+1)
	x0sort, _ := goNum.MinMaxSort(x0.Data)
	for i := 0; i < x0.Rows; i++ {
		p0.SetMatrix(i, 0, x0sort[i])
		p0.SetMatrix(i, 1, fun(x0sort[i]))
	}

	//迭代计算
	for i := 0; i < n; i++ {
		//准备系数
		h0 := p0.GetFromMatrix(0, 0) - p0.GetFromMatrix(2, 0)
		h1 := p0.GetFromMatrix(1, 0) - p0.GetFromMatrix(2, 0)
		c := p0.GetFromMatrix(2, 1)
		e0 := p0.GetFromMatrix(0, 1) - c
		e1 := p0.GetFromMatrix(1, 1) - c
		b := h1*h0*h0 - h0*h1*h1
		a := (e0*h1 - e1*h0) / b
		b = (e1*h0*h0 - e0*h1*h1) / b
		//求根
		z2 := b*b - 4.0*a*c
		if z2 < 0 {
			//panic("Error in goNum.Muller: There is complex values exist")
			z2 = 0
		}

		var z float64
		if b < 0 {
			z = -2.0 * c / (b - math.Sqrt(z2))
		}
		z = -2.0 * c / (b + math.Sqrt(z2))
		z = p0.GetFromMatrix(2, 0) + z

		//判断解
		if math.Abs(fun(z)) < tol {
			err = true
			sol = z
			return sol, err
		}

		//删除离z最远的点
		dis := []float64{
			z - p0.GetFromMatrix(0, 0),
			z - p0.GetFromMatrix(1, 0),
			z - p0.GetFromMatrix(2, 0)}
		_, deli, _ := goNum.MaxAbs(dis)
		for j := 0; j < 3; j++ {
			if deli == j {
				p0.SetMatrix(j, 0, z)
			}
		}
		x0sort, _ = goNum.MinMaxSort(p0.ColumnOfMatrix(0))
		for j := 0; j < 3; j++ {
			p0.SetMatrix(j, 0, x0sort[j])
			p0.SetMatrix(j, 1, fun(x0sort[j]))
		}
	}

	err = false
	return sol, err
}

func fun44(x float64) float64 {
	return x*x*x - 3.0*x + 2.0
}

func BenchmarkMuller(b *testing.B) {
	x44 := goNum.NewMatrix(3, 1, []float64{2.0, 3.0, 4.0})
	for i := 0; i < b.N; i++ {
		goNum.Muller(fun44, x44, 1e-12, 1e6)
	}
}
