// PDEDiffEllipticalLL5_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2019-01-07
版本   : 0.0.0
------------------------------------------------------
    求解椭圆型偏微分方程（Laplace）的差分解法（五点格式）
理论：
    对于椭圆型偏微分方程（Laplace方程）：
     d^2u     d^2u
    ------ + ------ = 0
     dx^2     dy^2

    u(x, 0) = fy0(x), u(x, b) = fyb(x)
    u(0, y) = fx0(y), u(a, y) = fxa(y)

    0 < x < a, 0 < y < b

    x分为n等份，y分为m等份

    hy^2[u_(i+1,j) + u_(i-1,j) - 2u_(i,j)] +
       hx^2[u_(i,j+1) + u_(i,j-1) - 2u_(i,j)] = 0

    解以上方程组可得解

    参考 John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 10.3.1.
------------------------------------------------------
输入   :
    funy0, funyb, funx0, funxa   边界函数
    x0      求解范围，2x2
    n, m    网格数量, 对应x和y
输出   :
    sol     解矩阵
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// PDEDiffEllipticalLL5 求解椭圆型偏微分方程（Laplace）的差分解法（五点格式）
func PDEDiffEllipticalLL5(funy0, funyb, funx0, funxa func(float64) float64,
	x0 goNum.Matrix, n, m int) (goNum.Matrix, bool) {
	/*
		求解椭圆型偏微分方程（Laplace）的差分解法（五点格式）
		输入   :
		    funy0, funyb, funx0, funxa   边界函数
		    x0      求解范围，2x2
		    n, m    网格数量, 对应x和y
		输出   :
		    sol     解矩阵
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断网格数量
	if (m < 1) || (n < 1) {
		panic("Error in goNum.PDEDiffEllipticalLL5: Grid numbers error")
	}
	//判断初值维数
	if (x0.Rows < 2) || (x0.Columns < 2) {
		panic("Error in goNum.PDEDiffEllipticalLL5: Initial values error")
	}

	var err bool = false
	sol := goNum.ZeroMatrix(m+1, n+1)                                    //行y变化，列x变化
	hx := (x0.GetFromMatrix(1, 0) - x0.GetFromMatrix(0, 0)) / float64(n) //x方向步长
	hy := (x0.GetFromMatrix(1, 1) - x0.GetFromMatrix(0, 1)) / float64(m) //y方向步长

	//边界框的解
	//第一行的值和最后一行的值，不包括第一个和最后一个
	for i := 1; i < n; i++ {
		sol.SetMatrix(0, i, funy0(x0.GetFromMatrix(0, 0)+hx*float64(i)))
		sol.SetMatrix(m, i, funyb(x0.GetFromMatrix(0, 0)+hx*float64(i)))
	}
	//第一列的值和最后一列的值，包括第一个和最后一个
	for j := 0; j < m+1; j++ {
		sol.SetMatrix(j, 0, funx0(x0.GetFromMatrix(0, 1)+hy*float64(j)))
		sol.SetMatrix(j, n, funxa(x0.GetFromMatrix(0, 1)+hy*float64(j)))
	}

	//求解中间点，主对角占优矩阵解法，利用高斯消去方法
	AA := goNum.ZeroMatrix((n-1)*(m-1), (n-1)*(m-1)) //系数矩阵A
	BA := goNum.ZeroMatrix((n-1)*(m-1), 1)           //值矩阵B

	//赋值系数矩阵和值矩阵
	//第一行, j = 1
	//第一个
	AA.SetMatrix(0, 0, -2.0*hx*hx-2.0*hy*hy)
	AA.SetMatrix(0, 1, hy*hy)
	AA.SetMatrix(0, n-1, hx*hx)
	tempBA := 0.0 - hy*hy*funx0(x0.GetFromMatrix(0, 1)+hy*1.0)
	tempBA = tempBA - hx*hx*funy0(x0.GetFromMatrix(0, 0)+hx*1)
	BA.SetMatrix(0, 0, tempBA)
	for i := 2; i < n-1; i++ {
		AA.SetMatrix(i-1, i-2, hy*hy)
		AA.SetMatrix(i-1, i-1, -2.0*hx*hx-2.0*hy*hy)
		AA.SetMatrix(i-1, i, hy*hy)
		AA.SetMatrix(i-1, (n-1)*1+i-1, hx*hx)
		tempBA = 0.0 - hx*hx*funy0(x0.GetFromMatrix(0, 0)+hx*float64(i))
		BA.SetMatrix((n-1)*0+i-1, 0, tempBA)
	}
	//最后一个
	AA.SetMatrix(n-1-1, n-1-2, hy*hy)
	AA.SetMatrix(n-1-1, n-1-1, -2.0*hx*hx-2.0*hy*hy)
	AA.SetMatrix(n-1-1, (n-1)*1+n-1-1, hx*hx)
	tempBA = 0.0 - hy*hy*funxa(x0.GetFromMatrix(0, 1)+hy*1.0)
	tempBA = tempBA - hx*hx*funy0(x0.GetFromMatrix(0, 0)+hx*float64(n-1))
	BA.SetMatrix(n-1-1, 0, tempBA)
	//中间行
	for j := 2; j < m-1; j++ {
		//第一个
		AA.SetMatrix((n-1)*(j-1), (n-1)*(j-1-1), hx*hx)
		AA.SetMatrix((n-1)*(j-1), (n-1)*(j-1), -2.0*hx*hx-2.0*hy*hy)
		AA.SetMatrix((n-1)*(j-1), (n-1)*(j-1)+1, hy*hy)
		AA.SetMatrix((n-1)*(j-1), (n-1)*(j-1)+n-1, hx*hx)
		tempBA = 0.0 - hy*hy*funx0(x0.GetFromMatrix(0, 1)+hy*float64(j))
		BA.SetMatrix((n-1)*(j-1), 0, tempBA)
		for i := 2; i < n-1; i++ {
			AA.SetMatrix((n-1)*(j-1)+i-1, (n-1)*(j-1-1)+i-1, hx*hx)
			AA.SetMatrix((n-1)*(j-1)+i-1, (n-1)*(j-1)+i-2, hy*hy)
			AA.SetMatrix((n-1)*(j-1)+i-1, (n-1)*(j-1)+i-1, -2.0*hx*hx-2.0*hy*hy)
			AA.SetMatrix((n-1)*(j-1)+i-1, (n-1)*(j-1)+i, hy*hy)
			AA.SetMatrix((n-1)*(j-1)+i-1, (n-1)*(j-1+1)+i-1, hx*hx)
			BA.SetMatrix((n-1)*(j-1)+i-1, 0, 0.0)
		}
		//最后一个
		AA.SetMatrix((n-1)*(j-1)+n-1-1, (n-1)*(j-1-1)+n-1-1, hx*hx)
		AA.SetMatrix((n-1)*(j-1)+n-1-1, (n-1)*(j-1)+n-1-2, hy*hy)
		AA.SetMatrix((n-1)*(j-1)+n-1-1, (n-1)*(j-1)+n-1-1, -2.0*hx*hx-2.0*hy*hy)
		AA.SetMatrix((n-1)*(j-1)+n-1-1, (n-1)*(j-1+1)+n-1-1, hx*hx)
		tempBA = 0.0 - hy*hy*funxa(x0.GetFromMatrix(0, 1)+hy*float64(j))
		BA.SetMatrix((n-1)*(j-1)+n-1-1, 0, tempBA)
	}
	//最后一行, j = m-1
	//第一个
	AA.SetMatrix((n-1)*(m-1-1), (n-1)*(m-1-1-1), hx*hx)
	AA.SetMatrix((n-1)*(m-1-1), (n-1)*(m-1-1), -2.0*hx*hx-2.0*hy*hy)
	AA.SetMatrix((n-1)*(m-1-1), (n-1)*(m-1-1)+1, hy*hy)
	tempBA = 0.0 - hy*hy*funx0(x0.GetFromMatrix(0, 1)+hy*float64(m-1))
	tempBA = tempBA - hx*hx*funyb(x0.GetFromMatrix(0, 0)+hx*1)
	BA.SetMatrix((n-1)*(m-1-1), 0, tempBA)
	for i := 2; i < n-1; i++ {
		AA.SetMatrix((n-1)*(m-1-1)+i-1, (n-1)*(m-1-1-1)+i-1, hx*hx)
		AA.SetMatrix((n-1)*(m-1-1)+i-1, (n-1)*(m-1-1)+i-2, hy*hy)
		AA.SetMatrix((n-1)*(m-1-1)+i-1, (n-1)*(m-1-1)+i-1, -2.0*hx*hx-2.0*hy*hy)
		AA.SetMatrix((n-1)*(m-1-1)+i-1, (n-1)*(m-1-1)+i, hy*hy)
		tempBA = 0.0 - hx*hx*funyb(x0.GetFromMatrix(0, 0)+hx*float64(i))
		BA.SetMatrix((n-1)*(m-1-1)+i-1, 0, tempBA)
	}
	//最后一个
	AA.SetMatrix((n-1)*(m-1-1)+n-1-1, (n-1)*(m-1-1-1)+n-1-1, hx*hx)
	AA.SetMatrix((n-1)*(m-1-1)+n-1-1, (n-1)*(m-1-1)+n-1-2, hy*hy)
	AA.SetMatrix((n-1)*(m-1-1)+n-1-1, (n-1)*(m-1-1)+n-1-1, -2.0*hx*hx-2.0*hy*hy)
	tempBA = 0.0 - hy*hy*funxa(x0.GetFromMatrix(0, 1)+hy*float64(m-1))
	tempBA = tempBA - hx*hx*funyb(x0.GetFromMatrix(0, 0)+hx*float64(n-1))
	BA.SetMatrix((n-1)*(m-1-1)+n-1-1, 0, tempBA)
	//求解矩阵方程
	tempp, temperr := goNum.LEs_ECPE(goNum.Matrix2ToSlices(AA), goNum.Matrix1ToSlices(BA))
	if temperr != true {
		panic("Error in goNum.PDEDiffEllipticalLL5: Solve error")
	}
	//解赋予sol
	ii := 1
	jj := 1
	for i := 0; i < len(tempp); i++ {
		sol.SetMatrix(ii, jj, tempp[i])
		if jj == n-1 {
			ii++
			jj = 0
		}
		jj++
	}

	err = true
	return sol, err
}

func fun55y0(x float64) float64 {
	return 20.0
}

func fun55yb(x float64) float64 {
	return 180.0
}

func fun55x0(y float64) float64 {
	return 80.0
}

func fun55xa(Y float64) float64 {
	return 0.0
}

func BenchmarkPDEDiffEllipticalLL5(b *testing.B) {
	x55 := goNum.NewMatrix(2, 2, []float64{0.0, 0.0, 4.0, 4.0})
	for i := 0; i < b.N; i++ {
		goNum.PDEDiffEllipticalLL5(fun55y0, fun55yb, fun55x0, fun55xa, x55, 4, 4)
	}
}
