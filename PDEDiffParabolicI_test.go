// PDEDiffParabolicI_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-14
版本   : 0.0.0
------------------------------------------------------
    求解抛物型偏微分方程的差分解法（隐式）
理论：
    对于抛物型偏微分方程：
     du       d^2u
    ---- = A ------ + B
     dt       dx^2

    u(x, 0) = p(x)
    u(0, t) = u1(t), u(L, t) = u2(t)

    0 < x < L, 0 < t < T

    则古典隐式差分格式为，x分为m等份，t分为n等份
    Au_(j+1) = uj + F_(j+1)

        |1+2l -l              |
        |-l   1+2l -l         |
    A = |      ..........     |
        |         -l 1+2l -l  |
        |            -l   1+2l|

    u_(j+1) = [u_(1,j+1),u_(2,j+1),...,u_(m-1,j+1)]'
    F_(j+1) = [lu1((j+1)*tau)+B*tau,B*tau,B*tau,...,B*tau,lu2((j+1)*tau)+B*tau]'
    V_(j+1) = uj + F_(j+1)
    j = 0,1,...,n-1

    u0 = [u_(1,0),u_(2,0),...,u_(m-1,0)]'
       = [p(h),p(2h),...,p((m-1)h)]'


    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 214-215.
------------------------------------------------------
输入   :
    funp, funu1, funu2   边界函数
    x0      求解范围，2x2
    A, B    常系数
    m, n    网格数量
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

func PDEDiffParabolicI(funp, funu1, funu2 func(float64) float64, x0 goNum.Matrix, A, B float64, m, n int) (goNum.Matrix, bool) {
	/*
	   求解抛物型偏微分方程的差分解法（隐式）
	   输入   :
	       funp, funu1, funu2   边界函数
	       x0      求解范围，2x2
	       A, B    常系数
	       m, n    网格数量
	   输出   :
	       sol     解矩阵
	       err     解出标志：false-未解出或达到步数上限；
	                        true-全部解出
	*/
	//判断网格数量
	if (m < 1) || (n < 1) {
		panic("Error in goNum.PDEDiffParabolicI: Grid numbers error")
	}

	var err bool = false
	sol := goNum.ZeroMatrix(m+1, n+1)
	hx := (x0.GetFromMatrix(1, 0) - x0.GetFromMatrix(0, 0)) / float64(m) //x方向步长
	ht := (x0.GetFromMatrix(1, 1) - x0.GetFromMatrix(0, 1)) / float64(n) //t方向步长

	//1. 计算t第零层上的值u_(i,0) i=0,1,...,m
	for i := 0; i < m+1; i++ {
		sol.SetMatrix(i, 0, funp(x0.GetFromMatrix(0, 0)+float64(i)*hx))
	}
	//2. 计算左右边界上的节点u_(0,j)和u_(m,j) j=1,2,...,n
	for j := 1; j < n+1; j++ {
		sol.SetMatrix(0, j, funu1(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //左边界
		sol.SetMatrix(m, j, funu2(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //右边界
	}

	l := A * ht / (hx * hx)
	//稳定性判断
	if l <= 0 {
		panic("Error in goNum.PDEDiffParabolicS: lambda less than or equal to zero")
	}
	//A赋值
	AA := goNum.ZeroMatrix(m-1, m-1)
	ui := goNum.ZeroMatrix(m-1, 1)
	Fi := goNum.ZeroMatrix(m-1, 1)
	AA.SetMatrix(0, 0, 1.0+2.0*l) //第零行
	AA.SetMatrix(0, 1, -1.0*l)
	ui.Data[0] = sol.GetFromMatrix(1, 0)
	for i := 1; i < m-2; i++ {
		AA.SetMatrix(i, i-1, -1.0*l)
		AA.SetMatrix(i, i, 1.0+2.0*l)
		AA.SetMatrix(i, i+1, -1.0*l)
		ui.Data[i] = sol.GetFromMatrix(i+1, 0)
		Fi.Data[i] = B * ht
	}
	AA.SetMatrix(m-2, m-3, -1.0*l) //第零行
	AA.SetMatrix(m-2, m-2, 1.0+2.0*l)
	ui.Data[m-2] = sol.GetFromMatrix(m-1, 0)
	//内部节点循环求解
	for j := 0; j < n; j++ {
		//F，每一步需要重新计算第一项和最后一项
		Fi.Data[0] = l*funu1(float64(j+1)*ht) + B*ht
		Fi.Data[m-2] = l*funu2(float64(j+1)*ht) + B*ht
		//
		ui1, errtemp := goNum.LEs_Chasing(AA, goNum.AddMatrix(ui, Fi))
		if errtemp != true {
			panic("Error in goNum.PDEDiffParabolicI: Chasing solved error")
		}
		for i := 0; i < m-1; i++ {
			ui.Data[i] = ui1.Data[i]
			sol.SetMatrix(i+1, j+1, ui1.Data[i])
		}
	}

	err = true
	return sol, err
}

func fun41_1_p(x float64) float64 {
	return 4.0 * x * (1.0 - x)
}

func fun41_1_u1(t float64) float64 {
	return 0.0
}

func fun41_1_u2(t float64) float64 {
	return 0.0
}

func BenchmarkPDEDiffParabolicI(b *testing.B) {
	x41_1 := goNum.NewMatrix(2, 2, []float64{0.0, 0.0, 1.0, 0.06})
	for i := 0; i < b.N; i++ {
		goNum.PDEDiffParabolicI(fun41_1_p, fun41_1_u1, fun41_1_u2, x41_1, 1.0, 0.0, 10, 36)
	}
}
