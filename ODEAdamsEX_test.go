// ODEAdamsEX_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-13
版本   : 0.0.0
------------------------------------------------------
    四步Adams外推公式，显式、线性
理论：
                    h
    y_(n+1) = yn + ----(55f(xn,yn) - 59f(x_(n-1),y_(n-1)) +
                    24
              37f(x_(n-2),y_(n-2)) - 9f(x_(n-3),y_(n-3)))

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 200-201.
------------------------------------------------------
输入   :
    fun     被积分函数
    x0      初值
    xend    积分终止点
    fn      方程个数
    n       迭代次数
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

func ODEAdamsEX(fun func(goNum.Matrix, int) float64, x0 goNum.Matrix, xend float64, fn, n int) (goNum.Matrix, bool) {
	/*
		四步Adams外推公式，显式、线性，单个方程
		输入   :
		    fun     被积分函数
		    x0      初值
		    xend    积分终止点
		    fn      方程个数
		    n       迭代次数
		输出   :
		    sol     解矩阵
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断方程个数是否对应初值个数
	if x0.Rows != fn+1 {
		panic("Error in goNum.ODEAdamsEX: Quantities of x0 and fn+1 are not equal")
	}

	sol := goNum.ZeroMatrix(fn+1, n+1)
	h := (xend - x0.GetFromMatrix(0, 0)) / float64(n)

	//把初值赋给sol
	for i := 0; i < fn+1; i++ {
		sol.SetMatrix(i, 0, x0.Data[i])
	}

	//前三个使用RK44计算，不包括已有的初值点
	xendRK := x0.GetFromMatrix(0, 0) + 3.0*h
	solRK, errRK := goNum.RK44(fun, x0, xendRK, fn, 3)

	if errRK != true {
		panic("Error in goNum.ODEAdamsEX: RK44 solving error")
	}

	//传递RK44计算的结果到sol
	for k := 0; k < fn+1; k++ { //fn个方程，fn+1个参数
		for i := 1; i < 4; i++ { //三个结果
			sol.SetMatrix(k, i, solRK.GetFromMatrix(k, i))
		}
	}

	//Adams外推公式, 4(即n+1)需要3,2,1,0四个
	for i := 4; i < n+1; i++ {
		sol.SetMatrix(0, i, sol.GetFromMatrix(0, i-1)+h) //xi
		//临时初值
		xyn := goNum.ZeroMatrix(fn+1, 1)
		xyn_1 := goNum.ZeroMatrix(fn+1, 1)
		xyn_2 := goNum.ZeroMatrix(fn+1, 1)
		xyn_3 := goNum.ZeroMatrix(fn+1, 1)
		for j := 0; j < fn+1; j++ {
			xyn.Data[j] = sol.GetFromMatrix(j, i-1)
			xyn_1.Data[j] = sol.GetFromMatrix(j, i-2)
			xyn_2.Data[j] = sol.GetFromMatrix(j, i-3)
			xyn_3.Data[j] = sol.GetFromMatrix(j, i-4)
		}
		//计算
		for j := 0; j < fn; j++ { //不包含xi的其他参数
			temp0 := 55.0*fun(xyn, j) - 59.0*fun(xyn_1, j) + 37.0*fun(xyn_2, j) - 9.0*fun(xyn_3, j)
			temp0 = xyn.Data[j+1] + temp0*h/24.0
			sol.SetMatrix(j+1, i, temp0) //yi
		}
	}

	return sol, true
}

func fun40(x0 goNum.Matrix, i int) float64 {
	switch i {
	case 0:
		return x0.Data[1] - 2.0*x0.Data[0]/x0.Data[1]
	default:
		return 0.0
	}

	return 0.0
}

func BenchmarkODEAdamsEX(b *testing.B) {
	x40 := goNum.NewMatrix(2, 1, []float64{0.0, 1.0})
	for i := 0; i < b.N; i++ {
		goNum.ODEAdamsEX(fun40, x40, 1.0, 1, 10)
	}
}
