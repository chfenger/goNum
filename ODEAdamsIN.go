// ODEAdamsIN
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-13
版本   : 0.0.0
------------------------------------------------------
    三次Adams内插公式，隐式、线性
理论：
                    h
    y_(n+1) = yn + ----(9f(x_(n+1),y_(n+1)) + 19f(xn,yn) -
                    24
              5f(x_(n-1),y_(n-1)) + f(x_(n-2),y_(n-2)))

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 201-202.
------------------------------------------------------
输入   :
    fun     被积分函数
    x0      初值
    xend    积分终止点
    tol     内迭代控制误差
    fn      方程个数
    n       迭代次数
输出   :
    sol     解矩阵
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// ODEAdamsIN 四步Adams外推公式，显式、线性，单个方程
func ODEAdamsIN(fun func(Matrix, int) float64, x0 Matrix, xend, tol float64, fn, n int) (Matrix, bool) {
	/*
		四步Adams外推公式，显式、线性，单个方程
		输入   :
		    fun     被积分函数
		    x0      初值
		    xend    积分终止点
		    tol     内迭代控制误差
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

	sol := ZeroMatrix(fn+1, n+1)
	h := (xend - x0.GetFromMatrix(0, 0)) / float64(n)

	//把初值赋给sol
	for i := 0; i < fn+1; i++ {
		sol.SetMatrix(i, 0, x0.Data[i])
	}

	//前三个使用RK44计算，不包括已有的初值点
	xendRK := x0.GetFromMatrix(0, 0) + 3.0*h
	solRK, errRK := RK44(fun, x0, xendRK, fn, 3)

	if errRK != true {
		panic("Error in goNum.ODEAdamsEX: RK44 solving error")
	}

	//传递RK44计算的结果到sol
	for k := 0; k < fn+1; k++ { //fn个方程，fn+1个参数
		for i := 1; i < 4; i++ { //三个结果
			sol.SetMatrix(k, i, solRK.GetFromMatrix(k, i))
		}
	}

	//三次Adams内插公式, 4(即n+1)需要3,2,1,0四个
	for i := 4; i < n+1; i++ {
		sol.SetMatrix(0, i, sol.GetFromMatrix(0, i-1)+h) //xi
		//临时初值
		xyn := ZeroMatrix(fn+1, 1)
		xyn_1 := ZeroMatrix(fn+1, 1)
		xyn_2 := ZeroMatrix(fn+1, 1)
		xyn_3 := ZeroMatrix(fn+1, 1)
		xyn10 := ZeroMatrix(fn+1, 1)
		for j := 0; j < fn+1; j++ {
			xyn.Data[j] = sol.GetFromMatrix(j, i-1)
			xyn_1.Data[j] = sol.GetFromMatrix(j, i-2)
			xyn_2.Data[j] = sol.GetFromMatrix(j, i-3)
			xyn_3.Data[j] = sol.GetFromMatrix(j, i-4)
		}
		xyn10.Data[0] = sol.GetFromMatrix(0, i) //x_(n+1)
		//内插公式隐式迭代初值，为4步Adams外推公式结果
		for j := 0; j < fn; j++ { //不包含xi的其他参数
			temp0 := 55.0*fun(xyn, j) - 59.0*fun(xyn_1, j) + 37.0*fun(xyn_2, j) - 9.0*fun(xyn_3, j)
			xyn10.Data[j+1] = xyn.Data[j+1] + temp0*h/24.0 //y_(n+1)0
		}
		//内插，对每个公式
		for j := 0; j < fn; j++ {
			//隐式迭代,误差控制
			yn1k := xyn10.Data[j+1]
			for {
				temp0 := 9.0*fun(xyn10, j) + 19.0*fun(xyn, j) - 5.0*fun(xyn_1, j) + fun(xyn_2, j)
				temp0 = xyn.Data[j+1] + h*temp0/24.0
				if math.Abs(temp0-yn1k) < tol {
					sol.SetMatrix(j+1, i, temp0)
					break //跳出无条件循环
				}
				yn1k = temp0
			}
		}
	}

	return sol, true
}
