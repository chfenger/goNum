// ODETrapezoid
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-13
版本   : 0.0.0
------------------------------------------------------
    常微分方程的梯形解法
理论：
    对于常微分方程
     dy
    ---- = f(x, y)
     dx
    y(x0) = y0, x0 <= x

    梯形解法：
                    h
    y_(n+1) = yn + ---(f(xn, yn)+f(x_(n+1), y_(n+1))), n = 0,1,2,3,...
                    2

    梯形法是无条件稳定的
    梯形法为二阶精度的方法

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 181.
------------------------------------------------------
输入   :
    fun     被积分函数
    x0, y0  初值
    h       积分步长
    tol     内循环控制误差
    n       迭代次数
输出   :
    sol     解矩阵，nx2
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// ODETrapezoid 常微分方程的梯形解法
func ODETrapezoid(fun func(float64, float64) float64, x0, y0, h, tol float64, n int) (Matrix, bool) {
	/*
		常微分方程的梯形解法
		输入   :
		    fun     被积分函数
		    x0, y0  初值
		    h       积分步长
		    tol     内循环控制误差
		    n       迭代次数
		输出   :
		    sol     解矩阵，nx2
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if n < 0 {
		panic("Error in goNum.ODETrapezoid: n is not a positive value")
	}

	sol := ZeroMatrix(n+1, 2)
	var err bool = false

	//初值
	sol.SetMatrix(0, 0, x0)
	sol.SetMatrix(0, 1, y0)

	for i := 0; i < n; i++ {
		xi := sol.GetFromMatrix(i, 0)
		yi := sol.GetFromMatrix(i, 1)
		xi10 := xi + h
		yi10 := yi + h*fun(xi, yi)
		//内循环
		yik := make([]float64, 0)
		yik = append(yik, yi10) //k=0
		var k int = 0
		for {
			yik = append(yik, yi+h*(fun(xi, yi)+fun(xi10, yik[k]))/2.0)
			if math.Abs(yik[k+1]-yik[k]) < tol {
				break
			}
			k++
		}

		sol.SetMatrix(i+1, 0, xi10)
		sol.SetMatrix(i+1, 1, yik[k+1])
	}

	err = true
	return sol, err
}
