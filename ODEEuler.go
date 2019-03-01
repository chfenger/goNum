// ODEEuler
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-13
版本   : 0.0.0
------------------------------------------------------
    常微分方程的Euler（欧拉）解法
理论：
    对于常微分方程
     dy
    ---- = f(x, y)
     dx
    y(x0) = y0, x0 <= x

    Euler（欧拉）解法：
    y_(n+1) = yn + hf(xn, yn), n = 0,1,2,3,...

    欧拉法是条件稳定的： 0 <= h <=-2.0/(y'/y)
    欧拉法为一阶精度的方法

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 179.
------------------------------------------------------
输入   :
    fun     被积分函数
    x0, y0  初值
    h       积分步长
    n       迭代次数
输出   :
    sol     解矩阵，nx2
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// ODEEuler 常微分方程的Euler（欧拉）解法
func ODEEuler(fun func(float64, float64) float64, x0, y0, h float64, n int) (Matrix, bool) {
	/*
		常微分方程的Euler（欧拉）解法
		输入   :
		    fun     被积分函数
		    x0, y0  初值
		    h       积分步长
		    n       迭代次数
		输出   :
		    sol     解矩阵，nx2
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if n < 0 {
		panic("Error in goNum.ODEEuler: n is not a positive value")
	}

	sol := ZeroMatrix(n+1, 2)
	var err bool = false

	//初值
	sol.SetMatrix(0, 0, x0)
	sol.SetMatrix(0, 1, y0)

	for i := 0; i < n; i++ {
		xi := sol.GetFromMatrix(i, 0)
		xi1 := xi + h
		yi1 := sol.GetFromMatrix(i, 1) + h*fun(xi, sol.GetFromMatrix(i, 1))
		sol.SetMatrix(i+1, 0, xi1)
		sol.SetMatrix(i+1, 1, yi1)
	}

	err = true
	return sol, err
}
