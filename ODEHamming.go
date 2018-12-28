// ODEHamming
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-26
版本   : 0.0.0
------------------------------------------------------
    Hamming预估校正方法
理论：
    预估：
                         4h
    p_(k+1) = y_(k-3) + ---(2f_(k-2)-f_(k-1)+2fk)
                         3
    校正）：
               -y_(k-2)+9yk     3h
    y_(k+1) = -------------- + ---(-f_(k-1)+2fk+f_(k+1))
                     8          8

    步长 h < 0.69/|fy(x,y)|

    四阶精度

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 9.6.6
------------------------------------------------------
输入   :
    fun     被积分函数
    x0      初值,2x4
    h       步长
    n       积分步数
输出   :
    sol     解矩阵
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

func ODEHamming(fun func(float64, float64) float64, x0 Matrix, h float64, n int) (Matrix, bool) {
	/*
		Hamming预估校正方法
		输入   :
		    fun     被积分函数
		    x0      初值,2x4
		    h       步长
		    n       积分步数
		输出   :
		    sol     解矩阵
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if n < 0 {
		panic("Error in goNum.Hamming: n is not a positive value")
	}
	//判断初值
	if (x0.Rows != 2) || (x0.Columns < 4) {
		panic("Error in goNum.Hamming: Initial values error")
	}

	sol := ZeroMatrix(2, n+1)
	p := ZeroMatrix(n+1, 1)
	var err bool = false

	//初值
	for i := 0; i < 4; i++ {
		sol.SetMatrix(0, i, x0.GetFromMatrix(0, i))
		sol.SetMatrix(1, i, x0.GetFromMatrix(1, i))
	}

	//计算
	for i := 4; i < n+1; i++ {
		sol.SetMatrix(0, i, sol.GetFromMatrix(0, i-1)+h) //xi
		//pi
		temp0 := fun(sol.GetFromMatrix(0, i-3), sol.GetFromMatrix(1, i-3)) //f_(i-3)
		temp1 := fun(sol.GetFromMatrix(0, i-2), sol.GetFromMatrix(1, i-2)) //f_(i-2)
		temp2 := fun(sol.GetFromMatrix(0, i-1), sol.GetFromMatrix(1, i-1)) //f_(i-1)
		soltemp := 2.0 * temp0
		soltemp += -1.0 * temp1
		soltemp += 2.0 * temp2
		p.SetMatrix(i, 0, sol.GetFromMatrix(1, i-4)+4.0*h*soltemp/3.0)
		//yi
		soltemp = -1.0 * temp1
		soltemp += 2.0 * temp2
		soltemp += fun(sol.GetFromMatrix(0, i), p.GetFromMatrix(i, 0)) //fi
		soltemp = 3.0 * h * soltemp / 8.0
		soltemp += (-1.0*sol.GetFromMatrix(1, i-3) + 9.0*sol.GetFromMatrix(1, i-1)) / 8.0
		sol.SetMatrix(1, i, soltemp)
	}

	err = true
	return sol, err
}
