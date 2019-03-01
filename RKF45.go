// RKF45
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-19
版本   : 0.0.0
------------------------------------------------------
    四级五阶变步长Runge-Kutta法求解常微分方程组
理论：

    参考 John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 9.5.4.
------------------------------------------------------
输入   :
    fun     第i个方程(计算变量值向量, i)
    x0      初值向量，(fn+1)x1，一个x，fn个因变量
    xend    终止x
    tol     步长控制误差
    fn      方程个数
    n       最大迭代步数
输出   :
    sol     解向量
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// RKF45 四级五阶变步长Runge-Kutta法求解常微分方程组
func RKF45(fun func(Matrix, int) float64, x0 Matrix,
	xend, tol float64, fn, n int) (Matrix, bool) {
	/*
		四级五阶变步长Runge-Kutta法求解常微分方程组
		输入   :
		    fun     第i个方程(计算变量值向量, i)
		    x0      初值向量，(fn+1)x1，一个x，fn个因变量
		    xend    终止x
		    tol     步长控制误差
		    fn      方程个数
		    n       最大迭代步数
		输出   :
		    sol     解向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断方程个数是否对应初值个数
	if x0.Rows != fn+1 {
		panic("Error in goNum.RKF45: Quantities of x0 and fn+1 are not equal")
	}
	//判断tol值
	if tol <= 0.0 {
		panic("Error in goNum.RKF45: tol less than or euqals to zero")
	}
	//判断xend值
	if xend <= x0.Data[0] {
		panic("Error in goNum.RKF45: xend less than or euqals to x0")
	}

	sol0 := ZeroMatrix(fn+1, n+1)
	var err bool = false
	h := 100.0 * (xend - x0.Data[0]) / float64(n) //初始步长，100倍最小步长，可修改

	//把初值赋给sol
	for i := 0; i < fn+1; i++ {
		sol0.SetMatrix(i, 0, x0.Data[i])
	}

	//nreal解矩阵实际长度
	var i, nreal int = 1, 1

	for sol0.GetFromMatrix(0, i-1) < xend { //最大迭代次数控制
		temp0 := ZeroMatrix(fn+1, 1)
		//给temp0赋i-1步值，每一步开始
		for j := 0; j < fn+1; j++ {
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1)
		}
		k1 := ZeroMatrix(fn, 1)
		k2 := ZeroMatrix(fn, 1)
		k3 := ZeroMatrix(fn, 1)
		k4 := ZeroMatrix(fn, 1)
		k5 := ZeroMatrix(fn, 1)
		k6 := ZeroMatrix(fn, 1)
		//1. k1
		for j := 0; j < fn; j++ { //微分方程迭代
			k1.Data[j] = h * fun(temp0, j)
		}
		//2. k2
		temp0.Data[0] = sol0.GetFromMatrix(0, i-1) + h/4.0 //xn+h/4
		for j := 1; j < fn+1; j++ {                        //yn+k1/4
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1) + k1.Data[j-1]/4.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k2.Data[j] = h * fun(temp0, j)
		}
		//3. k3
		temp0.Data[0] = sol0.GetFromMatrix(0, i-1) + 3.0*h/8.0 //xn+3h/8
		for j := 1; j < fn+1; j++ {                            //yn+3k1/32+9k2/32
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1) + 3.0*k1.Data[j-1]/32.0 +
				9.0*k2.Data[j-1]/32.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k3.Data[j] = h * fun(temp0, j)
		}
		//4. k4
		temp0.Data[0] = sol0.GetFromMatrix(0, i-1) + 12.0*h/13.0 //xn+12h/13
		for j := 1; j < fn+1; j++ {                              //yn+1932k1/2197-7200k2/2197+7296k3/2197
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1) + 1932.0*k1.Data[j-1]/2197.0 -
				7200.0*k2.Data[j-1]/2197.0 + 7296.0*k3.Data[j-1]/2197.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k4.Data[j] = h * fun(temp0, j)
		}
		//5. k5
		temp0.Data[0] = sol0.GetFromMatrix(0, i-1) + h //xn+h
		for j := 1; j < fn+1; j++ {                    //yn+439k1/216-8k2+3680k3/513-845k4/4104
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1) + 439.0*k1.Data[j-1]/216.0 -
				8.0*k2.Data[j-1] + 3680.0*k3.Data[j-1]/513.0 - 845.0*k4.Data[j-1]/4104.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k5.Data[j] = h * fun(temp0, j)
		}
		//6. k6
		temp0.Data[0] = sol0.GetFromMatrix(0, i-1) + h/2.0 //xn+h/2
		for j := 1; j < fn+1; j++ {                        //yn-8k1/27+2k2-3544k3/2565+1859k4/4104-11k5/40
			temp0.Data[j] = sol0.GetFromMatrix(j, i-1) - 8.0*k1.Data[j-1]/27.0 +
				2.0*k2.Data[j-1] - 3544.0*k3.Data[j-1]/2565.0 + 1859.0*k4.Data[j-1]/4104.0 -
				11.0*k5.Data[j-1]/40.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k6.Data[j] = h * fun(temp0, j)
		}

		//误差与步长
		errtemp := ZeroMatrix(fn, 1) //=ABS(z_(k+1)-y_(k+1))
		for j := 1; j < fn+1; j++ {
			errtemp.Data[j-1] = k1.Data[j-1]/360.0 - 128.0*k3.Data[j-1]/4275.0 -
				2197.0*k4.Data[j-1]/75240.0 + k5.Data[j-1]/50.0 + 2.0*k6.Data[j-1]/55.0
		}
		errtemp0, _, _ := MaxAbs(errtemp.Data)

		//正常推进
		if math.Abs(errtemp0) < tol {
			//i步值
			sol0.SetMatrix(0, i, sol0.GetFromMatrix(0, i-1)+h) //xi
			for j := 1; j < fn+1; j++ {
				soltemp1 := 25.0*k1.Data[j-1]/216.0 + 1408.0*k3.Data[j-1]/2565.0 +
					2197.0*k4.Data[j-1]/4104.0 - k5.Data[j-1]/5.0
				soltemp1 = sol0.GetFromMatrix(j, i-1) + soltemp1
				sol0.SetMatrix(j, i, soltemp1)
			}
			i++
			nreal = i
			continue
		}

		//最大步数强边界
		if i >= n {
			break
		}

		//变步长
		scale := tol * h / (2.0 * math.Abs(errtemp0))
		scale = math.Pow(scale, 0.25)
		if scale < 0.75 {
			h = h / 2.0
		} else if scale > 1.5 {
			h = h * 2.0
		}
	}

	//解矩阵缩减
	sol := ZeroMatrix(fn+1, nreal)
	for j := 0; j < nreal; j++ {
		for k := 0; k < fn+1; k++ {
			sol.SetMatrix(k, j, sol0.GetFromMatrix(k, j))
		}
	}

	err = true
	return sol, err
}
