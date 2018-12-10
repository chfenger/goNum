// RK22
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-8
版本   : 0.0.0
------------------------------------------------------
    二级二阶Runge-Kutta法求解常微分方程组
理论：


    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 192-199.
------------------------------------------------------
输入   :
    fun     第i个方程(计算变量值向量, i)
    x0      初值向量，(fn+1)x1，一个x，fn个因变量
    xend    终止x
    fn      方程个数
    n       最大迭代步数
输出   :
    B       解向量
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

func RK22(fun func(Matrix, int) float64, x0 Matrix, xend float64, fn, n int) (Matrix, bool) {
	/*
		二级二阶Runge-Kutta法求解常微分方程组
		输入   :
		    fun     第i个方程(计算变量值向量, i)
		    x0      初值向量，(fn+1)x1，一个x，fn个因变量
		    xend    终止x
		    fn      方程个数
		    n       最大迭代步数
		输出   :
		    B       解向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断方程个数是否对应初值个数
	if x0.Rows != fn+1 {
		panic("Error in goNum.RK22: Quantities of x0 and fn are not equal")
	}

	sol := ZeroMatrix(fn+1, n+1)
	var err bool = false
	h := (xend - x0.Data[0]) / float64(n) //步长

	//稳定性条件，建议迭代时即时判断，但此举会拖慢速度
	// if true {
	// 	lambda := ZeroMatrix(fn, 1)
	// 	for j := 0; j < fn; j++ { //微分方程迭代
	// 		lambda.Data[j] = fun(x0, j) / x0.Data[j+1]
	// 	}
	// 	maxl, _, _ := Max(lambda.Data)
	// 	stab := 1.0 + maxl*h + math.Pow(maxl*h, 2.0)/2.0
	// 	if math.Abs(stab) > 1 {
	// 		panic("Error in goNum.RK22: Step length too large or step number little less")
	// 	}
	// }

	//把初值赋给sol
	for i := 0; i < fn+1; i++ {
		sol.SetMatrix(i, 0, x0.Data[i])
	}

	for i := 1; i < n+1; i++ { //最大迭代次数迭代
		temp0 := ZeroMatrix(fn+1, 1)
		//给temp0赋i-1步值，每一步开始
		for j := 0; j < fn+1; j++ {
			temp0.Data[j] = sol.GetFromMatrix(j, i-1)
		}
		k1 := ZeroMatrix(fn, 1)
		k2 := ZeroMatrix(fn, 1)
		//1. k1
		for j := 0; j < fn; j++ { //微分方程迭代
			k1.Data[j] = h * fun(temp0, j)
		}
		//2. k2
		temp0.Data[0] = sol.GetFromMatrix(0, i-1) + 2.0*h/3.0 //xn+2h/3
		for j := 1; j < fn+1; j++ {                           //yn+2k1/3
			temp0.Data[j] = sol.GetFromMatrix(j, i-1) + 2.0*k1.Data[j-1]/3.0
		}
		for j := 0; j < fn; j++ { //微分方程迭代
			k2.Data[j] = h * fun(temp0, j)
		}

		//i步值
		sol.SetMatrix(0, i, sol.GetFromMatrix(0, i-1)+h) //xi
		for j := 1; j < fn+1; j++ {
			temp1 := sol.GetFromMatrix(j, i-1) + (k1.Data[j-1]+3.0*k2.Data[j-1])/4.0
			sol.SetMatrix(j, i, temp1)
		}
	}

	err = true
	return sol, err
}
