// InterpNewtonForward
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-6
版本   : 0.0.0
------------------------------------------------------
    计算x点n次Newton前向插值结果，拟合n+1个等距数据点
    Newton前向等距节点插值，满阶插值，即阶数为给定点数-1
理论：
                   ^y0           ^2y0
    f(x) = f(x0) + ---(x-x0)/h + -----(x-x0)(x-x1) +
                    h            2!h^2
    ... +
     ^ny0
    -------(x-x0)(x-x1)...(x-x_(n-1))
     n!h^n
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 107-110.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x2，第一列xi等距分布；第二列yi
    xq      插值点, xq!=xi
输出   :
    sol     xq点插值结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import "math"

//k阶差分
func difff_InterpNewtonForward(A Matrix, k int) float64 {
	sol := A.GetFromMatrix(k, 1) //yk
	for s := 1; s <= k; s++ {
		sol += math.Pow(-1.0, float64(s)) * float64(Cnm(k, s)) * A.GetFromMatrix(k-s, 1)
	}
	return sol
}

func InterpNewtonForward(A Matrix, xq float64) (float64, bool) {
	/*
		计算x点n次Newton前向插值结果，拟合n+1个等距数据点
		输入   :
		    A       数据点矩阵，(n+1)x2，第一列xi等距分布；第二列yi
		    xq      插值点, xq!=xi
		输出   :
		    sol     xq点插值结果
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断xq是否等于xi
	for i := 0; i < A.Rows; i++ {
		if math.Abs(xq-A.GetFromMatrix(i, 0)) < 1e-3 {
			return A.GetFromMatrix(i, 1), true
		}
	}
	//判断xi是否等距节点
	for i := 0; i < A.Rows; i++ {
		x0 := A.GetFromMatrix(0, 0)
		if math.Abs(A.GetFromMatrix(i, 0)-float64(i)*x0) < 1e-3 {
			panic("Error in goNum.InterpNewtonForward: xi is not in equidistance")
		}
	}

	var sol float64
	var err bool = false
	n := A.Rows - 1
	h := A.GetFromMatrix(n, 0) - A.GetFromMatrix(n-1, 0)
	BA := ZeroMatrix(n+1, 1)

	//计算
	BA.SetMatrix(0, 0, A.GetFromMatrix(0, 1)) //f(x0)
	sol = BA.GetFromMatrix(0, 0)
	for k := 1; k < n+1; k++ {
		//求差分
		BA.SetMatrix(k, 0, difff_InterpNewtonForward(A, k))
		//乘系数1/(k!h^k)
		BA.Data[k] = BA.Data[k] / (float64(Factorial(k)) * math.Pow(h, float64(k)))
		//求乘积
		for j := 0; j < k; j++ {
			BA.Data[k] = BA.Data[k] * (xq - A.GetFromMatrix(j, 0))
		}
		//累加
		sol += BA.Data[k]
	}

	err = true
	return sol, err
}
