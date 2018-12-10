// InterpNewton
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-6
版本   : 0.0.0
------------------------------------------------------
    计算x点n次Newton插值结果，拟合n+1个数据点
    满阶插值，即阶数为给定点数-1
理论：
    f(x) = f(x0) + f[x, x0](x-x0)
    f[x, x0] = f[x0, x1] + f[x, x0, x1](x-x1)
    ...
    f(x) = f(x0) + f[x0, x0](x-x0) +
           f[x0, x1, x2](x-x0)(x-x1) +
           ... +
           f[x0, x1, ..., xn](x-x0)(x-x1)...(x-x_(n-1))
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 101-105.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
    xq      插值点, xq!=xi
输出   :
    sol     xq点插值结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

//求差商
func diffq_InterpNewton(A Matrix, k int) float64 {
	var sol float64
	for j := 0; j <= k; j++ {
		xj := A.GetFromMatrix(j, 0)
		//为保证理论可读性，并不采取调用omega1_InterpLagrangeFunc函数的方式
		var temp0 float64 = 1.0
		for i := 0; i <= k; i++ {
			if i != j {
				temp0 = temp0 * (xj - A.GetFromMatrix(i, 0))
			}
		}
		sol += A.GetFromMatrix(j, 1) / temp0
	}
	return sol
}

func InterpNewton(A Matrix, xq float64) (float64, bool) {
	/*
		计算x点n次Newton插值结果，拟合n+1个数据点
		输入   :
		    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
		    xq      插值点, xq!=xi
		输出   :
		    sol     xq点插值结果
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	//判断xq是否等于xi
	for i := 0; i < A.Rows; i++ {
		if math.Abs(xq-A.GetFromMatrix(i, 0)) < 1e-3 {
			panic("Error in goNum.InterpNewton: xq equals about xi")
		}
	}

	var sol float64
	var err bool = false
	n := A.Rows - 1
	BA := ZeroMatrix(n+1, 1)

	//开始计算
	BA.SetMatrix(0, 0, A.GetFromMatrix(0, 1)) //f(x0)
	sol = BA.GetFromMatrix(0, 0)
	for k := 1; k < n+1; k++ {
		//求差商
		BA.SetMatrix(k, 0, diffq_InterpNewton(A, k))
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
