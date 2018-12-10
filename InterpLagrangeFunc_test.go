// InterpLagrangeFunc_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-4
版本   : 0.0.0
------------------------------------------------------
    求解n次拉格朗日Lagrange插值方程系数，拟合n+1个数据点
    满阶插值，即阶数为给定点数-1（因不能确定非满阶时所选取的
      插值点是否合理）
理论：
             n       omega0n+1(x)
    Ln(x) = Sum(----------------------)
            k=0  (x-xk)*omega1n+1(xk)

                      n
      omega0n+1(x) = Prod(x-xi)
                     i=0

                        n
      omega1n+1(xk) =  Prod  (xk-xi)
                     i=0,i!=k

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 94-100.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
输出   :
    B       插值系数矩阵，(n+1)x1，0～n
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

//计算分母
func omega1_InterpLagrangeFunc(A goNum.Matrix, k, n int) float64 {
	var sol float64 = 1.0
	for i := 0; i <= n; i++ {
		if i != k {
			sol = sol * (A.GetFromMatrix(k, 0) - A.GetFromMatrix(i, 0))
		}
	}
	return sol
}

//计算分子并由高阶到低阶排序
func omega0_InterpLagrangeFunc(A goNum.Matrix, k, n int) goNum.Matrix {
	B := goNum.ZeroMatrix(n+1, 1)
	//第零阶 x-x0
	if k == 0 { //如果k==0，则从x1循环
		B.SetMatrix(0, 0, -1.0*A.GetFromMatrix(1, 0))
		B.SetMatrix(1, 0, 1.0)
	}
	if k > 0 { //如果k>0，则从x0循环
		B.SetMatrix(0, 0, -1.0*A.GetFromMatrix(0, 0))
		B.SetMatrix(1, 0, 1.0)
	}
	//其他i!=k阶
	for i := 1; i <= n; i++ {
		if (i != k) && ((k > 0) || ((k == 0) && (i > 1))) {
			if k < i {
				CA := goNum.ZeroMatrix(i+1, 1) //实际i+1行
				CB := goNum.ZeroMatrix(i+1, 1) //实际i行
				//先用x乘以之前每一项，相当于给每一项提升一阶,i+1
				for ii := 1; ii < i+1; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CA.Data[ii] = B.Data[ii-1]
				}
				//再用-xi乘以B的每一有效项,i
				for ii := 0; ii < i; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CB.Data[ii] = -1.0 * A.GetFromMatrix(i, 0) * B.Data[ii]
				}
				//同阶相加赋予B
				for ii := 0; ii < i+1; ii++ {
					B.Data[ii] = CA.Data[ii] + CB.Data[ii]
				}
			} else { // k > i
				CA := goNum.ZeroMatrix(i+2, 1) //实际i+2行
				CB := goNum.ZeroMatrix(i+2, 1) //实际i+1行
				//先用x乘以之前每一项，相当于给每一项提升一阶,i+1
				for ii := 1; ii < i+2; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CA.Data[ii] = B.Data[ii-1]
				}
				//再用-xi乘以B的每一有效项,i+1
				for ii := 0; ii < i+1; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CB.Data[ii] = -1.0 * A.GetFromMatrix(i, 0) * B.Data[ii]
				}
				//同阶相加赋予B
				for ii := 0; ii < i+2; ii++ {
					B.Data[ii] = CA.Data[ii] + CB.Data[ii]
				}
			}

		}
	}
	return B
}

func InterpLagrangeFunc(A goNum.Matrix) (goNum.Matrix, bool) {
	/*
		求解n次拉格朗日Lagrange插值方程系数，拟合n+1个数据点
		输入   :
		    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
		输出   :
		    B       插值系数矩阵，(n+1)x1
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	//阶数
	n := A.Rows - 1

	B := goNum.ZeroMatrix(n+1, 1) //最终系数矩阵
	var err bool = false

	//计算系数矩阵
	for k := 0; k <= n; k++ {
		//1. 计算分母和系数乘积
		temp0 := A.GetFromMatrix(k, 1) / omega1_InterpLagrangeFunc(A, k, n)
		//2. 计算分子并乘以上一步的结果，由高阶到低阶排序
		temp1 := omega0_InterpLagrangeFunc(A, k, n)
		for i := 0; i < temp1.Rows; i++ {
			temp1.Data[i] = temp1.Data[i] * temp0
		}
		//3. 累加
		for i := 0; i < B.Rows; i++ {
			B.Data[i] += temp1.Data[i]
		}
	}

	err = true
	return B, err
}

func BenchmarkInterpLagrangeFunc(b *testing.B) {
	A23 := goNum.NewMatrix(3, 2, []float64{0.0, 1.0,
		1.0, 2.0,
		2.0, 3.0})
	for i := 0; i < b.N; i++ {
		InterpLagrangeFunc(A23)
	}
}
