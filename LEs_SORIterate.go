// LEs_SORIterate
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-22
版本   : 0.0.0
------------------------------------------------------
    解n阶线性方程组的SOR(逐次超松弛, successive over
       relaxation)迭代法
理论：
    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 68-72.
    收敛的条件：（B为变化后的系数矩阵）
       1. 系数矩阵A严格对角占优，且0 < omega <= 1，或者
       2. 系数矩阵A对称正定，且0 < omega < 2
------------------------------------------------------
输入   :
    A       系数矩阵
    b       常数值向量
    tol     最大容许误差
    omega   松弛因子，0 < omega < 2, omega = 1: Siedel,
            omega < 1: 低松弛, omega > 1: 超松弛
    n       最大迭代步数
输出   :
    sol     解向量
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import "math"

func LEs_SORIterate(A, b, x0 Matrix, tol, omega float64, n int) ([]float64, bool) {
	/*
		解n阶线性方程组的SOR(逐次超松弛, successive over relaxation)迭代法
		输入   :
		    A       系数矩阵
		    b       常数值向量
		    tol     最大容许误差
		    omega   松弛因子，0 < omega < 2, omega = 1: Siedel,
		            omega < 1: 低松弛, omega > 1: 超松弛
		    n       最大迭代步数
		输出   :
		    sol     解向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	x1 := ZeroMatrix(A.Rows, 1)
	sol := ZeroMatrix(A.Rows, 1)
	var err bool = false

	//求解
	for i := 0; i < n; i++ {
		for i0 := 0; i0 < A.Rows; i0++ {
			sum0 := 0.0
			for j := 0; j < i0; j++ {
				sum0 += A.GetFromMatrix(i0, j) * x1.GetFromMatrix(j, 0)
			}
			sum1 := 0.0
			for j := i0 + 1; j < A.Columns; j++ {
				sum1 += A.GetFromMatrix(i0, j) * x0.GetFromMatrix(j, 0)
			}
			x1.SetMatrix(i0, 0, (1-omega)*x0.GetFromMatrix(i0, 0)+omega*(b.Data[i0]-sum0-sum1)/A.GetFromMatrix(i0, i0))
		}

		//判断收敛
		sol = SubMatrix(x1, x0)
		max, _, _ := Max(sol.Data)
		if math.Abs(max) < tol {
			sol = x1
			err = true
			return sol.Data, err
		}

		//准备下次迭代
		for i0 := 0; i0 < x0.Rows; i0++ {
			x0.Data[i0] = x1.Data[i0]
		}
	}

	return make([]float64, A.Rows), err
}
