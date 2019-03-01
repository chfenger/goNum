// InterpLagrange
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-3
版本   : 0.0.0
------------------------------------------------------
    求解n次拉格朗日Lagrange插值法拟合n+1个数据点
    满阶插值，即阶数为给定点数-1
    内插/外插
理论：
              n       omega0n+1(xq)
    Ln(xq) = Sum(-----------------------)
             k=0  (xq-xk)*omega1n+1(xk)

                       n
      omega0n+1(xq) = Prod(xq-xi)
                      i=0

                        n
      omega1n+1(xk) =  Prod  (xk-xi)
                     i=0,i!=k

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 94-100.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
    xq      插值点
    n       最大插值阶数 1 <= ... <= n
输出   :
    sol     插值结果
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// InterpLagrange 求解n次拉格朗日Lagrange插值法拟合n+1个数据点
func InterpLagrange(A Matrix, xq float64) (float64, bool) {
	/*
		求解n次拉格朗日Lagrange插值法拟合n+1个数据点
		输入   :
		    A       数据点矩阵，(n+1)x2，第一列xi；第二列yi
		    xq      插值点
		    n       最大插值阶数 1 <= ... <= n
		输出   :
		    sol     插值结果
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	var sol float64
	var err bool = false
	n := A.Rows - 1

	//计算系数矩阵
	for k := 0; k <= n; k++ {
		//1. 计算分子
		var temp0 float64 = 1.0
		for i := 0; i <= n; i++ {
			temp0 = temp0 * (xq - A.GetFromMatrix(i, 0))
		}
		temp0 = temp0 / (xq - A.GetFromMatrix(k, 0))
		//2. 计算分母
		var temp1 float64 = 1.0
		for i := 0; i <= n; i++ {
			if i != k {
				temp1 = temp1 * (A.GetFromMatrix(k, 0) - A.GetFromMatrix(i, 0))
			}
		}
		//3. 求和
		sol += temp0 * A.GetFromMatrix(k, 1) / temp1
	}

	err = true
	return sol, err
}
