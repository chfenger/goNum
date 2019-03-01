// PDEDiffParabolicE
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-14
版本   : 0.0.0
------------------------------------------------------
    求解抛物型偏微分方程的差分解法（显式）
理论：
    对于抛物型偏微分方程：
     du       d^2u
    ---- = A ------ + B
     dt       dx^2

    u(x, 0) = p(x)
    u(0, t) = u1(t), u(L, t) = u2(t)

    0 < x < L, 0 < t < T

    则古典显式差分格式为，x分为m等份，t分为n等份
    u_(i,j+1) = lu_(i-1,j) + (1-2l)u_(i,j) + lu_(i+1,j) + B*tau
         A*tau
    l = -------
          h^2
    u_(i, 0) = p(ih), i=1,2,..,m-1
    u_(0, j) = u1(j*tau), u_(m, j) = u2(j, tau), j=0,1,...,n

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 214-215.
------------------------------------------------------
输入   :
    funp, funu1, funu2   边界函数
    x0      求解范围，2x2
    A, B    常系数
    m, n    网格数量
输出   :
    sol     解矩阵
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

// PDEDiffParabolicE 求解抛物型偏微分方程的差分解法（显式）
func PDEDiffParabolicE(funp, funu1, funu2 func(float64) float64, x0 Matrix,
	A, B float64, m, n int) (Matrix, bool) {
	/*
		求解抛物型偏微分方程的差分解法（显式）
		输入   :
		    funp, funu1, funu2   边界函数
		    x0      求解范围，2x2
		    A, B    常系数
		    m, n    网格数量
		输出   :
		    sol     解矩阵
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断网格数量
	if (m < 1) || (n < 1) {
		panic("Error in goNum.PDEDiffParabolicE: Grid numbers error")
	}

	var err bool = false
	sol := ZeroMatrix(m+1, n+1)
	hx := (x0.GetFromMatrix(1, 0) - x0.GetFromMatrix(0, 0)) / float64(m) //x方向步长
	ht := (x0.GetFromMatrix(1, 1) - x0.GetFromMatrix(0, 1)) / float64(n) //t方向步长

	//1. 计算t第零层上的值u_(i,0) i=0,1,...,m
	for i := 0; i < m+1; i++ {
		sol.SetMatrix(i, 0, funp(x0.GetFromMatrix(0, 0)+float64(i)*hx))
	}
	//2. 计算左右边界上的节点u_(0,j)和u_(m,j) j=1,2,...,n
	for j := 1; j < n+1; j++ {
		sol.SetMatrix(0, j, funu1(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //左边界
		sol.SetMatrix(m, j, funu2(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //右边界
	}
	//内部节点循环求解
	l := A * ht / (hx * hx)
	//稳定性判断
	if (l <= 0) || (l > 0.5) {
		panic("Error in goNum.PDEDiffParabolicS: lambda less than or equal to zero, or greater than 0.5")
	}
	for j := 1; j < n+1; j++ { //层循环, ti
		for i := 1; i < m; i++ { //列循环, xi
			uij := l * sol.GetFromMatrix(i-1, j-1)
			uij += (1 - 2.0*l) * sol.GetFromMatrix(i, j-1)
			uij += l * sol.GetFromMatrix(i+1, j-1)
			sol.SetMatrix(i, j, uij+B*ht)
		}
	}

	err = true
	return sol, err
}
