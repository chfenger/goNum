// PDEDiffHyperbolic1
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-17
版本   : 0.0.0
------------------------------------------------------
    求解双曲型偏微分方程的差分解法（第一种差分格式）
理论：
    对于抛物型偏微分方程：
     d^2u       d^2u
    ------ = A ------ + B
     dt^2       dx^2

    u(x, 0) = phi(x), (du/dt)_(t=0) = psi(x)
    u(0, t) = u1(t), u(L, t) = u2(t)

    0 < x < L, 0 < t < T

    则差分格式为，x分为m等份，t分为n等份
    u_(i,j+1) = lu_(i+1,j) + 2(1-l)u_(i,j) + lu_(i-1,j) -
                u_(i,j-1) + B*ht^2

    初值需要计算第零层和第一层、左右边界
    第零层：u_(i,0) = phi(i*hx), i=1,2,...,m-1
    第一层：u_(i,1) = u_(i,0) + ht*psi(i*hx)
    左边界：u_(0,j) = u1(j*ht)
    右边界：u_(m,j) = u2(j*ht), j=0,1,2,...,n

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 226-228.
------------------------------------------------------
输入   :
    funphi, funpsi, funu1, funu2   边界函数
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

// PDEDiffHyperbolic1 求解双曲型偏微分方程的差分解法（第一种差分格式）
func PDEDiffHyperbolic1(funphi, funpsi, funu1, funu2 func(float64) float64,
	x0 Matrix, A, B float64, m, n int) (Matrix, bool) {
	/*
		求解双曲型偏微分方程的差分解法（第一种差分格式）
		输入   :
		    funphi, funpsi, funu1, funu2   边界函数
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
		panic("Error in goNum.PDEDiffHyperbolic1: Grid numbers error")
	}

	var err bool = false
	sol := ZeroMatrix(m+1, n+1)
	hx := (x0.GetFromMatrix(1, 0) - x0.GetFromMatrix(0, 0)) / float64(m) //x方向步长
	ht := (x0.GetFromMatrix(1, 1) - x0.GetFromMatrix(0, 1)) / float64(n) //t方向步长

	//1. 计算t第零层上的值u_(i,0) i=,1,...,m-1
	for i := 1; i < m; i++ {
		sol.SetMatrix(i, 0, funphi(x0.GetFromMatrix(0, 0)+float64(i)*hx))
	}
	//2. 计算x左右边界上的节点u_(0,j)和u_(m,j) j=0,1,2,...,n
	for j := 0; j < n+1; j++ {
		sol.SetMatrix(0, j, funu1(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //左边界
		sol.SetMatrix(m, j, funu2(x0.GetFromMatrix(0, 1)+float64(j)*ht)) //右边界
	}
	//lambda及稳定性判断
	l := A * ht * ht / (hx * hx)
	if l > 1 {
		panic("Error in goNum.PDEDiffHyperbolic1: lambda greater than one")
	}
	//3. 计算t第一层上的值u_(i,1) i=,1,...,m-1
	for i := 1; i < m; i++ {
		sol.SetMatrix(i, 1, sol.GetFromMatrix(i, 0)+ht*funpsi(x0.GetFromMatrix(0, 0)+float64(i)*hx))
	}
	//4. 2～n层
	for j := 2; j < n+1; j++ {
		for i := 1; i < m; i++ {
			temp0 := l * sol.GetFromMatrix(i+1, j-1)
			temp0 += 2.0 * (1.0 - l) * sol.GetFromMatrix(i, j-1)
			temp0 += l * sol.GetFromMatrix(i-1, j-1)
			temp0 -= sol.GetFromMatrix(i, j-2)
			temp0 += B * ht * ht
			sol.SetMatrix(i, j, temp0)
		}
	}

	err = true
	return sol, err
}
