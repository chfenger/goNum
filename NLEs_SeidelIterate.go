// NLEs_SeidelIterate
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-20
版本   : 0.0.0
------------------------------------------------------
    多元非线性方程组Seidel迭代
理论：
    Pk = x0
    Fk = [f1, f2,..., fn]'

         |df1/dx1 df1/dx2 ... df1/dxn|
         |df2/dx1 df2/dx2 ... df2/dxn|
    Jk = |...     ...     ... ...    |
         |dfn/dx1 dfn/dx2 ... dfn/dxn|

    Jk*dPk = -Fk
    P_(k+1) = Pk+dPk

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 3.7
------------------------------------------------------
输入   :
    funs    方程组，nx1
    J       Joccobi矩阵，nxn
    x0      初值x
    tol     控制误差
    n       最大迭代次数
输出   :
    sol     解，nx1
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// NLEs_SeidelIterate 多元非线性方程组Seidel迭代
func NLEs_SeidelIterate(funs, J func(Matrix) Matrix, x0 Matrix,
	tol float64, n int) (Matrix, bool) {
	/*
		多元非线性方程组Seidel迭代
		输入   :
		    funs    方程组，nx1
		    J       Joccobi矩阵，nxn
		    x0      初值x
		    tol     控制误差
		    n       最大迭代次数
		输出   :
		    sol     解，nx1
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断x维数
	if x0.Columns != 1 {
		panic("Error in goNum.NLEs_SeidelIterate: x0 is not a vector")
	}

	sol := ZeroMatrix(x0.Rows, 1)  //解向量
	xold := ZeroMatrix(x0.Rows, 1) //Pk
	var err bool = false

	//将x0赋予xold
	for i := 0; i < x0.Rows; i++ {
		xold.Data[i] = x0.Data[i]
		sol.Data[i] = x0.Data[i]
	}
	//循环迭代
	y := NumProductMatrix(funs(xold), -1.0)
	for i := 0; i < n; i++ {
		ja := J(xold)
		dx, dxerr := LEs_ECPE(Matrix2ToSlices(ja), y.Data)
		if dxerr != true {
			panic("Error in goNum.NLEs_SeidelIterate: Solve error")
		}
		//求解新值
		for i := 0; i < x0.Rows; i++ {
			sol.Data[i] = xold.Data[i] + dx[i]
			xold.Data[i] = sol.Data[i]
		}
		y = NumProductMatrix(funs(xold), -1.0)
		//判断误差
		maxy, _, _ := MaxAbs(y.Data)
		if math.Abs(maxy) < tol {
			err = true
			return sol, err
		}
	}

	return sol, err
}
