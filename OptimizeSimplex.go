// OptimizeSimplex
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-25
版本   : 0.0.0
------------------------------------------------------
    Nelder-Mead单纯形法求解多自变量函数极小值
理论：
    对于函数z=f(x0,x1,...,xn)，取三个相异的点构成三角形，并按函数值
    从小到大排序为B、G、W，依下列方法进行操作：
    0. 取BG中点M = (B+G)/2；
    1. 取反射点R = M+(M-W)；
    2. 取延伸点E = R+(R-M)；
    3. 收缩点C = zMin(C1=M+(W-M)/2, C2=M+(M-W)/2)；
    4. 收缩点S = (B+W)/2。
    1~4每一步计算函数值并置换排序BGW
    n个x需要n+1个初始点

       参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 8.2.1
------------------------------------------------------
输入   :
    fun     函数表达式
    x0      初始点，nx(n+1)，第一行x0，第二行x1,...
    tol     控制误差
    Nn      最大迭代步数
输出   :
    sol     解,(n+1)x1
    |xPath  自变量变化历程，二维浮点，可使用Slices2ToMatrix函数转换为Matrix类型
    |fxPath 函数值变化历程，一维浮点，可使用Slices1ToMatrix函数转换为Matrix类型
    |errPath 误差绝对值历程，一维浮点，可使用Slices1ToMatrix函数转换为Matrix类型
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func OptimizeSimplex(fun func(Matrix) float64, x0 Matrix, tol float64, Nn int) (Matrix, bool) {
	/*
		Nelder-Mead单纯形法求解多自变量函数极小值
		输入   :
		    fun     函数表达式
		    x0      初始点，nx(n+1)，第一行x0，第二行x1,...
		    tol     控制误差
		    Nn      最大迭代步数
		输出   :
		    sol     解,(n+1)x1
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断x0大小
	n := x0.Rows           //xi
	if x0.Columns != n+1 { //初始点个数等于自变量个数加一
		panic("Error in goNum.OptimizeSimplex: Initial values error")
	}
	//判断N
	if Nn < 1 {
		panic("Error in goNum.OptimizeSimplex: Iteration number error")
	}

	sol := ZeroMatrix(n+1, 1)
	xPath := make([][]float64, 0)
	fxPath := make([]float64, 0)
	errPath := make([]float64, 0)
	var err bool = false

	//计算f(x)
	for i := 0; i < n+1; i++ {
		sol.Data[i] = fun(NewMatrix(n, 1, x0.ColumnOfMatrix(i)))
	}

	//取最大、最小、次大和次小序号
	_, l0, _ := Min(sol.Data) //最小
	_, h0, _ := Max(sol.Data) //最大
	l1 := h0                  //次小
	h1 := l0                  //次大
	for i := 0; i < n+1; i++ {
		if (i != l0) && (i != h0) && (sol.Data[i] < sol.Data[l1]) {
			l1 = i
		}
		if (i != l0) && (i != h0) && (sol.Data[i] > sol.Data[h1]) {
			h1 = i
		}
	}
	xPath = append(xPath, x0.ColumnOfMatrix(l0))
	fxPath = append(fxPath, sol.Data[l0])
	errPath = append(errPath, math.Abs(sol.Data[h0]-sol.Data[l0]))

	//迭代
	for i := 0; i < Nn; i++ {
		//中点M = (Sum-W)/n
		temp0 := ZeroMatrix(n, 1)
		for j := 0; j < n+1; j++ {
			temp0 = AddMatrix(temp0, NewMatrix(n, 1, x0.ColumnOfMatrix(j)))
		}
		mm := NumProductMatrix(SubMatrix(temp0, NewMatrix(n, 1, x0.ColumnOfMatrix(h0))), 1.0/float64(n))
		//反射点R = 2M-W
		rr := SubMatrix(NumProductMatrix(mm, 2.0), NewMatrix(n, 1, x0.ColumnOfMatrix(h0)))
		fr := fun(rr)
		//判断
		if fr < sol.Data[h1] { //fr<fh1, case1
			if fr > sol.Data[l1] { //R-->W
				for j := 0; j < n; j++ {
					x0.SetMatrix(j, h0, rr.Data[j])
				}
				sol.Data[h0] = fr
			} else { //延伸E
				ee := SubMatrix(NumProductMatrix(rr, 2.0), mm)
				fe := fun(ee)
				if fe < sol.Data[l1] { //E-->W
					for j := 0; j < n; j++ {
						x0.SetMatrix(j, h0, ee.Data[j])
					}
					sol.Data[h0] = fe
				} else { //R-->W
					for j := 0; j < n; j++ {
						x0.SetMatrix(j, h0, rr.Data[j])
					}
					sol.Data[h0] = fr
				}
			}
		} else { //case 2
			if fr < sol.Data[h0] {
				for j := 0; j < n; j++ {
					x0.SetMatrix(j, h0, rr.Data[j])
				}
				sol.Data[h0] = fr
			}
			//C1 = (W+M)/2, C2 = (R+M)/2，默认C=C1
			cc := NumProductMatrix(AddMatrix(NewMatrix(n, 1, x0.ColumnOfMatrix(h0)), mm), 0.5)
			fc := fun(cc)
			c2 := NumProductMatrix(AddMatrix(rr, mm), 0.5)
			fc2 := fun(c2)
			//判断获得C
			if fc > fc2 {
				for j := 0; j < n; j++ {
					cc.Data[j] = c2.Data[j]
				}
				fc = fc2
			}
			if fc < sol.Data[h0] {
				for j := 0; j < n; j++ {
					x0.SetMatrix(j, h0, cc.Data[j])
				}
				sol.Data[h0] = fc
			} else { //xj = (xj+x0)/2
				for j := 0; j < n+1; j++ {
					if j != l0 {
						temp1 := NumProductMatrix(AddMatrix(NewMatrix(n, 1, x0.ColumnOfMatrix(j)),
							NewMatrix(n, 1, x0.ColumnOfMatrix(l0))), 0.5)
						for k := 0; k < n; k++ {
							x0.SetMatrix(k, j, temp1.Data[k])
						}
						sol.Data[j] = fun(temp1)
					}
				}
			}
		}
		//下一步
		_, l0, _ = Min(sol.Data) //最小
		_, h0, _ = Max(sol.Data) //最大
		l1 = h0                  //次小
		h1 = l0                  //次大
		for j := 0; j < n+1; j++ {
			if (j != l0) && (j != h0) && (sol.Data[j] < sol.Data[l1]) {
				l1 = j
			}
			if (j != l0) && (j != h0) && (sol.Data[j] > sol.Data[h1]) {
				h1 = j
			}
		}
		//记录历程
		xPath = append(xPath, x0.ColumnOfMatrix(l0))
		fxPath = append(fxPath, sol.Data[l0])
		errPath = append(errPath, math.Abs(sol.Data[h0]-sol.Data[l0]))
		//判断满足精度否
		if errPath[i+1] < tol {
			//将所有数据赋予sol，前n项为x，最后一项为f(x)
			sol.Data[n] = sol.Data[l0]
			for j := 0; j < n; j++ {
				sol.Data[j] = x0.GetFromMatrix(j, l0)
			}
			err = true
			return sol, err
		}
	}

	return sol, err //,xPath,fxPath,errPath
}
