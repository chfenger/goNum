// InterpSpline12
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-8
版本   : 0.0.0
------------------------------------------------------
    用节点处的一阶导数表示的三次样条插值函数,
    二阶导数边界条件
    n+1个点, n个区间
理论：
    区间[x(i-1), xi]上的三次样条函数表达为：
            (x-xi)^2 * [hi+2(x-x(i-1))]
   Si(x) = -----------------------------y(i-1) +
                       hi^3
            (x-x(i-1))^2 * [hi+2(xi-x)]
           -----------------------------yi +
                       hi^3
            (x-xi)^2 * (x-x(i-1))
           -----------------------m(i-1) +
                     hi^2
            (x-x(i-1))^2 * (x-xi)
           -----------------------mi
                     hi^2

    令 lambdai = h(i+1)/(hi+h(i+1))
       Mi = 1-lambdai = hi/(hi+h(i+1))
                 y(i+1)-yi           yi-y(i-1)
       fi = 3(Mi---------- + lambdai-----------)
                  h(i+1)                hi
       (i = 1,...,n-1)

    则mi可由n+1阶线性方程组求得（利用LEs_Chasing）：
    |2  1                  ||  m0  |   |  f0  |
    |l1 2  M1              ||  m1  |   |  f1  |
    |   l2 2  M2           ||  m2  | = |  f2  |
    |      ........        || ...  |   | ...  |
    |       l(n-1) 2 M(n-1)||m(n-1)|   |f(n-1)|
    |              1 2     ||  mn  |   |  fn  |

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 116-123.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；
            第三列y''i，且y''i只需给出y''0和y''n
输出   :
    B       插值方程系数结果矩阵，从前到后对应从0到3阶，4xn
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import "math"

// InterpSpline12 用节点处的一阶导数表示的三次样条插值函数，二阶导数边界条件
func InterpSpline12(A Matrix) (Matrix, bool) {
	/*
		用节点处的一阶导数表示的三次样条插值函数，二阶导数边界条件
		输入   :
		    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；
		            第三列y'i，且y'i只需给出y'0和y'n
		输出   :
		    B       插值方程系数结果矩阵，从前到后对应从0到3阶，4xn
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	var err bool = false
	n := A.Rows - 1
	sol := ZeroMatrix(4, n)
	BA := ZeroMatrix(n+1, n+1) //对角占优的三对角矩阵
	BB := ZeroMatrix(n+1, 1)   //解向量
	BC := ZeroMatrix(n+1, 1)   //值向量

	//1解插值函数的一阶导数mi
	//1.0.1第一行
	if true { //限制变量使用范围
		h1 := A.GetFromMatrix(1, 0) - A.GetFromMatrix(0, 0)
		y0 := A.GetFromMatrix(0, 1)
		y1 := A.GetFromMatrix(1, 1)

		f0 := 3.0*(y1-y0)/h1 - h1*A.GetFromMatrix(0, 2)/2.0

		BA.SetMatrix(0, 0, 2.0)
		BA.SetMatrix(0, 1, 1.0)
		BC.Data[0] = f0
	}
	//1.0.2其它行
	for i := 1; i < n; i++ {
		yi_1 := A.GetFromMatrix(i-1, 0)
		yi := A.GetFromMatrix(i, 0)
		yi1 := A.GetFromMatrix(i+1, 0)

		hi := A.GetFromMatrix(i, 0) - A.GetFromMatrix(i-1, 0)
		hi1 := A.GetFromMatrix(i+1, 0) - A.GetFromMatrix(i, 0)

		lambdai := hi1 / (hi + hi1)
		Mi := 1.0 - lambdai
		fi := 3.0 * (Mi*(yi1-yi)/hi1 + lambdai*(yi-yi_1)/hi)
		//赋予BA
		BA.SetMatrix(i, i-1, lambdai)
		BA.SetMatrix(i, i, 2.0)
		BA.SetMatrix(i, i+1, Mi)
		BC.Data[i] = fi
	}
	//1.0.3最后一行
	if true { //i=n
		hn := A.GetFromMatrix(n, 0) - A.GetFromMatrix(n-1, 0)
		yn_1 := A.GetFromMatrix(n-1, 1)
		yn := A.GetFromMatrix(n, 1)

		fn := 3.0*(yn-yn_1)/hn + hn*A.GetFromMatrix(n, 2)/2.0

		BA.SetMatrix(n, n-1, 1.0)
		BA.SetMatrix(n, n, 2.0)
		BC.Data[n] = fn
	}
	//1.1求解
	soltemp, errtemp := LEs_Chasing(BA, BC)
	if errtemp != true {
		panic("Error in goNum.InterpSpline12: Solve Error with goNum.LEs_Chasing")
	}
	for i := 0; i < n+1; i++ {
		BB.Data[i] = soltemp.Data[i]
	}

	//2求解Si(x)
	S0 := ZeroMatrix(4, 1)
	S1 := ZeroMatrix(4, 1)
	S2 := ZeroMatrix(4, 1)
	S3 := ZeroMatrix(4, 1)
	for i := 1; i < n+1; i++ {
		xi_1 := A.GetFromMatrix(i-1, 0)
		xi := A.GetFromMatrix(i, 0)
		yi_1 := A.GetFromMatrix(i-1, 1)
		yi := A.GetFromMatrix(i, 1)
		mi_1 := BB.Data[i-1]
		mi := BB.Data[i]
		hi := xi - xi_1

		temp0 := ZeroMatrix(4, 1)
		temp1 := ZeroMatrix(4, 1)
		//2.1 S0
		temp0.Data[2] = 1.0
		temp0.Data[1] = -2.0 * xi
		temp0.Data[0] = xi * xi
		for j := 3; j > 0; j-- {
			temp0.Data[j] = 2.0 * temp0.Data[j-1]
			temp1.Data[j-1] = (hi - 2.0*xi_1) * temp0.Data[j-1]
			S0.Data[j] = (temp0.Data[j] + temp1.Data[j]) * yi_1 / math.Pow(hi, 3.0)
		}
		S0.Data[0] = temp1.Data[0] * yi_1 / math.Pow(hi, 3.0)
		//2.1 S1
		temp0 = ZeroMatrix(4, 1)
		temp1 = ZeroMatrix(4, 1)
		temp0.Data[2] = 1.0
		temp0.Data[1] = -2.0 * xi_1
		temp0.Data[0] = xi_1 * xi_1
		for j := 3; j > 0; j-- {
			temp0.Data[j] = -2.0 * temp0.Data[j-1]
			temp1.Data[j-1] = (hi + 2.0*xi) * temp0.Data[j-1]
			S1.Data[j] = (temp0.Data[j] + temp1.Data[j]) * yi / math.Pow(hi, 3.0)
		}
		S1.Data[0] = temp1.Data[0] * yi / math.Pow(hi, 3.0)
		//2.2 S2
		temp0 = ZeroMatrix(4, 1)
		temp1 = ZeroMatrix(4, 1)
		temp0.Data[2] = 1.0
		temp0.Data[1] = -2.0 * xi
		temp0.Data[0] = xi * xi
		for j := 3; j > 0; j-- {
			temp0.Data[j] = temp0.Data[j-1]
			temp1.Data[j-1] = -1.0 * xi_1 * temp0.Data[j-1]
			S2.Data[j] = (temp0.Data[j] + temp1.Data[j]) * mi_1 / math.Pow(hi, 2.0)
		}
		S2.Data[0] = temp1.Data[0] * mi_1 / math.Pow(hi, 2.0)
		//2.3 S3
		temp0 = ZeroMatrix(4, 1)
		temp1 = ZeroMatrix(4, 1)
		temp0.Data[2] = 1.0
		temp0.Data[1] = -2.0 * xi_1
		temp0.Data[0] = xi_1 * xi_1
		for j := 3; j > 0; j-- {
			temp0.Data[j] = temp0.Data[j-1]
			temp1.Data[j-1] = -1.0 * xi * temp0.Data[j-1]
			S3.Data[j] = (temp0.Data[j] + temp1.Data[j]) * mi / math.Pow(hi, 2.0)
		}
		S3.Data[0] = temp1.Data[0] * mi / math.Pow(hi, 2.0)
		//2.4 Si(x)
		for j := 0; j < 4; j++ {
			sol.SetMatrix(j, i-1, S0.Data[j]+S1.Data[j]+S2.Data[j]+S3.Data[j])
		}
	}

	err = true
	return sol, err
}
