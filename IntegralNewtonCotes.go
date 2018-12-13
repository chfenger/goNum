// IntegralNewtonCotes
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-11
版本   : 0.0.0
------------------------------------------------------
    1-8级Newton-Cotes求积分公式
理论：
    对于积分
    b           n
    |f(x)dx ~= Sum Ak*f(xk)
    a          k=0

               (n)
    Ak = (b-a)C
               k

    (n)   (-1)^(n-k)  n
   C   = ------------ |t(t-1)(t-2)...(t-(k-1))(t-(k+1))...(t-n)dt
    k     k!(n-k)!n   0

    特别的，n=1为梯形公式；
           n=2为Simpson（辛浦生）公式；
           n=4为Cotes（科特斯）公式

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 145-153.
------------------------------------------------------
输入   :
    fun     被积分函数
    a, b    积分范围
    n       Newton-Cotes公式级数
输出   :
    sol     解
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
注意   ：
    由于误差得不到有效控制，稳定性无法保证，故而并不是n值越
    大越好，实际应用中很少使用n值较大的Newton-Cotes公式
------------------------------------------------------
*/

package goNum

func IntegralNewtonCotes(fun func(float64) float64, a, b float64, n int) (float64, bool) {
	/*
		1-8级Newton-Cotes求积分公式
		输入   :
		    fun     被积分函数
		    a, b    积分范围
		    n       Newton-Cotes公式级数
		输出   :
		    sol     解
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断n
	if (n < 1) || (n > 8) {
		panic("Error in goNum.IntegralNewtonCotes: n is not correct")
	}
	//判断a, b
	if a == b {
		return 0.0, true
	}

	var sol float64
	var err bool = false

	//计算xi
	xi := ZeroMatrix(n+1, 1)
	for i := 0; i < n+1; i++ {
		xi.Data[i] = a + (b-a)*float64(i)/float64(n)
	}

	//系数切片
	coeff := [][]float64{
		{1.0, 1.0},
		{1.0, 4.0, 1.0},
		{1.0, 3.0, 3.0, 1.0},
		{7.0, 32.0, 12.0, 32.0, 7.0},
		{19.0, 75.0, 50.0, 50.0, 75.0, 19.0},
		{41.0, 216.0, 27.0, 272.0, 27.0, 216.0, 41.0},
		{751.0, 3577.0, 1323.0, 2989.0, 2989.0, 1323.0, 3577.0, 751.0},
		{989.0, 5888.0, -928.0, 10496.0, -4540.0, 10496.0, -928.0, 5888.0, 989.0},
	}

	//计算积分值
	switch n {
	case 1:
		for i := 0; i < n+1; i++ {
			sol += coeff[0][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 2.0
		return sol, true
	case 2:
		for i := 0; i < n+1; i++ {
			sol += coeff[1][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 6.0
		return sol, true
	case 3:
		for i := 0; i < n+1; i++ {
			sol += coeff[2][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 8.0
		return sol, true
	case 4:
		for i := 0; i < n+1; i++ {
			sol += coeff[3][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 90.0
		return sol, true
	case 5:
		for i := 0; i < n+1; i++ {
			sol += coeff[4][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 288.0
		return sol, true
	case 6:
		for i := 0; i < n+1; i++ {
			sol += coeff[5][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 840.0
		return sol, true
	case 7:
		for i := 0; i < n+1; i++ {
			sol += coeff[6][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 17280.0
		return sol, true
	case 8:
		for i := 0; i < n+1; i++ {
			sol += coeff[7][i] * fun(xi.Data[i])
		}
		sol = sol * (b - a) / 28350.0
		return sol, true
	default:
		return 0.0, err
	}

	return sol, err
}
