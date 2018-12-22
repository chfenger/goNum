// Norm
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-21
版本   : 0.0.0
------------------------------------------------------
    求向量p范数
理论：

------------------------------------------------------
输入   :
    A       向量,nx1
    p       指定范数
输出   :
    sol     范数值
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
注释 p :
    1       1 modulus
    2       2 modulus
    p       p modulus
    -1      infinite modulus
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func Norm(A Matrix, p float64) (float64, bool) {
	/*
		求向量p范数
		输入   :
		    A       向量,nx1
		    p       指定范数
		输出   :
		    sol     范数值
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//A的维数
	if A.Columns != 1 {
		panic("Error in goNum.Norm: A is not a vector")
	}
	//判断p的值
	if (p < (-1.0)) || ((p > (-1.0)) && (p <= 0.0)) {
		panic("Error in goNum.Norm: p is wrong")
	}

	var sol float64
	var err bool = false
	switch {
	case p == 1.0: //1范数
		for i := 0; i < A.Rows; i++ {
			sol += math.Abs(A.Data[i])
		}
	case p == 2.0: //2范数
		for i := 0; i < A.Rows; i++ {
			sol += A.Data[i] * A.Data[i]
		}
		sol = math.Sqrt(sol)
	case p == -1.0: //无穷范数
		sol, _, _ = MaxAbs(A.Data)
		sol = math.Abs(sol)
	default: //p范数
		for i := 0; i < A.Rows; i++ {
			sol += math.Pow(A.Data[i], p)
		}
		sol = math.Pow(sol, 1.0/p)
	}

	err = true
	return sol, err
}
