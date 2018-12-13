// PowInt
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-12
版本   : 0.0.0
------------------------------------------------------
    计算整数或浮点数的整数次幂
理论：

------------------------------------------------------
输入   :
    a, n    a^n
输出   :
    sol     解
------------------------------------------------------
*/

package goNum

import (
	"math"
)

//浮点数的整数次幂
func PowFInt(a float64, n int) float64 {
	/*
	   计算浮点数的整数次幂
	   输入   :
	       a, n    a^n
	   输出   :
	       sol     解
	*/
	if n < 0 {
		panic("Error in goNum.PowFInt: n less than zero")
	} else if n == 0 {
		return 1.0
	} else if n == 1 {
		return a
	}
	sol := a
	for i := 2; i < n+1; i++ {
		sol = sol * a
	}
	return sol
}

//整数的浮点数次幂
func PowIF(a int, n float64) float64 {
	/*
	   计算整数的浮点数次幂
	   输入   :
	       a, n    a^n
	   输出   :
	       sol     解
	*/
	return math.Pow(float64(a), n)
}

//整数的整数次幂，输出整数
func PowIInt(a, n int) int {
	/*
	   计算整数的整数次幂，输出整数
	   输入   :
	       a, n    a^n
	   输出   :
	       sol     解
	*/
	if n < 0 {
		panic("Error in goNum.PowIInt: n less than zero")
	} else if n == 0 {
		return 1
	} else if n == 1 {
		return a
	}
	sol := a
	for i := 2; i < n+1; i++ {
		sol = sol * a
	}
	return sol
}

//整数的整数次幂，输出浮点
func PowIIntF(a, n int) float64 {
	/*
	   计算整数的整数次幂，输出浮点
	   输入   :
	       a, n    a^n
	   输出   :
	       sol     解
	*/
	if n < 0 {
		panic("Error in goNum.PowIInt: n less than zero")
	} else if n == 0 {
		return 1.0
	} else if n == 1 {
		return float64(a)
	}
	sol := a
	for i := 2; i < n+1; i++ {
		sol = sol * a
	}
	return float64(sol)
}
