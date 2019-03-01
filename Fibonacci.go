// Fibonacci
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-24
版本   : 0.0.0
------------------------------------------------------
    求Fibonacci数列
理论：

------------------------------------------------------
输入   :
    n       Fibonacci数列参数
输出   :
    sol     解
------------------------------------------------------
*/

package goNum

// Fibonacci 求Fibonacci数列
func Fibonacci(n int) int {
	/*
	   求Fibonacci数列
	   输入   :
	       n       Fibonacci数列参数
	   输出   :
	       sol     解
	*/
	//判断n
	F := make([]int, n+1)
	if n == 0 {
		F[0] = 0
		return 0
	} else if n == 1 {
		F[1] = 1
		return 1
	}

	F[0] = 0
	F[1] = 1
	for i := 2; i < n+1; i++ {
		F[i] = F[i-1] + F[i-2]
	}

	return F[n]
}
