// IntegralRumberg_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-12
版本   : 0.0.0
------------------------------------------------------
    Rumberg(龙贝格)求积分公式
理论：
    对于积分
    b
    |f(x)dx
    a

          b-a
    T1 = -----(f(a)+f(b))
           2

           1       b-a N              b-a
    T2N = ---TN + -----Sum f(a+(2j-1)------)
           2       2N  j=1             2N
    N=2^(k-1), k=1,2,3,...

                1            4T2N-TN
    SN = T2N + ---(T2N-TN) = --------
                3              4-1

                1              4^2S2N-SN
    CN = S2N + ----(S2N-SN) = -----------
                15               4^2-1

                1              4^3C2N-CN
    RN = C2N + ----(C2N-CN) = -----------
                63              4^3-1

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 162-164.
------------------------------------------------------
输入   :
    fun     被积分函数
    a, b    积分范围
    tol     控制误差
    Nn      最大循环步数
输出   :
    sol     解
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"math"
	"testing"

	"github.com/chfenger/goNum"
)

// IntegralRumberg Rumberg(龙贝格)求积分公式
func IntegralRumberg(fun func(float64) float64, a, b, tol float64, Nn int) (float64, bool) {
	/*
		Rumberg(龙贝格)求积分公式
		输入   :
		    fun     被积分函数
		    a, b    积分范围
		    tol     控制误差
		    Nn      最大循环步数
		输出   :
		    sol     解
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/

	T := make([]float64, 0) //梯形序列
	S := make([]float64, 0) //辛浦生序列
	C := make([]float64, 0) //柯特斯序列
	R := make([]float64, 0) //龙贝格序列
	//第一步
	temp0 := (b - a) * (fun(a) + fun(b)) / 2.0
	T = append(T, temp0) //T[0]=T1
	//第二步, k=1
	temp0 = 0.0
	for j := 1; j < goNum.PowIInt(2, 0)+1; j++ {
		temp0 += fun(a + (2.0*float64(j)-1.0)*(b-a)/2.0)
	}
	temp0 = T[0]/2.0 + temp0*(b-a)/2.0
	T = append(T, temp0) //T[1]=T2
	temp1 := T[1] + (T[1]-T[0])/3.0
	S = append(S, temp1) //S[0]=S1
	//第三步, k=2
	temp0 = 0.0
	for j := 1; j < goNum.PowIInt(2, 1)+1; j++ {
		temp0 += fun(a + (2.0*float64(j)-1.0)*(b-a)/(2.0*2.0))
	}
	temp0 = T[1]/2.0 + temp0*(b-a)/(2.0*2.0)
	T = append(T, temp0) //T[2]=T4
	temp1 = T[2] + (T[2]-T[1])/3.0
	S = append(S, temp1) //S[1]=S2
	temp2 := S[1] + (S[1]-S[0])/15.0
	C = append(C, temp2) //C[0]=C1
	//第四步, k=3
	temp0 = 0.0
	for j := 1; j < goNum.PowIInt(2, 2)+1; j++ {
		temp0 += fun(a + (2.0*float64(j)-1.0)*(b-a)/(2.0*4.0))
	}
	temp0 = T[2]/2.0 + temp0*(b-a)/(2.0*4.0)
	T = append(T, temp0) //T[3]=T8
	temp1 = T[3] + (T[3]-T[2])/3.0
	S = append(S, temp1) //S[2]=S4
	temp2 = S[2] + (S[2]-S[1])/15.0
	C = append(C, temp2) //C[1]=C2
	temp3 := C[1] + (C[1]-C[0])/63.0
	R = append(R, temp3) //R[0]=R1
	//进入Rumberg循环
	for i := 1; i < Nn; i++ {
		temp0 = 0.0
		for j := 1; j < goNum.PowIInt(2, i+2)+1; j++ {
			temp0 += fun(a + (2.0*float64(j)-1.0)*(b-a)/(2.0*goNum.PowIIntF(2, i+2)))
		}
		temp0 = T[i+2]/2.0 + temp0*(b-a)/(2.0*goNum.PowIIntF(2, i+2))
		T = append(T, temp0) //T[i+3]
		temp1 = T[i+3] + (T[i+3]-T[i+2])/3.0
		S = append(S, temp1) //S[i+2]
		temp2 = S[i+2] + (S[i+2]-S[i+1])/15.0
		C = append(C, temp2) //C[i+1]
		temp3 = C[i+1] + (C[i+1]-C[i])/63.0
		R = append(R, temp3) //R[i]

		if math.Abs(R[i]-R[i-1]) < tol {
			return R[i], true
		}
	}
	return 0.0, false
}

func fun37(x float64) float64 {
	return 4.0 / (1 + x*x)
}

func BenchmarkIntegralRumberg(b *testing.B) {
	for i := 0; i < b.N; i++ {
		goNum.IntegralRumberg(fun37, 0.0, 1.0, 1e-6, 1e3) //3.141592653...
	}
}
