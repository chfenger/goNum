// InterpHermiteFunc_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-8
版本   : 0.0.0
------------------------------------------------------
    计算不高于2n+1次Hermite插值方程，拟合n+1个函数值数据
        点和对应的n+1个一阶导数点
    满阶插值，即阶数不高于2n+1
理论：
                n
    H2n+1(x) = Sum (alphaj(x)*yj+betaj(x)*mj)
               j=0

    yj, mj分别为函数值和一阶导数值

                            n         1
    alphaj(x) = (1-2(x-xj)*Sum     -------)lj^2(x)
                          k=0,k!=j  xj-xk

    betaj(x) = (x-xj)lj^2(x)

              (x-x0)(x-x1)...(x-xn)
    lj(x) = --------------------------, (被减数不含xj项)
             (xj-x0)(xj-x1)...(xj-xn)

    参考 李信真, 车刚明, 欧阳洁, 等. 计算方法. 西北工业大学
       出版社, 2000, pp 111-113.
------------------------------------------------------
输入   :
    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；第三列y'i
输出   :
    B       插值方程系数结果，从前到后对应从0到2n+1阶，(2n+2)x1
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

//求解lj(x)
func ljx_InterpHermiteFunc(A goNum.Matrix, j int) goNum.Matrix {
	Bljx := goNum.ZeroMatrix(A.Rows, 1) //出去j，加常数项总共n+1
	xj := A.GetFromMatrix(j, 0)
	// j!=0
	temp0 := (xj - A.GetFromMatrix(0, 0)) //x0
	Bljx.Data[0] = -1.0 * A.GetFromMatrix(0, 0)
	Bljx.Data[1] = 1.0
	// j==0
	if j == 0 {
		temp0 = (xj - A.GetFromMatrix(1, 0)) //x1
		Bljx.Data[0] = -1.0 * A.GetFromMatrix(1, 0)
	}
	//其它
	for i := 1; i < A.Rows; i++ {
		if (i != j) && (((j == 0) && (i > 1)) || (j > 0)) {
			if i < j {
				CA := goNum.ZeroMatrix(i+2, 1) //实际i+2行
				CB := goNum.ZeroMatrix(i+2, 1) //实际i+1行
				//先用x乘以之前每一项，相当于给每一项提升一阶,i+1
				//再用-xi乘以B的每一有效项,i+1
				CB.Data[0] = -1.0 * A.GetFromMatrix(i, 0) * Bljx.Data[0]
				for ii := 1; ii < i+1; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CA.Data[ii] = Bljx.Data[ii-1]
					CB.Data[ii] = -1.0 * A.GetFromMatrix(i, 0) * Bljx.Data[ii]
				}
				CA.Data[i+1] = Bljx.Data[i]
				//同阶相加赋予B
				for ii := 0; ii < i+2; ii++ {
					Bljx.Data[ii] = CA.Data[ii] + CB.Data[ii]
				}
			} else { //i>j
				CA := goNum.ZeroMatrix(i+1, 1) //实际i+1行
				CB := goNum.ZeroMatrix(i+1, 1) //实际i行
				//先用x乘以之前每一项，相当于给每一项提升一阶,i+1
				//再用-xi乘以B的每一有效项,i
				CB.Data[0] = -1.0 * A.GetFromMatrix(i, 0) * Bljx.Data[0]
				for ii := 1; ii < i; ii++ {
					//单列可以这样，否则只能用SetMatrix和GetFromMatrix方法
					CA.Data[ii] = Bljx.Data[ii-1]
					CB.Data[ii] = -1.0 * A.GetFromMatrix(i, 0) * Bljx.Data[ii]
				}
				CA.Data[i] = Bljx.Data[i-1]
				//同阶相加赋予B
				for ii := 0; ii < i+1; ii++ {
					Bljx.Data[ii] = CA.Data[ii] + CB.Data[ii]
				}
			}
			temp0 = temp0 * (xj - A.GetFromMatrix(i, 0))
		}
	}
	for i := 0; i < Bljx.Rows; i++ {
		Bljx.Data[i] = Bljx.Data[i] / temp0
	}
	return Bljx
}

//求解ljx^2
func ljx2_InterpHermiteFunc(A goNum.Matrix, j int) goNum.Matrix {
	ljx := ljx_InterpHermiteFunc(A, j)      //n+1 rows
	BA := goNum.ZeroMatrix(2*ljx.Rows-1, 1) //2n+1 rows
	for i := ljx.Rows - 1; i >= 0; i-- {
		for j := ljx.Rows - 1; j >= 0; j-- {
			BA.Data[i+j] = BA.Data[i+j] + ljx.Data[i]*ljx.Data[j]
		}
	}
	return BA
}

//求解alphajx和betajx，合并是为了减少对ljx2_InterpHermiteFunc的调用
func alphabetajx_InterpHermiteFunc(A goNum.Matrix, j int) (goNum.Matrix, goNum.Matrix) {
	var temp0 float64
	ljx2 := ljx2_InterpHermiteFunc(A, j)        //2n+1 rows
	alphajx := goNum.ZeroMatrix(ljx2.Rows+1, 1) //2n+2 rows
	betajx := goNum.ZeroMatrix(ljx2.Rows+1, 1)  //2n+2 rows
	xj := A.GetFromMatrix(j, 0)

	//计算alphajx中的求和
	for k := 0; k < A.Rows; k++ {
		if k != j {
			temp0 += 1.0 / (xj - A.GetFromMatrix(k, 0))
		}
	}
	temp0 = 2.0 * temp0

	//alphajx = (temp0*xj+1 - temp0*x) * ljx2
	//betajx = (x-xj) * ljx2
	//2n+1阶,2n+2行
	alphajx.Data[alphajx.Rows-1] = -1.0 * temp0 * ljx2.Data[ljx2.Rows-1]
	betajx.Data[betajx.Rows-1] = ljx2.Data[ljx2.Rows-1]
	//其它非零阶, alphajx.Rows-2 == betajx.Rows-2 == ljx2.Rows-1
	for i := alphajx.Rows - 2; i > 0; i-- {
		alphajx.Data[i] = (temp0*xj+1.0)*ljx2.Data[i] - 1.0*temp0*ljx2.Data[i-1]
		betajx.Data[i] = -1.0*xj*ljx2.Data[i] + ljx2.Data[i-1]
	}
	//零阶
	alphajx.Data[0] = (temp0*xj + 1.0) * ljx2.Data[0]
	betajx.Data[0] = -1.0 * xj * ljx2.Data[0]

	return alphajx, betajx
}

// InterpHermiteFunc 计算不高于2n+1次Hermite插值方程，拟合n+1个函数值数据点和对应的n+1个一阶导数点
func InterpHermiteFunc(A goNum.Matrix) (goNum.Matrix, bool) {
	/*
		计算不高于2n+1次Hermite插值方程，拟合n+1个函数值数据点和对应的n+1个一阶导数点
		输入   :
		    A       数据点矩阵，(n+1)x3，第一列xi；第二列yi；第三列y'i
		输出   :
		    B       插值方程系数结果，从前到后对应从0到2n+1阶，(2n+2)x1
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//判断A列数是否为3
	if A.Columns != 3 {
		panic("Error in goNum.InterpHermite: give me xi, yi and y'i")
	}

	var err bool = false
	n := A.Rows - 1
	BA := goNum.ZeroMatrix(2*n+2, 1)

	for j := 0; j <= n; j++ {
		alphajx, betajx := alphabetajx_InterpHermiteFunc(A, j)
		for i := 0; i < alphajx.Rows; i++ {
			BA.Data[i] = BA.Data[i] + alphajx.Data[i]*A.GetFromMatrix(j, 1)
			BA.Data[i] = BA.Data[i] + betajx.Data[i]*A.GetFromMatrix(j, 2)
		}
	}

	err = true
	return BA, err
}

func BenchmarkInterpHermiteFunc(b *testing.B) {
	//2.4x^3+1.5x^2+0.3x-1.63
	A26 := goNum.NewMatrix(4, 3, []float64{-10.0, -2254.63, 690.3,
		-4.0, -132.43, 103.5,
		4.0, 177.17, 127.5,
		10.0, 2551.37, 750.3})
	// A27 := goNum.NewMatrix(3, 3, []float64{0.0, 0.0, 1.0,
	// 	1.0, 1.0, 1.0,
	// 	2.0, 2.0, 1.0})
	for i := 0; i < b.N; i++ {
		InterpHermiteFunc(A26)
	}
}
