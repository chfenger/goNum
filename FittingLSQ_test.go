// FittingLSQ_test
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-23
版本   : 0.0.0
------------------------------------------------------
    线性最小二乘拟合
理论：
    设对N个数据对的线性拟合表示为
    y = Ax + B

       N            N        N
    A*Sum xi^2 + B*Sum xi = Sum xiyi
      i=1          i=1      i=1
       N             N
    A*Sum xi + NB = Sum yi
      i=1           i=1
    解此二元线性方程组即可得A、B

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. ss 5.1
------------------------------------------------------
输入   :
    XY      数据对，nx2，x-y
输出   :
    sol     解，2x1
    err     解出标志：false-未解出或达到边界；
                     true-全部解出
------------------------------------------------------
*/

package goNum_test

import (
	"testing"

	"github.com/chfenger/goNum"
)

// FittingLSQ 线性最小二乘拟合
func FittingLSQ(XY goNum.Matrix) (goNum.Matrix, bool) {
	/*
		线性最小二乘拟合
		输入   :
		    XY      数据对，nx2，x-y
		输出   :
		    sol     解，2x1
		    err     解出标志：false-未解出或达到边界；
		                     true-全部解出
	*/
	//判断XY的维数
	if XY.Columns < 2 {
		panic("Error in goNum.FittingLSQ: At least 2 columns of XY needed")
	}
	sol := goNum.ZeroMatrix(2, 1)
	AS := goNum.ZeroMatrix(2, 2)
	BS := goNum.ZeroMatrix(2, 1)
	var err bool = false
	var sx2, sx, sxy, sy float64
	n := XY.Rows

	//求累加和
	for i := 0; i < n; i++ {
		sx2 += XY.GetFromMatrix(i, 0) * XY.GetFromMatrix(i, 0)
		sx += XY.GetFromMatrix(i, 0)
		sxy += XY.GetFromMatrix(i, 0) * XY.GetFromMatrix(i, 1)
		sy += XY.GetFromMatrix(i, 1)
	}
	AS.SetMatrix(0, 0, sx2)
	AS.SetMatrix(0, 1, sx)
	AS.SetMatrix(1, 0, sx)
	AS.SetMatrix(1, 1, float64(n))
	BS.SetMatrix(0, 0, sxy)
	BS.SetMatrix(1, 0, sy)

	//解二元线性方程组
	soltemp, errtemp := goNum.LEs_ECPE(goNum.Matrix2ToSlices(AS), goNum.Matrix1ToSlices(BS))
	if errtemp != true {
		panic("Error in goNum.FittingLSQ: Solve error")
	}
	sol.SetMatrix(0, 0, soltemp[1])
	sol.SetMatrix(1, 0, soltemp[0])

	err = true
	return sol, err
}

func BenchmarkFittingLSQ(b *testing.B) {
	xy47 := goNum.NewMatrix(8, 2, []float64{
		-1.0, 10.0,
		0.0, 9.0,
		1.0, 7.0,
		2.0, 5.0,
		3.0, 4.0,
		4.0, 3.0,
		5.0, 0.0,
		6.0, -1.0})
	for i := 0; i < b.N; i++ {
		goNum.FittingLSQ(xy47)
	}
}
