// ErrorEvaluation
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-23
版本   : 0.0.0
------------------------------------------------------
    误差估计方法
理论：
    0 最大误差：
    E = max(Abs(f(xk)-y(xk)))

    1 平均误差：
         1   N
    E = --- Sum (Abs(f(xk)-y(xk)))
         N  k=1

    2 均方根误差：
              1
    E = Sqrt(--- Sum (f(xk)-y(xk))^2)
              N

    参考：John H. Mathews and Kurtis D. Fink. Numerical
         methods using MATLAB, 4th ed. Pearson
         Education, 2004. pp. 196
------------------------------------------------------
输入   :
    FY      数据对，nx2，f(xk)---y(xk)
输出   :
    sol     误差结果
------------------------------------------------------
*/

package goNum

import (
	"math"
)

// MaxError 最大误差
func MaxError(FY Matrix) float64 {
	//最大误差
	//判断FY的维数
	if FY.Columns < 2 {
		panic("Error in goNum.MaxError: FY is at least 2 columns")
	}

	errs := ZeroMatrix(FY.Rows, 1)
	var maxE float64
	for i := 0; i < FY.Rows; i++ {
		errs.Data[i] = math.Abs(FY.GetFromMatrix(i, 1) - FY.GetFromMatrix(i, 0))
	}
	maxE, _, _ = Max(errs.Data)
	return maxE
}

// MeanError 平均误差
func MeanError(FY Matrix) float64 {
	//平均误差
	//判断FY的维数
	if FY.Columns < 2 {
		panic("Error in goNum.MaxError: FY is at least 2 columns")
	}

	var meanE float64
	for i := 0; i < FY.Rows; i++ {
		meanE += math.Abs(FY.GetFromMatrix(i, 1) - FY.GetFromMatrix(i, 0))
	}
	meanE = meanE / float64(FY.Rows)
	return meanE
}

// RMSError 均方根误差
func RMSError(FY Matrix) float64 {
	//均方根误差
	//判断FY的维数
	if FY.Columns < 2 {
		panic("Error in goNum.MaxError: FY is at least 2 columns")
	}

	var rmsE float64
	for i := 0; i < FY.Rows; i++ {
		temp0 := FY.GetFromMatrix(i, 1) - FY.GetFromMatrix(i, 0)
		rmsE += temp0 * temp0
	}
	rmsE = math.Sqrt(rmsE / float64(FY.Rows))
	return rmsE
}
