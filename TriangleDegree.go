// TriangleDegree
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-13
版本   : 0.0.0
------------------------------------------------------
    以角度为输入的三角函数计算
理论：

------------------------------------------------------
输入   :
    x      角度值
    y      数值
输出   :
    ***    数值或角度值
------------------------------------------------------
*/

package goNum

import (
	"math"
)

//三角函数
//角度的正弦
func Sind(x float64) float64 {
	return math.Sin(x * math.Pi / 180.0)
}

//角度的余弦
func Cosd(x float64) float64 {
	return math.Cos(x * math.Pi / 180.0)
}

//角度的正切
func Tand(x float64) float64 {
	return math.Tan(x * math.Pi / 180.0)
}

//反三角函数
//反正弦的角度
func Asind(y float64) float64 {
	return 180.0 * math.Asin(y) / math.Pi
}

//反余弦的角度
func Cosd(y float64) float64 {
	return 180.0 * math.Acos(y) / math.Pi
}

//反正切的角度
func Tand(y float64) float64 {
	return 180.0 * math.Atan(y) / math.Pi
}
