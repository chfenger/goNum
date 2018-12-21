// VectorRotation
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-12-20
版本   : 0.0.0
------------------------------------------------------
    向量在三维空间的旋转
理论：

------------------------------------------------------
输入   :
    u       初始向量，3x1
    angle   旋转角度，3x1，按绕x、y、z顺序，弧度
    旋转顺序 []int
输出   :
    sol     解,向量
    err     解出标志：false-未解出或达到步数上限；
                     true-全部解出
------------------------------------------------------
*/

package goNum

import (
	"math"
)

func VectorRotation(u, angle Matrix, seq []int) (Matrix, bool) {
	/*
		向量在三维空间的旋转
		输入   :
		    u       初始向量，3x1
		    angle   旋转角度，3x1，按绕x、y、z顺序，弧度
		    旋转顺序 []int
		输出   :
		    sol     解,向量
		    err     解出标志：false-未解出或达到步数上限；
		                     true-全部解出
	*/
	//向量大小
	if u.Rows != 3 {
		panic("Error in goNum.VectorRotation: Vector length is not right")
	}
	if angle.Rows != 3 {
		panic("Error in goNum.VectorRotation: angles number is not right")
	}
	if len(seq) != 3 {
		panic("Error in goNum.VectorRotation: seq length is not right")
	}
	//判断角度大小
	for i := 0; i < 3; i++ {
		if angle.Data[i] > math.Pi/2.0 {
			panic("Error in goNum.VectorRotation: angle value is not right")
		}
	}
	sol := ZeroMatrix(3, 1)
	var err bool = false
	//matrix around x
	Rx := NewMatrix(3, 3, []float64{
		1.0, 0.0, 0.0,
		0.0, math.Cos(angle.Data[0]), -1.0 * math.Sin(angle.Data[0]),
		0.0, math.Sin(angle.Data[0]), math.Cos(angle.Data[0])})
	//matrix around y
	Ry := NewMatrix(3, 3, []float64{
		math.Cos(angle.Data[1]), 0.0, math.Sin(angle.Data[1]),
		0.0, 1.0, 0.0,
		-1.0 * math.Sin(angle.Data[1]), 0.0, math.Cos(angle.Data[1])})
	//matrix around z
	Rz := NewMatrix(3, 3, []float64{
		math.Cos(angle.Data[2]), -1.0 * math.Sin(angle.Data[2]), 0.0,
		math.Sin(angle.Data[2]), math.Cos(angle.Data[2]), 0.0,
		0.0, 0.0, 1.0})
	W := IdentityE(3)
	for i := 0; i < 3; i++ {
		switch seq[i] {
		case 1:
			W = DotPruduct(Rx, W)
		case 2:
			W = DotPruduct(Ry, W)
		case 3:
			W = DotPruduct(Rz, W)
		default:
			panic("Error in goNum.VectorRotation: sequence number is not right")
		}
	}
	sol = DotPruduct(W, u)
	err = true
	return sol, err
}
