// Matrix
/*
------------------------------------------------------
作者   : Black Ghost
日期   : 2018-11-20
版本   : 0.0.0
         0.0.1 2018-12-11 增加切片与矩阵转换
         0.0.2 2018-12-26 增加错误报告
------------------------------------------------------
    矩阵的创建及其操作创建及其简单操作/运算
理论：
    参考 OneThin // http://outofmemery.cn/code-snippet
         /16991/go-language-matrix-operation
    进行了主要运算和结构的补充与修改
------------------------------------------------------
注意事项：
    1. r, c 是从零开始算的
------------------------------------------------------
*/

package goNum

import (
	"fmt"
	"strconv"
)

//数据结构定义----------------------------------------+
//定义Matrix数据类型
type Matrix struct {
	Rows, Columns int       //行数和列数
	Data          []float64 //将矩阵中所有元素作为一维切片
}

//矩阵操作-------------------------------------------+
//通过行列号寻找指定矩阵位置在一维切片中的编号
func findIndex(r, c int, A *Matrix) int {
	//r E [0, n), c E [0, n)
	return r*A.Columns + c
}

//设置指定行列的值
func (A *Matrix) SetMatrix(r, c int, val float64) {
	if (r >= A.Rows) || (c >= A.Columns) {
		panic("Error in goNum.(*Matrix).SetMatrix: Out of range")
	}
	A.Data[findIndex(r, c, A)] = val
}

//获取指定行列的值
func (A *Matrix) GetFromMatrix(r, c int) float64 {
	if (r >= A.Rows) || (c >= A.Columns) {
		panic("Error in goNum.(*Matrix).GetFromMatrix: Out of range")
	}
	return A.Data[findIndex(r, c, A)]
}

//获取指定行的值的切片
func (A *Matrix) RowOfMatrix(i int) []float64 {
	if i >= A.Rows {
		panic("Error in goNum.(*Matrix).RowOfMatrix: Out of range")
	}
	return A.Data[findIndex(i, 0, A):findIndex(i, A.Columns, A)]
}

//获取指定列的值的切片
func (A *Matrix) ColumnOfMatrix(j int) []float64 {
	if j >= A.Columns {
		panic("Error in goNum.(*Matrix).ColumnOfMatrix: Out of range")
	}
	col := make([]float64, A.Rows)
	for i := 0; i < A.Rows; i++ {
		col[i] = A.RowOfMatrix(i)[j]
	}
	return col
}

//矩阵转置
func (A *Matrix) Transpose() Matrix {
	B := ZeroMatrix(A.Columns, A.Rows)
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Columns; j++ {
			B.SetMatrix(j, i, A.GetFromMatrix(i, j))
		}
	}
	return B
}

//格式输出
func (A *Matrix) PrintMatrix() {
	//求出最长字符
	colwidstr := make([]string, A.Columns)
	for i := range colwidstr {
		var maxLen int
		thisColumn := A.ColumnOfMatrix(i)
		for j := range thisColumn {
			thisLen := len(strconv.FormatFloat(thisColumn[j], 'f', -1, 64))
			if thisLen > maxLen {
				maxLen = thisLen
			}
		}
	}
	for i := 0; i < A.Rows; i++ {
		thisRow := A.RowOfMatrix(i)
		fmt.Printf("[")
		for j := range thisRow {
			var format string
			if j == 0 {
				format = "%" + colwidstr[j] + "s"
			} else {
				format = " %" + colwidstr[j] + "s"
			}
			fmt.Printf(format, strconv.FormatFloat(thisRow[j], 'f', -1, 64))
		}
		fmt.Printf("]\n")
	}
}

//矩阵初始化-----------------------------------------+
//r行c列零矩阵
func ZeroMatrix(r, c int) Matrix {
	return Matrix{r, c, make([]float64, r*c)}
}

//n阶单位矩阵
func IdentityE(n int) Matrix {
	A := ZeroMatrix(n, n)
	for i := 0; i < len(A.Data); i += (n + 1) {
		A.Data[i] = 1.0
	}
	return A
}

//以已有数据创建r行c列矩阵
func NewMatrix(r, c int, data []float64) Matrix {
	if len(data) != r*c {
		panic("goNum.Matrix.New: Length of data does not matched r rows and c columns")
	}
	A := ZeroMatrix(r, c)
	A.Data = data
	return A
}

//一维切片转为矩阵(列向量)
func Slices1ToMatrix(s []float64) Matrix {
	A := ZeroMatrix(len(s), 1)
	for i := 0; i < A.Rows; i++ {
		A.Data[i] = s[i]
	}
	return A
}

//二维切片转为矩阵
func Slices2ToMatrix(s [][]float64) Matrix {
	row := len(s)
	col := len(s[0])
	A := ZeroMatrix(row, col)
	for i := 0; i < row; i++ {
		for j := 0; j < col; j++ {
			A.SetMatrix(i, j, s[i][j])
		}
	}
	return A
}

//列向量转为一维切片
func Matrix1ToSlices(A Matrix) []float64 {
	s := make([]float64, A.Rows)
	for i := 0; i < A.Rows; i++ {
		s[i] = A.Data[i]
	}
	return s
}

//二维矩阵转为二维切片
func Matrix2ToSlices(A Matrix) [][]float64 {
	s := make([][]float64, A.Rows)
	for i := 0; i < A.Rows; i++ {
		s[i] = make([]float64, A.Columns)
		for j := 0; j < A.Columns; j++ {
			s[i][j] = A.GetFromMatrix(i, j)
		}
	}
	return s
}

//矩阵运算------------------------------------------+
//矩阵相加
func AddMatrix(A, B Matrix) Matrix {
	if (A.Rows != B.Rows) || (A.Columns != B.Columns) {
		panic("goNum.Matrix.Add: A and B does not matched")
	}
	AaddB := ZeroMatrix(A.Rows, A.Columns)
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Columns; j++ {
			AaddB.SetMatrix(i, j, A.GetFromMatrix(i, j)+B.GetFromMatrix(i, j))
		}
	}
	return AaddB
}

//矩阵相减
func SubMatrix(A, B Matrix) Matrix {
	if (A.Rows != B.Rows) || (A.Columns != B.Columns) {
		panic("goNum.Matrix.Sub: A and B does not matched")
	}
	AsubB := ZeroMatrix(A.Rows, A.Columns)
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < A.Columns; j++ {
			AsubB.SetMatrix(i, j, A.GetFromMatrix(i, j)-B.GetFromMatrix(i, j))
		}
	}
	return AsubB
}

//矩阵数乘
func NumProductMatrix(A Matrix, c float64) Matrix {
	cA := ZeroMatrix(A.Rows, A.Columns)
	for i := 0; i < len(cA.Data); i++ {
		cA.Data[i] = c * A.Data[i]
	}
	return cA
}

//矩阵点乘
func DotPruduct(A, B Matrix) Matrix {
	if A.Columns != B.Rows {
		panic("goNum.Matrix.DotPruduct: A and B does not matched")
	}
	AdotB := ZeroMatrix(A.Rows, B.Columns)
	for i := 0; i < A.Rows; i++ {
		for j := 0; j < B.Columns; j++ {
			for k := 0; k < A.Columns; k++ {
				AdotB.Data[B.Columns*i+j] += A.GetFromMatrix(i, k) * B.GetFromMatrix(k, j)
			}
		}
	}
	return AdotB
}

//向量叉乘，得到垂直于两个向量所在平面的向量
func CrossVector(a, b []float64) []float64 {
	if (len(a) != 3) || (len(b) != 3) {
		panic("goNum.Matrix.CrossVector: vector a or b length is not 3")
	}
	acrossb := make([]float64, 3)
	acrossb[0] = a[1]*b[2] - a[2]*b[1]
	acrossb[1] = a[2]*b[0] - a[0]*b[2]
	acrossb[2] = a[0]*b[1] - a[1]*b[0]
	return acrossb
}
