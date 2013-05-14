/*
 * gonum.go, part of gochem.
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation; either version 2.1 of the 
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General 
 * Public License along with this program.  If not, see 
 * <http://www.gnu.org/licenses/>.
 * 
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

//Package chem provides atom and molecule structures, facilities for reading and writing some
//files used in computational chemistry and some functions for geometric manipulations and shape
//indicators.
package chem


import "github.com/skelterjohn/go.matrix"


/*Here I make a -very incomplete- implementation of the gonum api backed by go.matrix, which will enable me to port gochem to gonum. 
 * Since the agreement in the gonum community was NOT to build a temporary implementation, I just build the functions that
 * gochem uses, only the type Float64, and do not export any of the functions.
 * all the names will start with gn (i.e. RandomFunc becomes gnRandomFunc) so its latter easy to use search and replace to set the 
 * correct import path when gonum is implemented (such as gonum.RandomFunc)*/
 
 
 
type CoordMatrix struct {
	*matrix.DenseMatrix
	}
 
 
func NewCoordMatrix(data []float64,rows,cols int) *CoordMatrix{
	return &CoordMatrix{matrix.MakeDenseMatrix(data, rows, cols)}

}
 
//returns and empty, but not nil, coordmatrix
func NewEmptyCoordMatrix() *CoordMatrix{
	var a  *matrix.DenseMatrix
	return &CoordMatrix{a}

} 
 
func Zeros(rows, cols int) *CoordMatrix{
	return &CoordMatrix{matrix.Zeros(rows,cols)}
}

func Eye(ar,ac int) *CoordMatrix{
	A:=&CoordMatrix{matrix.Zeros(ar,ac)}
	for i:=0;i<ar;i++{
		for j:=0;j<ac;j++ {
			if i==j{
				A.Set(i,j,1.0)
			}
		}
		
	}
	return A
}



func (F *CoordMatrix) Row(i int) []float64{
	c,r:=F.Dims()
	if i>=r{
		panic("Requested row out of bounds")
	}
	a:=make([]float64,c,c)
	for j:=0;j<c;j++{
		a[j]=F.At(i,j)
	}
	return a
} 

func (F *CoordMatrix) SomeRows(A *CoordMatrix,clist []int){
	ar,ac:=A.Dims()
	fr,fc:=F.Dims()
	if ac!=fc || fr!=len(clist) || ar<len(clist){
		panic(ErrRowLength)
	}
	r:=make([]float64,len(clist)*ac)
	for key,val:=range(clist){
		for j:=0;j<ac;j++{
			F.Set(key,j,A.Get(val,j))
		}
	} 

//same as before but returns an error instead of panicking
(F *CoordMatrix) SafeSomeRows(A *CoordMatrix,clist []int) error{
	return gnMaybe(F.SomeRows(A,clist))
}

func (F *CoordMatrix) View(A *CoordMatrix,i,j,rows,cols int) {
	*F=CoordMatrix{A.GetMatrix(i,j,rows,cols)}
	}

func (F *CoordMatrix) RowView(i int) *CoordMatrix{
	a:=NewEmptyCoordMatrix()
	a.View(F,i,0,1,3)
	return a
}  



  
  
func (F *CoordMatrix)  Norm(i int)(float64){  
	//temporary hack
	if i != 2 {
		panic("only 2-norm is implemented")
		}
	return F.TwoNorm()
}
 
  

func (F *CoordMatrix)  Clone(A *CoordMatrix){
	ar,ac:=A.Dims()
	/*
	if F==nil{
		F=&CoordMatrix{A.Copy()}
		return
	*/

	fr,fc:=F.Dims()
	if ac != fc || ar != fr {
		panic("receiver must has same shape as to-be-cloned matrix")
	}
	
	for i:=0;i<ar;i++{
		for j:=0;j<ac;j++ {
			F.Set(i,j,A.At(i,j))
			}
		
		}
 
}


func (F *CoordMatrix)  Dims()(int, int){
	return F.Rows(),F.Cols()
}

func (F *CoordMatrix)  At(A, B int)(float64){
	return F.Get(A,B)
}

func (F *CoordMatrix)  MulElem(A, B *CoordMatrix) {
	arows,acols:=A.Dims()
	brows,bcols:=B.Dims()
	if arows != brows || acols != bcols{
		panic("Matrices need to have the same shape for ElemMul")
		}
		for i:=0;i<arows;i++{
			for j:=0;j<acols;j++ {
				F.Set(i,j,A.At(i,j)*B.At(i,j))
				}
			
			}
	}


//Sum returns the sum of all elements in matrix A.
func (F *CoordMatrix)  Sum()  float64 {
	Rows,Cols:=F.Dims()
	var sum float64
	for i := 0; i < Cols; i++ {
		for j := 0; j < Rows; j++ {
			sum += F.Get(j, i)
		}
	}
	return sum
}

//A slightly modified version of John Asmuth's ParalellProduct function. 
func (F *CoordMatrix)  Mul(A, B *CoordMatrix) {
	if A.Cols() != B.Rows() {
		panic("Wrong dimenstions for matrix multiplication")
	}
	Arows,Acols:=A.Dims()
	_,Bcols:=B.Dims()
	
	if F==nil{
	F = Zeros(Arows, Bcols) //I don't know if the final API will allow this.
	}
	
	in := make(chan int)
	quit := make(chan bool)

	dotRowCol := func() {
		for {
			select {
			case i := <-in:
				sums := make([]float64, Bcols)
				for k := 0; k < Acols; k++ {
					for j := 0; j < Bcols; j++ {
						sums[j] += A.At(i, k) * B.At(k, j)
					}
				}
				for j := 0; j < Bcols; j++ {
					F.Set(i, j, sums[j])
				}
			case <-quit:
				return
			}
		}
	}

	threads := 2

	for i := 0; i < threads; i++ {
		go dotRowCol()
	}

	for i := 0; i < Arows; i++ {
		in <- i
	}

	for i := 0; i < threads; i++ {
		quit <- true
	}

	return
}
 
 
func (F *CoordMatrix) Stack(A,B *CoordMatrix)  {
	Arows,Acols:=A.Dims()
	Brows,Bcols:=B.Dims()
	Frows,Fcols:=F.Dims()
	
	if Acols != Bcols || Acols!=Fcols || Arows+Brows!=Frows {
		panic("Mismatched matrices for stacking")
	}

	/*if F==nil{
	F = Zeros(Arows+Brows, Acols) //I don't know if the final API will allow this.
	}*/
	for i:=0;i<Arows+Brows;i++{
		for j:=0;j<Acols;j++{
			if i<Arows{
				F.Set(i,j,A.At(i,j))
			}else{
				F.Set(i,j,B.At(i,j))
			}
		}
	}
	
	return
}
 

//Dot returns the dot product between 2 vectors or matrices
func (F *CoordMatrix) Dot(B *CoordMatrix) float64 {
	var err error
	if F.Cols() != B.Cols() || F.Rows() != B.Rows() {
		panic("Dot: Matrices must have the same dimmension to obtain dot product")
	}
	a,b:=F.Dims()
	A := Zeros(a,b)
	A.MulElem(F,B)
	if err != nil {
		panic(err.Error())
	}
	return F.Sum()
}
 
//Transpose 
func (F *CoordMatrix) T(A *CoordMatrix)  {
	/*if F==nil{
		*F=CoordMatrix{A.Transpose()}
		return
	}*/
	arows,acols:=A.Dims()
	frows,fcols:=F.Dims()
	if arows != fcols || acols != frows{
		panic("Mismatched matrices for Transposing")
	}
	for i:=0;i<arows;i++{
		for j:=0;j<acols;j++ {
			F.Set(i,j,A.At(j,i))
		}
		
	}
} 
 
 
 /**These are from the current proposal for gonum, will be taken out when gonum is implemented**/
 
 // A Panicker is a function that may panic.
type Panicker func()
 
 
 // Maybe will recover a panic with a type matrix.Error from fn, and return this error.
// Any other error is re-panicked.
func gnMaybe(fn Panicker) (err error) {
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			if err, ok = r.(Error); ok {
				return
			}
			panic(r)
		}
	}()
	fn()
	return
}
 
 // Type Error represents matrix package errors. These errors can be recovered by Maybe wrappers.
type gnError string

func (err gnError) gnError() string { return string(err) }
 
 
 const (
	ErrIndexOutOfRange = Error("matrix: index out of range")
	ErrZeroLength      = Error("matrix: zero length in matrix definition")
	ErrRowLength       = Error("matrix: row length mismatch")
	ErrColLength       = Error("matrix: col length mismatch")
	ErrSquare          = Error("matrix: expect square matrix")
	ErrNormOrder       = Error("matrix: invalid norm order for matrix")
	ErrSingular        = Error("matrix: matrix is singular")
	ErrShape           = Error("matrix: dimension mismatch")
	ErrIllegalStride   = Error("matrix: illegal stride")
	ErrPivot           = Error("matrix: malformed pivot list")
)
 
 
 
 
 
 
