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


import (
	"github.com/skelterjohn/go.matrix"
	"math"
)


/*Here I make a -very incomplete- implementation of the gonum api backed by go.matrix, which will enable me to port gochem to gonum. 
 * Since the agreement in the gonum community was NOT to build a temporary implementation, I just build the functions that
 * gochem uses, on my own type (which should implement all the relevant gonum interfaces).
 * all the gonum-owned names will start with gn (i.e. RandomFunc becomes gnRandomFunc) so its latter easy to use search and replace to set the 
 * correct import path when gonum is implemented (such as gonum.RandomFunc)*/
 
 
 
type CoordMatrix struct {
	*matrix.DenseMatrix
	}
 
 
func NewCoord(data []float64,rows,cols int) *CoordMatrix{
	return &CoordMatrix{matrix.MakeDenseMatrix(data, rows, cols)}

}
 
//returns and empty, but not nil, coordmatrix
func EmptyCoord() *CoordMatrix{
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

//Methods


func (F *CoordMatrix) Add(A, B *CoordMatrix)  {
	ar,ac := A.Dims()
	br,bc := B.Dims()
	fr,fc :=F.Dims()
	if ac != bc || br != ar || ac != fc || ar != fr {
		panic(gnErrShape)
	}
	for i:=0;i<fr;i++{
		for j:=0;j<fc;j++ {
			F.Set(i,j,A.At(i,j)+B.At(i,j))
		}
	}	
	
}


//AddRow adds the row vector row to each row of the matrix A, putting 
//the result on the receiver. Panics if matrices are mismatched.
func (F *CoordMatrix) AddRow(A, row *CoordMatrix)  {
	ar,ac := A.Dims()
	rr,rc :=row.Dims()
	fr,fc :=F.Dims()
	if ac != rc || rr != 1 || ac != fc || ar != fr {
		panic(gnErrShape)
	}
	j:=EmptyCoord()
	for i := 0; i < ar; i++ {
		j.RowView(A,i)
		j.Add(j,row)
	}
}

//AddRow subtracts the row vector row to each row of the matrix A, putting 
//the result on the receiver. Panics if matrices are mismatched.
func (F *CoordMatrix) SubRow(A, row *CoordMatrix)  {
	row.Scale(-1,row)
	F.AddRow(A, row)
	row.Scale(-1,row)
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
		panic(gnErrShape)
	}
	for key,val:=range(clist){
		for j:=0;j<ac;j++{
			F.Set(key,j,A.At(val,j))
		}
	} 
}	
	
	

//same as before but returns an error instead of panicking.
func (F *CoordMatrix) SomeRowsSafe(A *CoordMatrix,clist []int) (err error) {
	f:=func(){F.SomeRows(A,clist)}
	return gnMaybe(gnPanicker(f))
}


func (F *CoordMatrix) SetRows(A *CoordMatrix,clist []int){
	ar,ac:=A.Dims()
	fr,fc:=F.Dims()
	if ac!=fc || fr<len(clist) || ar<len(clist){
		panic(gnErrShape)
	}
	for key,val:=range(clist){
		for j:=0;j<ac;j++{
			F.Set(val,j,A.Get(key,j))
		}
	} 
} 




func (F *CoordMatrix) View(A *CoordMatrix,i,j,rows,cols int) {
	*F=CoordMatrix{A.GetMatrix(i,j,rows,cols)}
	}

//Puts a view of the matrix on the receiver
func (F *CoordMatrix) RowView(A *CoordMatrix, i int){
	F.View(A,i,0,1,3)
}  

//Unit takes a vector and divides it by its norm
//thus obtaining an unitary vector pointing in the same direction as 
//vector.
func (F *CoordMatrix) Unit(A *CoordMatrix) {
	norm := 1.0/A.Norm(2)
	F.Clone(A)
	F.Scale(norm,F)
}


func  (F *CoordMatrix) scaleAux(factor float64) {
	fr,fc:=F.Dims()
	for i:=0;i<fr;i++{
		for j:=0;j<fc;j++ {
			F.Set(i,j,F.At(i,j)*factor)
		}
		
	}
}

func (F *CoordMatrix) Scale(i float64, A *CoordMatrix) {
	if A == F{ //if A and F points to the same object.
		F.scaleAux(i)
	}else{
		F.Clone(A)
		F.scaleAux(i)
	}
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
	fr,fc:=F.Dims()
	if ac != fc || ar != fr {
		panic(gnErrShape)
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
		panic(gnErrShape)
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
		panic(gnErrShape)
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
		panic(gnErrShape)
	}


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
		panic(gnErrShape)
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
	arows,acols:=A.Dims()
	frows,fcols:=F.Dims()
	if arows != fcols || acols != frows{
		panic("gnErrShape")
	}
	for i:=0;i<arows;i++{
		for j:=0;j<acols;j++ {
			F.Set(i,j,A.At(j,i))
		}
		
	}
} 
 
func (F *CoordMatrix) Pow(A *CoordMatrix, exp float64)  {
	ar,ac:=A.Dims()
	fr,fc:=F.Dims()
	if ar!=fr || ac!=fc{
		panic(gnErrShape)
		}
	for i:=0;i<ar;i++{
		for j:=0;j<ac;j++ {
			F.Set(i,j,math.Pow(A.At(i,j),exp))
		}
		
	}
}
 
 
 /**These are from the current proposal for gonum, by Dan Kortschak. It will be taken out
  * from here when gonum is implemented**/
 
 // A Panicker is a function that may panic.
type gnPanicker func()
 
 
 // Maybe will recover a panic with a type matrix.Error from fn, and return this error.
// Any other error is re-panicked.
func gnMaybe(fn gnPanicker) (err error) {
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			if err, ok = r.(gnError); ok {
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

func (err gnError) Error() string { return string(err) }
 
 
 const (
	gnErrIndexOutOfRange = gnError("matrix: index out of range")
	gnErrZeroLength      = gnError("matrix: zero length in matrix definition")
	gnErrRowLength       = gnError("matrix: row length mismatch")
	gnErrColLength       = gnError("matrix: col length mismatch")
	gnErrSquare          = gnError("matrix: expect square matrix")
	gnErrNormOrder       = gnError("matrix: invalid norm order for matrix")
	gnErrSingular        = gnError("matrix: matrix is singular")
	gnErrShape           = gnError("matrix: dimension mismatch")
	gnErrIllegalStride   = gnError("matrix: illegal stride")
	gnErrPivot           = gnError("matrix: malformed pivot list")
)
 
 
 
 
 
 
