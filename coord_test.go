/*
 * coord_test.go
 * 
 * Copyright 2013 Raul Mera <rmera@zinc>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */



package chem

import "testing"
import "fmt"

func TestGeo(Te *testing.T) {
	a:=[]float64{1.0,2.0,3,4,5,6,7,8,9}
	A:=NewCoordMatrix(a,3,3)
	ar,ac:=A.Dims()
	T:=Zeros(ar,ac)
	T.T(A)
	B:=Eye(ar,ac)
	//B.Clone(A)
	T.Mul(A,B)
	E:=Zeros(ar,ac)
	E.MulElem(A,B)
	fmt.Println(T,"\n",T,"\n",A,"\n",B,"\n",ar,ac,A.Sum())
	//View:=Zeros(1,1)
	View:=A.AtomCoords(0)
	View.Set(0,0,100)
	fmt.Println("View\n",A,"\n",View)
	
}
