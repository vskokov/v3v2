// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#ifndef __ZHEEVV3_H
#define __ZHEEVV3_H

#include <complex>
#include <blitz/blitz.h>
#include <blitz/array.h>



using namespace std;
typedef complex<double> cd;
typedef blitz::TinyMatrix<cd,3,3>  colorMat;

int zheevv3(colorMat& A, colorMat& Q, blitz::TinyVector<double,3>& w);

#endif
