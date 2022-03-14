/*
	File: MC33.cpp
	Programmed by: David Vega - dvega@uc.edu.ve
	August 2019
	February 2020
	July 2020
	August 2020
	October 2020
	April 2021
	June 2021
	August 2021
	December 2021
	January 2022
	March 2022
*/

#ifdef _MSC_VER
#pragma warning( push )
// disable warning when a double value is assigned to float variable
#pragma warning( disable : 4244 )
#endif

#include <cstring>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "../include/MC33.h"

#if defined (__SSE__) || ((defined (_M_IX86) || defined (_M_X64)) && !defined(_CHPE_ONLY_))
// https://stackoverflow.com/questions/59644197/inverse-square-root-intrinsics
// faster than 1.0f/std::sqrt, but with little accuracy.
#include <immintrin.h>
inline float invSqrt(float f) {
	__m128 temp = _mm_set_ss(f);
	temp = _mm_rsqrt_ss(temp);
	return _mm_cvtss_f32(temp);
}
#else
#include <cmath>
inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}
#endif

using namespace std;

MC33::MC33() : F(0), table
#include "MC33_LookUpTable.h"
{}

MC33::~MC33() {
	clear_temp_isosurface();
}

void MC33::set_default_surface_color(unsigned char *color) {
	DefaultColor = *(reinterpret_cast<int*>(color));
}

/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function returns a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test should not be applied. The
return value of this function is the the sum of all six results.*/
int MC33::face_tests(int *face, int i) const {
	if (i&0x80) { //vertex 0
		face[0] = ((i&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((i&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((i&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//0xA0 = 10100000, vertices 0 and 2
	} else {
		face[0] = ((i&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((i&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((i&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if (i&0x02) { //vertex 6
		face[1] = ((i&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((i&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((i&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//0x0A = 00001010, vertices 4 and 6
	} else {
		face[1] = ((i&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((i&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((i&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}

/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int MC33::face_test1(int face) const {
	switch (face) {
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	default:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
}

// Silence dereferencing type-punned pointer warning in GCC
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif
// an ugly signbit:
#if MC33_double_precision
inline unsigned int signbf(double x) {
	return ((*(reinterpret_cast<unsigned long long int*>(&x))>>32)&0x80000000);
}
#else
inline unsigned int signbf(float x) {
	return (*(reinterpret_cast<unsigned int*>(&x))&0x80000000);
}
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different from 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flag13 must be 1, and the function returns 2 if one
of the vertices 0, 1, 2 or 3 is joined to the center point of the cube (case
13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are not joined (case 13.5.1)*/
int MC33::interior_test(int i, int flag13) const {
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	MC33_real At = v[4] - v[0], Bt = v[5] - v[1],
				Ct = v[6] - v[2], Dt = v[7] - v[3];
	MC33_real t = At*Ct - Bt*Dt;//the "a" value.
	if (signbf(t)) {
		if (i&0x01) return 0;
	} else {
		if (!(i&0x01) || t == 0) return 0;
	}
	t = 0.5f*(v[3]*Bt + v[1]*Dt - v[2]*At - v[0]*Ct)/t;//t = -b/2a

	if (t > 0 && t < 1) {
		At = v[0] + At*t;
		Bt = v[1] + Bt*t;
		Ct = v[2] + Ct*t;
		Dt = v[3] + Dt*t;
		Ct *= At;
		Dt *= Bt;
		if (i&0x01) {
			if (Ct < Dt && signbf(Dt) == 0)
				return (signbf(Bt) == signbf(v[i])) + flag13;
		} else {
			if (Ct > Dt && signbf(Ct) == 0)
				return (signbf(At) == signbf(v[i])) + flag13;
		}
	}
	return 0;
}

/******************************************************************
Auxiliary function that calculates the normal if a vertex
of the cube lies on the isosurface.
*/
unsigned int MC33::surfint(unsigned int x, unsigned int y, unsigned int z, MC33_real *r) {
	r[0] = x; r[1] = y; r[2] = z;
	if (x == 0)
		r[3] = F[z][y][0] - F[z][y][di];
	else if (x == nx) {
		x *= di;
		r[3] = F[z][y][x - di] - F[z][y][x];
	} else {
		x *= di;
		r[3] = 0.5f*(F[z][y][x - di] - F[z][y][x + di]);
	}
	if (y == 0)
		r[4] = F[z][0][x] - F[z][1][x];
	else if (y == ny)
		r[4] = F[z][y - 1][x] - F[z][y][x];
	else
		r[4] = 0.5f*(F[z][y - 1][x] - F[z][y + 1][x]);
	if (z == 0)
		r[5] = F[0][y][x] - F[1][y][x];
	else if (z == nz)
		r[5] = F[z - 1][y][x] - F[z][y][x];
	else
		r[5] = 0.5f*(F[z - 1][y][x] - F[z + 1][y][x]);
	return store_point(r);
}

/******************************************************************
This function find the MC33 case (using the index i, and the face and interior
tests). The correct triangle pattern is obtained from the arrays contained in
the MC33_LookUpTable.h file. The necessary vertices (intersection points) are
also calculated here (using trilinear interpolation).
       _____2_____
     /|          /|
   11 |<-3     10 |
   /____6_____ /  1     z
  |   |       |   |     |
  |   |_____0_|___|     |____y
  7  /        5  /     /
  | 8         | 9     x
  |/____4_____|/

The temporary matrices: Dx, Dy, Ux, Uy, and Lz are filled
and used here.*/
#define FF 0xFFFFFFFF
void MC33::find_case(unsigned int x, unsigned int y, unsigned int z, unsigned int i) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	unsigned int ti[3];//for vertex indices of a triangle
	union { // memory saving
		int f[6];//for the face tests
		MC33_real r[6];//for intercept and normal coordinates
	};
	const unsigned short int *pcase = table;
	unsigned int k, c, m, n;
	MC33_real t;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
		n = !m;
	} else {
		c = pcase[i];
		n = (c&0x800) == 0;
		m = !n;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&face_test1(k>>2)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (interior_test(k,0)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&face_test1(k%6))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (interior_test(k/6,0)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (f[4] + f[5] < 0? (f[0] + f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (f[4] + f[5] < 0? 983: (f[0] + f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (interior_test(k>>1,0)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2:
				if (k&2? interior_test(0,0): interior_test(0,0)||interior_test(k? 1: 3,0))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k&2? interior_test(1,0): interior_test(2,0)||interior_test(k? 3: 1,0))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2: //12.1
				pcase += (interior_test((0xDA010C>>(k<<1))&3,0)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (interior_test((0xA7B7E5>>(k<<1))&3,0)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(face_tests(f, 165))) {
			case 0:
				k = ((f[1] < 0)<<1)|(f[5] < 0);
				if (f[0]*f[1] == f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = interior_test(k, 1); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((f[0] < 0? f[2] > 0: 12 + (f[2] < 0)) + (f[1] < 0? f[3] < 0: 6 + (f[3] > 0)));
				if (f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*f[0] + 4*f[1] + 3*f[2] + 2*f[3] + f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) { // the vertices r[3] and normals (r + 3)[3] are calculated here
				case 0:
					if (z || x)
						ti[--k] = p[0] = Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = Dx[y][0];
							else if (y? signbf(S->iso - F[0][y - 1][0]): 0)
								p[0] = Dy[y - 1][0];
							else
								p[0] = surfint(0,y,0,r);
						} else if (v[1] == 0) {
							if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = (p[1] != FF? p[1]: surfint(0,y + 1,0,r));
						} else {
							t = v[0]/(v[0] - v[1]);
							r[0] = r[2] = 0;
							r[1] = y + t;
							r[3] = (v[4] - v[0])*(1 - t) + (v[5] - v[1])*t;
							r[4] = v[1] - v[0];
							r[5] = (v[3] - v[0])*(1 - t) + (v[2] - v[1])*t;
							p[0] = store_point(r);
						}
						Dy[y][0] = ti[--k] = p[0];
					}
					break;
				case 1:
					if (x)
						ti[--k] = p[1] = Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = Dy[y][0];
							//else if (z && signbf(v[5]))
							//	p[1] = Dx[y + 1][0];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][0]): 0)
								p[1] = Dy[y + 1][0];
							else if (z? signbf(S->iso - F[z - 1][y + 1][0]): 0)
								p[1] = Lz[y + 1][0]; // value of previous slice
							else
								p[1] = surfint(0,y + 1,z,r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = (p[2] != FF? p[2]: surfint(0,y + 1,z + 1,r));
						} else {
							t = v[1]/(v[1] - v[2]);
							r[0] = 0; r[1] = y + 1;
							r[2] = z + t;
							r[3] = (v[5] - v[1])*(1 - t) + (v[6] - v[2])*t;
							r[4] = (y + 1 < ny? 0.5f*((F[z][y][0] - F[z][y + 2][0])*(1 - t)
										+ (F[z + 1][y][0] - F[z + 1][y + 2][0])*t):
										(v[1] - v[0])*(1 - t) + (v[2] - v[3])*t);
							r[5] = v[2] - v[1];
							p[1] = store_point(r);
						}
						Lz[y + 1][0] = ti[--k] = p[1];
					}
					break;
				case 2:
					if (x)
						ti[--k] = p[2] = Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = Ux[y][0];
							else if (y? signbf(S->iso - F[z + 1][y - 1][0]): 0)
								p[2] = Uy[y - 1][0];
							else
								p[2] = surfint(0,y,z + 1,r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = (p[1] != FF? p[1]: surfint(0,y + 1,z + 1,r));
						} else {
							t = v[3]/(v[3] - v[2]);
							r[0] = 0; r[2] = z + 1;
							r[1] = y + t;
							r[3] = (v[7] - v[3])*(1 - t) + (v[6] - v[2])*t;
							r[4] = v[2] - v[3];
							r[5] = (z + 1 < nz? 0.5f*((F[z][y][0] - F[z + 2][y][0])*(1 - t)
										+ (F[z][y + 1][0] - F[z + 2][y + 1][0])*t):
										(v[3] - v[0])*(1 - t) + (v[2] - v[1])*t);
							p[2] = store_point(r);
						}
						Uy[y][0] = ti[--k] = p[2];
					}
					break;
				case 3:
					if (y || x)
						ti[--k] = p[3] = Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = Dx[0][0];
							else if (z? signbf(S->iso - F[z - 1][0][0]): 0)
								p[3] = Lz[0][0]; // value of previous slice
							else
								p[3] = surfint(0,0,z,r);
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else
								p[3] = (p[11] != FF? p[11]: surfint(0,0,z + 1,r));
						} else {
							t = v[0]/(v[0] - v[3]);
							r[0] = r[1] = 0;
							r[2] = z + t;
							r[3] = (v[4] - v[0])*(1 - t) + (v[7] - v[3])*t;
							r[4] = (v[1] - v[0])*(1 - t) + (v[2] - v[3])*t;
							r[5] = v[3] - v[0];
							p[3] = store_point(r);
						}
						Lz[0][0] = ti[--k] = p[3];
					}
					break;
				case 4:
					if (z)
						ti[--k] = p[4] = Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[7]))
								p[4] = Lz[y][x + 1];
							else if (y && signbf(v[0]))
								p[4] = Dx[y][x];
							else if (y? signbf(S->iso - F[0][y - 1][di*(x + 1)]): 0)
								p[4] = Dy[y - 1][x + 1];
							else if (y && x + 1 < nx? signbf(S->iso - F[0][y][di*(x + 2)]): 0)
								p[4] = Dx[y][x + 1];
							else
								p[4] = surfint(x + 1,y,0,r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = (p[9] != FF? p[9]: surfint(x + 1,y + 1,0,r));
						} else {
							t = v[4]/(v[4] - v[5]);
							r[0] = x + 1; r[2] = 0;
							r[1] = y + t;
							r[3] = (x + 1 < nx? 0.5f*((F[0][y][di*x] - F[0][y][di*(x + 2)])*(1 - t)
										+ (F[0][y + 1][di*x] - F[0][y + 1][di*(x + 2)])*t):
										(v[4] - v[0])*(1 - t) + (v[5] - v[1])*t);
							r[4] = v[5] - v[4];
							r[5] = (v[7] - v[4])*(1 - t) + (v[6] - v[5])*t;
							p[4] = store_point(r);
						}
						Dy[y][x + 1] = ti[--k] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = Dx[y + 1][x];
							else if (x + 1 < nx? signbf(S->iso - F[z][y + 1][di*(x + 2)]): 0)
								p[5] = Dx[y + 1][x + 1];
							else if (y + 1 < ny? signbf(S->iso - F[z][y + 2][di*(x + 1)]): 0)
								p[5] = Dy[y + 1][x + 1];
							else if (signbf(S->iso - F[z - 1][y + 1][di*(x + 1)])) {
								p[5] = Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = surfint(x + 1,y + 1,z,r);
						} else
							p[5] = surfint(x + 1,y + 1,0,r);
					} else if (v[6] == 0) {
						p[5] = surfint(x + 1,y + 1,z + 1,r);
					} else {
						t = v[5]/(v[5] - v[6]);
						r[0] = x + 1; r[1] = y + 1;
						r[2] = z + t;
						c = di*x;
						r[3] = (x + 1 < nx? 0.5f*((F[z][y + 1][c] - F[z][y + 1][c + 2*di])*(1 - t)
									+ (F[z + 1][y + 1][c] - F[z + 1][y + 1][c + 2*di])*t):
									(v[5] - v[1])*(1 - t) + (v[6] - v[2])*t);
						r[4] = (y + 1 < ny? 0.5f*((F[z][y][c + di] - F[z][y + 2][c + di])*(1 - t)
									+ (F[z + 1][y][c + di] - F[z + 1][y + 2][c + di])*t):
									(v[5] - v[4])*(1 - t) + (v[6] - v[7])*t);
						r[5] = v[6] - v[5];
						p[5] = store_point(r);
					}
					ti[--k] = Lz[y + 1][x + 1] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = Lz[y][x + 1];
							else if (signbf(S->iso - F[z + 1][y - 1][di*(x + 1)]))
								p[6] = Uy[y - 1][x + 1];
							else if (x + 1 < nx? signbf(S->iso - F[z + 1][y][di*(x + 2)]): 0)
								p[6] = Ux[y][x + 1];
							else
								p[6] = surfint(x + 1,y,z + 1,r);
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = surfint(x + 1,0,z + 1,r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[6] = p[5];
						else
							p[6] = (p[10] != FF? p[10]: surfint(x + 1,y + 1,z + 1,r));
					} else {
						t = v[7]/(v[7] - v[6]);
						r[0] = x + 1;
						r[1] = y + t; r[2] = z + 1;
						c = di*x;
						r[3] = (x + 1 < nx? 0.5f*((F[z + 1][y][c] - F[z + 1][y][c + 2*di])*(1 - t)
									+ (F[z + 1][y + 1][c] - F[z + 1][y + 1][c + 2*di])*t):
									(v[7] - v[3])*(1 - t) + (v[6] - v[2])*t);
						r[4] = v[6] - v[7];
						r[5] = (z + 1 < nz? 0.5f*((F[z][y][c + di] - F[z + 2][y][c + di])*(1 - t)
										+ (F[z][y + 1][c + di] - F[z + 2][y + 1][c + di])*t):
									(v[7] - v[4])*(1 - t) + (v[6] - v[5])*t);
						p[6] = store_point(r);
					}
					ti[--k] = Uy[y][x + 1] = p[6];
					break;
				case 7:
					if (y)
						ti[--k] = p[7] = Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = Dx[0][x];
							//else if (z && signbf(v[5]))
							//	p[7] = Dy[0][x + 1];
							else if (z && x + 1 < nx? signbf(S->iso - F[z][0][di*(x + 2)]): 0)
								p[7] = Dx[0][x + 1];
							else if (z? signbf(S->iso - F[z - 1][0][di*(x + 1)]): 0)
								p[7] = Lz[0][x + 1]; // value of previous slice
							else
								p[7] = surfint(x + 1,0,z,r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else
								p[7] = (p[11] != FF? p[11]: surfint(x + 1,0,z + 1,r));
						} else {
							t = v[4]/(v[4] - v[7]);
							r[0] = x + 1; r[1] = 0;
							r[2] = z + t;
							r[3] = (x + 1 < nx? 0.5f*((F[z][0][di*x] - F[z][0][di*(x + 2)])*(1 - t)
										+ (F[z + 1][0][di*x] - F[z + 1][0][di*(x + 2)])*t):
										(v[4] - v[0])*(1 - t) + (v[7] - v[3])*t);
							r[4] = (v[5] - v[4])*(1 - t) + (v[6] - v[7])*t;
							r[5] = v[7] - v[4];
							p[7] = store_point(r);
						}
						Lz[0][x + 1] = ti[--k] = p[7];
					}
					break;
				case 8:
					if (z || y)
						ti[--k] = p[8] = Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = Dy[0][x];
							else if (x? signbf(S->iso - F[0][0][di*(x - 1)]): 0)
								p[8] = Dx[0][x - 1];
							else
								p[8] = surfint(x,0,0,r);
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else
								p[8] = (p[7] != FF? p[7]: surfint(x + 1,0,0,r));
						} else {
							t = v[0]/(v[0] - v[4]);
							r[1] = r[2] = 0;
							r[0] = x + t;
							r[3] = v[4] - v[0];
							r[4] = (v[1] - v[0])*(1 - t) + (v[5] - v[4])*t;
							r[5] = (v[3] - v[0])*(1 - t) + (v[7] - v[4])*t;
							p[8] = store_point(r);
						}
						Dx[0][x] = ti[--k] = p[8];
					}
					break;
				case 9:
					if (z)
						ti[--k] = p[9] = Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = Lz[y + 1][x];
							else if (x? signbf(S->iso - F[0][y + 1][di*(x - 1)]): 0)
								p[9] = Dx[y + 1][x - 1];
							else
								p[9] = surfint(x,y + 1,0,r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else
								p[9] = (p[4] != FF? p[4]: surfint(x + 1,y + 1,0,r));
						} else {
							t = v[1]/(v[1] - v[5]);
							r[1] = y + 1; r[2] = 0;
							r[0] = x + t;
							r[3] = v[5] - v[1];
							r[4] = (y + 1 < ny? 0.5f*((F[0][y][di*x] - F[0][y + 2][di*x])*(1 - t)
										+ (F[0][y][di*(x + 1)] - F[0][y + 2][di*(x + 1)])*t):
										(v[1] - v[0])*(1 - t) + (v[5] - v[4])*t);
							r[5] = (v[2] - v[1])*(1 - t) + (v[6] - v[5])*t;
							p[9] = store_point(r);
						}
						Dx[y + 1][x] = ti[--k] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = Uy[y][x];
							else if (signbf(S->iso - F[z + 1][y + 1][di*(x - 1)]))
								p[10] = Ux[y + 1][x - 1];
							else
								p[10] = surfint(x,y + 1,z + 1,r);
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = surfint(0,y + 1,z + 1,r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[10] = p[5];
						else
							p[10] = (p[6] != FF? p[6]: surfint(x + 1,y + 1,z + 1,r));
					} else {
						t = v[2]/(v[2] - v[6]);
						r[0] = x + t;
						r[1] = y + 1; r[2] = z + 1;
						r[3] = v[6] - v[2];
						c = di*x;
						r[4] = (y + 1 < ny? 0.5f*((F[z + 1][y][c] - F[z + 1][y + 2][c])*(1 - t)
									+ (F[z + 1][y][c + di] - F[z + 1][y + 2][c + di])*t):
									(v[2] - v[3])*(1 - t) + (v[6] - v[7])*t);
						r[5] = (z + 1 < nz? 0.5f*((F[z][y + 1][c] - F[z + 2][y + 1][c])*(1 - t)
									+ (F[z][y + 1][c + di] - F[z + 2][y + 1][c + di])*t):
									(v[2] - v[1])*(1 - t) + (v[6] - v[5])*t);
						p[10] = store_point(r);
					}
					ti[--k] = Ux[y + 1][x] = p[10];
					break;
				case 11:
					if (y)
						ti[--k] = p[11] = Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = Uy[0][x];
							else if (x? signbf(S->iso - F[z + 1][0][di*(x - 1)]): 0)
								p[11] = Ux[0][x - 1];
							else
								p[11] = surfint(x,0,z + 1,r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else
								p[11] = (p[7] != FF? p[7]: surfint(x + 1,0,z + 1,r));
						} else {
							t = v[3]/(v[3] - v[7]);
							r[1] = 0; r[2] = z + 1;
							r[0] = x + t;
							r[3] = v[7] - v[3];
							r[4] = (v[2] - v[3])*(1 - t) + (v[6] - v[7])*t;
							r[5] = (z + 1 < nz? 0.5f*((F[z][0][di*x] - F[z + 2][0][di*x])*(1 - t)
										+ (F[z][0][di*(x + 1)] - F[z + 2][0][di*(x + 1)])*t):
										(v[3] - v[0])*(1 - t) + (v[7] - v[4])*t);
							p[11] = store_point(r);
						}
						Ux[0][x] = ti[--k] = p[11];
					}
				break;
				default:
					r[0] = x + 0.5f; r[1] = y + 0.5f; r[2] = z + 0.5f;
					r[3] = v[4] + v[5] + v[6] + v[7] - v[0] - v[1] - v[2] - v[3];
					r[4] = v[1] + v[2] + v[5] + v[6] - v[0] - v[3] - v[4] - v[7];
					r[5] = v[2] + v[3] + v[6] + v[7] - v[0] - v[1] - v[4] - v[5];
					ti[--k] = p[12] = store_point(r);
				}
			} else
				ti[--k] = p[c];//now ti contains the vertex indices of the triangle
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2]) { //to avoid zero area triangles
			if (!(S->nT&0x0FFF)) {
				try {
					S->T.resize(S->nT + 0x1000);
				}
				catch (...) {
					memoryfault = 1;
					if (S->nT)
						S->nT = 1;
					return;
				}
			}
			unsigned int *vp = S->T[S->nT++].v;
#ifndef MC_Normal_neg
			*vp = ti[n]; *(++vp) = ti[m]; *(++vp) = ti[2];
#else
			*vp = ti[m]; *(++vp) = ti[n]; *(++vp) = ti[2];
#endif
		}
	}
}

void MC33::case_count(unsigned int x, unsigned int y, unsigned int z, unsigned int i) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	union {
		int f[6];
		unsigned int ti[3];
	};
	const unsigned short int *pcase = table;
	unsigned int c, k, m = i>>7;
	c = pcase[m? i^0xFF: i];
	m = (m^((c&0xFFF)>>11));
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&face_test1(k>>2)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (interior_test(k,0)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&face_test1(k%6))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (interior_test(k/6,0)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (f[4] + f[5] < 0? (f[0] + f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (f[4] + f[5] < 0? 983: (f[0] + f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (interior_test(k>>1,0)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2:
				if (k&2? interior_test(0,0): interior_test(0,0)||interior_test(k? 1: 3,0))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k&2? interior_test(1,0): interior_test(2,0)||interior_test(k? 3: 1,0))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (face_tests(f,(m? i: i^0xFF))) {
			case -2: //12.1
				pcase += (interior_test((0xDA010C>>(k<<1))&3,0)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (interior_test((0xA7B7E5>>(k<<1))&3,0)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(face_tests(f,(m? 90: 165)))) {
			case 0:
				k = ((f[1] < 0)<<1)|(f[5] < 0);
				if (f[0]*f[1] == f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = interior_test(k, 1); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((f[0] < 0? f[2] > 0: 12 + (f[2] < 0)) + (f[1] < 0? f[3] < 0: 6 + (f[3] > 0)));
				if (f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*f[0] + 4*f[1] + 3*f[2] + 2*f[3] + f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) {
				case 0:
					if (z || x)
						p[0] = Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = Dx[y][0];
							else if (y? signbf(S->iso - F[0][y - 1][0]): 0)
								p[0] = Dy[y - 1][0];
							else
								p[0] = S->nV++;;
						} else if (v[1] == 0) {
							if (p[1] != FF)
								p[0] = p[1];
							else if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = S->nV++;;
						} else
							p[0] = S->nV++;;
						Dy[y][0] = p[0];
					}
					break;
				case 1:
					if (x)
						p[1] = Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = Dy[y][0];
							else if (z && signbf(v[5]))
								p[1] = Dx[y + 1][0];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][0]): 0)
								p[1] = Dy[y + 1][0];
							else if (z? signbf(S->iso - F[z - 1][y + 1][0]): 0)
								p[1] = Lz[y + 1][0];
							else
								p[1] = S->nV++;
						} else if (v[2] == 0) {
							if (p[2] != FF)
								p[1] = p[2];
							else if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = S->nV++;
						} else
							p[1] = S->nV++;
						Lz[y + 1][0] = p[1];
					}
					break;
				case 2:
					if (x)
						p[2] = Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = Ux[y][0];
							else if (y? signbf(S->iso - F[z + 1][y - 1][0]): 0)
								p[2] = Uy[y - 1][0];
							else
								p[2] = S->nV++;
						} else if (v[2] == 0) {
							if (p[1] != FF)
								p[2] = p[1];
							else if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = S->nV++;
						} else
							p[2] = S->nV++;
						Uy[y][0] = p[2];
					}
					break;
				case 3:
					if (y || x)
						p[3] = Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = Dx[0][0];
							else if (z? signbf(S->iso - F[z - 1][0][0]): 0)
								p[3] = Lz[0][0];
							else
								p[3] = S->nV++;
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else if (p[11] != FF)
								p[3] = p[11];
							else
								p[3] = S->nV++;
						} else
							p[3] = S->nV++;
						Lz[0][0] = p[3];
					}
					break;
				case 4:
					if (z)
						p[4] = Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = Lz[y][x + 1];
							else if (y? signbf(S->iso - F[0][y - 1][di*(x + 1)]): 0)
								p[4] = Dy[y - 1][x + 1];
							else if (y && x + 1 < nx? signbf(S->iso - F[0][y][di*(x + 2)]): 0)
								p[4] = Dx[y][x + 1];
							else
								p[4] = S->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else if (p[9] != FF)
								p[4] = p[9];
							else
								p[4] = S->nV++;
						} else
							p[4] = S->nV++;
						Dy[y][x + 1] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (signbf(v[4]))
							p[5] = p[4] = (z? Dy[y][x + 1]: S->nV++);
						else if (signbf(v[1]))
							p[5] = p[9] = (z? Dx[y + 1][x]: S->nV++);
						else {
							if (z && x + 1 < nx? signbf(S->iso - F[z][y + 1][di*(x + 2)]): 0)
								p[5] = Dx[y + 1][x + 1];
							else if (z && y + 1 < ny? signbf(S->iso - F[z][y + 2][di*(x + 1)]): 0)
								p[5] = Dy[y + 1][x + 1];
							else if (z? signbf(S->iso - F[z - 1][y + 1][di*(x + 1)]): 0)
								p[5] = Lz[y + 1][x + 1]; // value of previous slice
							else
								p[5] = S->nV++;
						}
					} else
						p[5] = S->nV++;
					Lz[y + 1][x + 1] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = Lz[y][x + 1];
							else if (signbf(S->iso - F[z + 1][y - 1][di*(x + 1)]))
								p[6] = Uy[y - 1][x + 1];
							else if (x + 1 < nx? signbf(S->iso - F[z + 1][y][di*(x + 2)]): 0)
								p[6] = Ux[y][x + 1];
							else
								p[6] = S->nV++;
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = S->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[6] = (p[10] == FF? S->nV++: p[10]);
						else
							p[6] = p[5];
					} else
						p[6] = S->nV++;
					Uy[y][x + 1] = p[6];
					break;
				case 7:
					if (y)
						p[7] = Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = Dy[0][x + 1];
							else if (z && x + 1 < nx? signbf(S->iso - F[z][0][di*(x + 2)]): 0)
								p[7] = Dx[0][x + 1];
							else if (z? signbf(S->iso - F[z - 1][0][di*(x + 1)]): 0)
								p[7] = Lz[0][x + 1];
							else
								p[7] = S->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else if (p[11] != FF)
								p[7] = p[11];
							else
								p[7] = S->nV++;
						} else
							p[7] = S->nV++;
						Lz[0][x + 1] = p[7];
					}
					break;
				case 8:
					if (z || y)
						p[8] = Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = Dy[0][x];
							else if (x? signbf(S->iso - F[0][0][di*(x - 1)]): 0)
								p[8] = Dx[0][x - 1];
							else
								p[8] = S->nV++;
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else if (p[7] != FF)
								p[8] = p[7];
							else
								p[8] = S->nV++;
						} else
							p[8] = S->nV++;
						Dx[0][x] = p[8];
					}
					break;
				case 9:
					if (z)
						p[9] = Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = Lz[y + 1][x];
							else if (x? signbf(S->iso - F[0][y + 1][di*(x - 1)]): 0)
								p[9] = Dx[y + 1][x - 1];
							else
								p[9] = S->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else if (p[4] != FF)
								p[9] = p[4];
							else
								p[9] = S->nV++;
						} else
							p[9] = S->nV++;
						Dx[y + 1][x] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = Uy[y][x];
							else if (signbf(S->iso - F[z + 1][y + 1][di*(x - 1)]))
								p[10] = Ux[y + 1][x - 1];
							else
								p[10] = S->nV++;
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = S->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[10] = (p[6] == FF? S->nV++: p[6]);
						else
							p[10] = p[5];
					} else
						p[10] = S->nV++;
					Ux[y + 1][x] = p[10];
					break;
				case 11:
					if (y)
						p[11] = Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = Uy[0][x];
							else if (x? signbf(S->iso - F[z + 1][0][di*(x - 1)]): 0)
								p[11] = Ux[0][x - 1];
							else
								p[11] = S->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else if (p[7] != FF)
								p[11] = p[7];
							else
								p[11] = S->nV++;
						} else
							p[11] = S->nV++;
						Ux[0][x] = p[11];
					}
				break;
				default:
					p[12] = S->nV++;
				}
			}
			ti[--k] = p[c];
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2])
			++S->nT;
	}
}
#undef FF

void MC33::free_temp_D_U() {
	delete[] Dx; delete[] Ux; delete[] Dy; delete[] Uy;
	delete[] Lz;
}

void MC33::clear_temp_isosurface() {
	if (F) {
		for (unsigned int y = 0; y != ny; ++y) {
			delete[] Dx[y]; delete[] Ux[y]; delete[] Dy[y]; delete[] Uy[y];
			delete[] Lz[y];
		}
		delete[] Dx[ny]; delete[] Ux[ny];
		delete[] Lz[ny];
		free_temp_D_U();
		F = 0;
	}
}

#define code_in_define_part01\
			unsigned int nv = S->nV++;\
			if (!(nv&0x0FFF)) {\
				try {\
					S->V.resize(nv + 0x1000);\
					S->N.resize(nv + 0x1000);\
				}\
				catch (...) {\
					memoryfault = 1;\
					S->nV = (nv? 1: 0);\
					return 0u;\
				}\
			}

#ifndef MC_Normal_neg
#define code_in_define_part02\
			MC33_real t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);\
			float *q = S->N[nv].v;\
			*q = t * *r; *(++q) = t * *(++r); *(++q) = t * *(++r);\
			return nv;
#else //MC_Normal_neg reverses the direction of the normal
#define code_in_define_part02\
			MC33_real t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);\
			float *q = S->N[nv].v;\
			*q = t * *r; *(++q) = t * *(++r); *(++q) = t * *(++r);\
			return nv;
#endif

int MC33::set_grid3d(grid3d &G) {
	clear_temp_isosurface();
	nx = G.N[0];
	ny = G.N[1];
	nz = G.N[2];
#ifndef GRD_orthogonal
	if (G.nonortho) {
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			multAbf(_A,r,r,0);
			for (int i = 0; i != 3; ++i)
				S->V[nv].v[i] = *(r++) + MC_O[i];
			multAbf(A_,r,r,1);
			code_in_define_part02
		};
		for (int j = 0; j != 3; ++j)
			for (int i = 0; i != 3; ++i) {
				_A[j][i] = G._A[j][i]*G.d[i]; // true transformation matrices
				A_[j][i] = G.A_[j][i]/G.d[j];
			}
	} else
#endif
	if (G.d[0] != G.d[1] || G.d[1] != G.d[2]) {
		ca = G.d[2]/G.d[0];
		cb = G.d[2]/G.d[1];
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; ++i)
				S->V[nv].v[i] = *(r++)*MC_D[i] + MC_O[i];
			r[0] *= ca; // normal[0]
			r[1] *= cb; // normal[1]
			code_in_define_part02
		};
	} else if (G.d[0] == 1 && G.r0[0] == 0 && G.r0[1] == 0 && G.r0[2] == 0)
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; ++i)
				S->V[nv].v[i] = *(r++);
			code_in_define_part02
		};
	else
		store_point = [this] (MC33_real* r) {
			code_in_define_part01
			for (int i = 0; i != 3; ++i)
				S->V[nv].v[i] = *(r++)*MC_D[i] + MC_O[i];
			code_in_define_part02
		};

	for (int j = 0; j != 3; ++j) {
		MC_O[j] = G.r0[j];
		MC_D[j] = G.d[j];
	}
	F = const_cast<const GRD_data_type***>(G.F);
	if (!F)
		return 0;
	Lz = new (nothrow) unsigned int*[ny + 1]; //edges 1, 3, 5 (only write) and 7
	Dy = new (nothrow) unsigned int*[ny]; //edges 0 and 4
	Uy = new (nothrow) unsigned int*[ny]; //edges 2 and 6 (only write)
	Dx = new (nothrow) unsigned int*[ny + 1]; //edges 8 and 9
	Ux = new (nothrow) unsigned int*[ny + 1]; //edges 10 (only write) and 11
	if (!Ux) {
		free_temp_D_U();
		return -1;
	}
	for (unsigned int j = 0; j != ny; ++j) {
		Dx[j] = new (nothrow) unsigned int[nx];
		Ux[j] = new (nothrow) unsigned int[nx];
		Lz[j] = new (nothrow) unsigned int[nx + 1];
		Dy[j] = new (nothrow) unsigned int[nx + 1];
		Uy[j] = new (nothrow) unsigned int[nx + 1];
	}
	if (Uy[ny - 1]) {
		Dx[ny] = new (nothrow) unsigned int[nx];
		Ux[ny] = new (nothrow) unsigned int[nx];
		Lz[ny] = new (nothrow) unsigned int[nx + 1];
		di = G.x_data > 1? G.x_data: 1;
		if (Lz[ny])
			return 0;
	} else
		Dx[ny] = Ux[ny] = Lz[ny] = 0;
	clear_temp_isosurface();
	return -1;
}
#undef code_in_define_part01
#undef code_in_define_part02

int MC33::set_grid3d(grid3d *G) {
	if (!G)
		return -1;
	return set_grid3d(*G);
}

int MC33::calculate_isosurface(surface &Sf, MC33_real iso) {
	const GRD_data_type ***FG = F;
	if (!FG)//The set_grid3d function was not executed
		return -2;
	if (Sf.nV)
		Sf.clear();
	S = &Sf;
	Sf.iso = iso;
	memoryfault = 0;
	unsigned int d = di, Nx = nx;
	MC33_real Vt[12];
	MC33_real *v2 = Vt;
	v = Vt + 4;
	for (unsigned int z = 0; z != nz; ++z) {
		const GRD_data_type **F0 = *FG;
		const GRD_data_type **F1 = *(++FG);
		for (unsigned int y = 0; y != ny; ++y) {
			const GRD_data_type *V00 = *F0;
			const GRD_data_type *V01 = *(++F0);
			const GRD_data_type *V10 = *F1;
			const GRD_data_type *V11 = *(++F1);
			v2[0] = iso - *V00;//the difference was inverted to use the signbit function
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			//the eight least significant bits of i correspond to vertex indices. (x...x01234567)
			//If the bit is 1 then the vertex value is greater than zero.
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[0])) i |= 8;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[2])) i |= 2;
			for (unsigned int x = 0; x != Nx; ++x) {
				swap(v2, v);
				V00 += d;
				v2[0] = iso - *V00;//*(++V00);
				V01 += d;
				v2[1] = iso - *V01;//*(++V01);
				V11 += d;
				v2[2] = iso - *V11;//*(++V11);
				V10 += d;
				v2[3] = iso - *V10;//*(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[0])) i |= 8;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[2])) i |= 2;
				if (i && i^0xFF) { //i is different from 0 and 0xFF
					if (v2 < v) {MC33_real *t = v2; MC33_real *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}// memcpy(v2 + 8, v2, 4*sizeof(MC33_real));
					find_case(x, y, z, i);
				}
			}
		}
		swap(Dx, Ux);
		swap(Dy, Uy);
	}
	try {
		Sf.color = vector<int>(Sf.nV, DefaultColor);
	}
	catch (...) {
		memoryfault = 1;
	}
	if (memoryfault) {
		Sf.clear();
		return -1;
	}
	return 0;
}

size_t MC33::size_of_isosurface(MC33_real iso, unsigned int &nV, unsigned int &nT) {
	const GRD_data_type ***FG = F;
	if (!FG)//The set_grid3d function was not executed
		return -2;
	surface Sf;
	S = &Sf;
	Sf.iso = iso;
	memoryfault = 0;
	unsigned int d = di, Nx = nx;
	MC33_real Vt[12];
	MC33_real *v2 = Vt;
	v = Vt + 4;
	for (unsigned int z = 0; z != nz; ++z) {
		const GRD_data_type **F0 = *FG;
		const GRD_data_type **F1 = *(++FG);
		for (unsigned int y = 0; y != ny; ++y) {
			const GRD_data_type *V00 = *F0;
			const GRD_data_type *V01 = *(++F0);
			const GRD_data_type *V10 = *F1;
			const GRD_data_type *V11 = *(++F1);
			v2[0] = iso - *V00;
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[2])) i |= 2;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[0])) i |= 8;
			for (unsigned int x = 0; x != Nx; ++x) {
				swap(v, v2);
				V00 += d;
				v2[0] = iso - *V00;//*(++V00);
				V01 += d;
				v2[1] = iso - *V01;//*(++V01);
				V11 += d;
				v2[2] = iso - *V11;//*(++V11);
				V10 += d;
				v2[3] = iso - *V10;//*(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[0])) i |= 8;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[2])) i |= 2;
				if (i && i^0xFF) {
					if (v > v2) {MC33_real *t = v2; MC33_real *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}
					case_count(x, y, z, i);
				}
			}
		}
		swap(Dx, Ux);
		swap(Dy, Uy);
	}
	nV = Sf.nV;
	nT = Sf.nT;
	// number of vertices * (size of vertex and normal + size of color ) + number of triangle * size of triangle + size of class surface
	return nV*(3*(sizeof(MC33_real) + sizeof(float)) + sizeof(int)) + nT*(3*sizeof(int)) + sizeof(surface);
}

size_t MC33::size_of_isosurface(MC33_real iso) {
	unsigned int nV, nT;
	return size_of_isosurface(iso, nV, nT);
}

#ifdef _MSC_VER
#pragma warning( pop )
#endif

