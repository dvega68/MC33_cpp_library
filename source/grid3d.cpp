/*
	File: grid3d.cpp
	Programmed by: David Vega - dvega@uc.edu.ve
	August 2019
	August 2020
	June 2021
	August 2021
	December 2021
*/

#ifdef _MSC_VER
#pragma warning( disable : 4244 )
#endif

#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

#include "../include/MC33.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

#ifndef GRD_orthogonal
void setIdentMat3x3d(double (*A)[3]) {
	for (double *d = A[0] + 8; --d != A[0];)
		d[0] = 0.0;
	for (int i = 0; i != 3; ++i)
		A[i][i] = 1.0;
}
#endif

//******************************************************************
const unsigned int* grid3d::get_N() {return N;}
const float* grid3d::get_L() {return L;}
const double* grid3d::get_r0() {return r0;}
const double* grid3d::get_d() {return d;}
const char* grid3d::get_title() {return title;}

GRD_data_type grid3d::get_grid_value(unsigned int i, unsigned int j, unsigned int k) {
	if (F && i <= N[0] && j <= N[1] && k <= N[2])
		return F[k][j][x_data > 1? i*x_data: i];
	return sqrt(-1); // NaN value
}

GRD_data_type grid3d::interpolated_value(double x, double y, double z) {
	double r[3] = {x - r0[0], y - r0[1], z - r0[2]};
#ifndef GRD_orthogonal
	if(nonortho)
		multAb(A_,r,r,0);
#endif
	return (this->*interpolation)(r);
}

GRD_data_type grid3d::bad_value(double*) const {
	return sqrt(-1); // NaN value
}

int grid3d::set_interpolation(int i) {
	if (!F || x_data > 0)
		i = 0;
	switch (i) {
		case 1:
			interpolation = &grid3d::trilinear;
			break;
		case 3:
			interpolation = &grid3d::tricubic;
			break;
		default:
			i = 0;
			interpolation = &grid3d::bad_value;
			break;
	}
	return (i? 0: -1);
}


GRD_data_type grid3d::trilinear(double *r) const {
	unsigned int i[3];
	double t[3], s[3];
	for (int j = 0; j < 3; j++) {
		if ((periodic>>j)&0x01) {
			t[j] = (r[j] - L[j]*floor(r[j]/L[j]))/d[j];
			i[j] = (int)t[j];
			t[j] -= i[j];
			i[j] = (i[j]? i[j] - 1: N[j] - 1);
		} else {
			t[j] = r[j]/d[j];
			if (t[j] < 0 || t[j] > N[j])
				return sqrt(-1);
			i[j] = (int)t[j];
			t[j] -= i[j];
			if (i[j] == N[j]) {
				t[j] = 1.0;
				--i[j];
			}
		}
		s[j] = 1.0 - t[j];
	}
	return s[0]*s[1]*s[2]*F[i[2]][i[1]][i[0]] + t[0]*s[1]*s[2]*F[i[2]][i[1]][i[0] + 1] + s[0]*t[1]*s[2]*F[i[2]][i[1] + 1][i[0]] + s[0]*s[1]*t[2]*F[i[2] + 1][i[1]][i[0]] + t[0]*t[1]*s[2]*F[i[2]][i[1] + 1][i[0] + 1] + s[0]*t[1]*t[2]*F[i[2] + 1][i[1] + 1][i[0]] + t[0]*s[1]*t[2]*F[i[2] + 1][i[1]][i[0] + 1] + t[0]*t[1]*t[2]*F[i[2] + 1][i[1] + 1][i[0] + 1];
}

//from https://facyt-quimicomp.neocities.org/Vega_en.html#c_library
GRD_data_type grid3d::tricubic(double *r) const {
	unsigned int iR[3], x, y, z, i, j, k;
	double f = 0, t, w[3][4];
	for (i = 0; i < 3; i++)
	{
		if ((periodic>>i)&0x01)
		{
			t = (r[i] - L[i]*floor(r[i]/L[i]))/d[i];
			iR[i] = (int)t;
			t -= iR[i];
			iR[i] = (iR[i]? iR[i] - 1: N[i] - 1);
		}
		else
		{
			t = r[i]/d[i];
			if (t < 0 || t > N[i])
				return sqrt(-1);
			iR[i] = (int)t;
			t -= iR[i];
			if (!iR[i])
				t--;
			else if (iR[i] > N[i] - 2)
			{
				z = iR[i] - N[i] + 2;
				iR[i] -= z + 1;
				t += z;
			}
			else
				iR[i]--;
		}
		w[i][0] = (-2.0 + (3.0 - t)*t)*t*(1.0/6.0);
		w[i][1] = (2.0 + (-1.0 + (-2.0 + t)*t)*t)*0.5;
		w[i][2] = (2.0 + (1.0 - t)*t)*t*0.5;
		w[i][3] = (t*t - 1.0)*t*(1.0/6.0);
	}
	z = iR[2];
	for (k = 0; k < 4; k++)
	{
		y = iR[1];
		for (j = 0; j < 4; j++)
		{
			x = iR[0];
			for (i = 0; i < 4; i++)
			{
				f += w[0][i]*w[1][j]*w[2][k]*F[z][y][x];
				if ((++x) > N[0]) x = 1;
			}
			if ((++y) > N[1]) y = 1;
		}
		if ((++z) > N[2]) z = 1;
	}
	return f;
}


#ifndef GRD_orthogonal
const float* grid3d::get_Ang() {return Ang;}
const double (*grid3d::get__A())[3] {return _A;}
const double (*grid3d::get_A_())[3] {return A_;}
int grid3d::isnotorthogonal() {return nonortho;}

void grid3d::update_matrices() {
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90) {
		nonortho = 1;
		double ca = cos(Ang[0]*(M_PI/180.0));
		double cb = cos(Ang[1]*(M_PI/180.0));
		double aux1 = Ang[2]*(M_PI/180.0);
		double sg = sin(aux1);
		double cg = cos(aux1);
		aux1 = ca - cb*cg;
		double aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		_A[0][0] = A_[0][0] = 1.0;
		_A[0][1] = cg;
		_A[0][2] = cb;
		_A[1][1] = sg;
		A_[1][1] = cb = 1.0/sg;
		A_[0][1] = -cg*cb;
		_A[1][2] = aux1*cb;
		_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		A_[1][2] = -aux1*cb*aux2;
		A_[2][2] = sg*aux2;
		_A[1][0] = _A[2][0] = _A[2][1] = 0.0;
		A_[1][0] = A_[2][0] = A_[2][1] = 0.0;
		multAbf = mult_TSAbf;
	} else {
		nonortho = 0;
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
}

void grid3d::set_Ang(float angle_bc, float angle_ca, float angle_ab) {
	Ang[0] = angle_bc; Ang[1] = angle_ca; Ang[2] = angle_ab;
	update_matrices();
}
#endif

void grid3d::set_ratio_aspect(double rx, double ry, double rz) {
	if (F) {
		d[0] = rx; d[1] = ry; d[2] = rz;
		for (int i = 0; i != 3; ++i)
			L[i] = N[i]*d[i];
	}
}

void grid3d::set_r0(double x, double y, double z) {
	r0[0] = x; r0[1] = y; r0[2] = z;
}

void grid3d::set_title(const char *s) {
	char *t = title, *te = title + sizeof(title) - 1;
	while (*s && t < te)
		*(t++) = *(s++);
	*t = (char)0;
}

void grid3d::set_grid_value(unsigned int i, unsigned int j, unsigned int k, GRD_data_type value) {
	if (!x_data && F && i <= N[0] && j <= N[1] && k <= N[2])
		F[k][j][i] = value;
}

//******************************************************************
void grid3d::free_F() {
	if (F) {
		clear_subgrid();
		if (x_data)
			for (unsigned int k = 0; k <= N[2]; ++k)
				delete[] F[k];
		else
			for (unsigned int k = 0; k <= N[2]; ++k) {
				if (F[k]) {
					for (unsigned int j = 0; j <= N[1]; ++j)
						delete[] F[k][j];
				} else {
					k = N[2];
					break;
				}
				delete[] F[k];
			}
		delete[] F;
		F = 0;
		interpolation = &grid3d::bad_value;
	}
}

grid3d::grid3d() : F(0), subgrid(0), nsg(0), maxnsg(0) {}

grid3d::grid3d(const grid3d &G) {
	if (!G.F)
		return;
	memcpy(this, &G, sizeof (grid3d));
	F = 0;
	subgrid = 0;
	nsg = maxnsg = 0;
	if (alloc_F())
		return;
	for (unsigned int k = 0; k <= N[2]; ++k)
		for (unsigned int j = 0; j <= N[1]; ++j) {
			if (G.x_data > 1)
				for (unsigned int i = 0; i <= N[0]; ++i)
					F[k][j][i] = G.F[k][j][i*G.x_data];
			else
				memcpy(F[k][j], G.F[k][j], (N[0] + 1)*sizeof(GRD_data_type));
		}
}

grid3d::~grid3d() {
	free_F();
}

int grid3d::alloc_F() {
	unsigned int j, k;
	F = new (nothrow) GRD_data_type**[N[2] + 1];
	if (!F)
		return -1;
	for (k = 0; k <= N[2]; ++k) {
		F[k] = new (nothrow) GRD_data_type*[N[1] + 1];
		if (!F[k])
			return -1;
		for (j = 0; j <= N[1]; ++j) {
			F[k][j] = new (nothrow) GRD_data_type[N[0] + 1];
			if (!F[k][j]) {
				while (j)
					delete[] F[k][--j];
				delete[] F[k];
				F[k] = 0;
				return -1;
			}
		}
	}
	x_data = 0;
	interpolation = &grid3d::trilinear;
	return 0;
}

void grid3d::delete_grid_data() {
	free_F();
	N[0] = N[1] = N[2] = 0;
}

int grid3d::set_grid_dimensions(unsigned int Nx, unsigned int Ny, unsigned int Nz) {
	free_F();
	if (Nx == 0 || Ny == 0 || Nz == 0) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (int i = 0; i != 3; ++i) {
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_orthogonal
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return 0;
}

int grid3d::set_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data) {
	free_F();
	if (!data || Nx == 0 || Ny == 0 || Nz == 0) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	F = new (nothrow) GRD_data_type**[Nz];
	if (!F)
		return -2;
	for (unsigned int k = 0; k < Nz; ++k) {
		F[k] = new (nothrow) GRD_data_type*[Ny];
		if (!F[k]) {
			while (k)
				delete[] F[--k];
			delete[] F;
			F = 0;
			return -2;
		}
		for (unsigned int j = 0; j < Ny; ++j)
			F[k][j] = data + j*Nx;
		data += Ny*Nx;
	}
	for (int i = 0; i != 3; ++i) {
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_orthogonal
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	x_data = -1;
	return 0;
}

void grid3d::clear_subgrid() {
	if (subgrid) {
		for (unsigned int i = 0; i != nsg; i++)
			delete subgrid[i];
		nsg = maxnsg = 0;
		delete[] subgrid;
		subgrid = 0;
	}
}

grid3d *grid3d::get_subgrid(unsigned int i) {
	return (i < nsg? subgrid[i]: 0);
}

void grid3d::del_subgrid(unsigned int i) {
	if (i < nsg) {
		delete subgrid[i];
		while (++i < nsg)
			subgrid[i - 1] = subgrid[i];
		--nsg;
	}
}

unsigned int grid3d::subgrid_size() {
	return nsg;
}

int grid3d::add_subgrid(unsigned int Oi, unsigned int Oj, unsigned int Ok,
					unsigned int Ni, unsigned int Nj, unsigned int Nk,
					unsigned int Si, unsigned int Sj, unsigned int Sk) {
	if (!F || x_data > 0 || !Ni || !Nj || !Nk || !Si || !Sj || !Sk)
		return -1;
	if (Oi + Si*(Ni - 1) > N[0] || Oj + Sj*(Nj - 1) > N[1] || Ok + Sk*(Nk - 1) > N[2])
		return -1;
	if (nsg == maxnsg) {
			grid3d **t = new (nothrow) grid3d*[maxnsg + 8];
			if (!t)
				return -2;
			if (subgrid) {
				memcpy(t, subgrid, nsg*sizeof(void*));
				delete[] subgrid;
		}
		subgrid = t;
		maxnsg += 8;
	}
	grid3d *G = new (nothrow) grid3d();
	if (!G)
		return -2;

	G->N[0] = Ni - 1;
	G->N[1] = Nj - 1;
	G->N[2] = Nk - 1;
	G->F = new (nothrow) GRD_data_type**[Nk];
	if (!G->F) {
		delete G;
		return -2;
	}
	for (unsigned int k = 0; k < Nk; ++k) {
		G->F[k] = new (nothrow) GRD_data_type*[Nj];
		if (!G->F[k]) {
			while (++k < Nk)
				G->F[k] = 0;
			delete G;
			return -2;
		}
		for (unsigned int j = 0; j < Nj; ++j)
			G->F[k][j] = F[Ok + k*Sk][Oj + j*Sj] + Oi;
	}
	G->d[0] = Si*d[0];
	G->d[1] = Sj*d[1];
	G->d[2] = Sk*d[2];
	G->title[0] = 0;
	double d0[3] = {Oi*d[0], Oj*d[1], Ok*d[2]};
#ifndef GRD_orthogonal
	memcpy(G->A_, A_, sizeof A_);
	memcpy(G->_A, _A, sizeof _A);
	memcpy(G->Ang, Ang, sizeof Ang);
	G->nonortho = nonortho;
	if (nonortho)
		multAb(_A, d0, d0, 0);
#endif
	for (int i = 0; i != 3; i++) {
		G->r0[i] = r0[i] + d0[i];
		G->L[i] = G->N[i]*G->d[i];
	}

	G->x_data = Si;
	G->periodic = 0;
	if (Si*Ni == N[0] + 1)
		G->periodic = periodic&1;
	if (Sj*Nj == N[1] + 1)
		G->periodic |= periodic&2;
	if (Sk*Nk == N[2] + 1)
		G->periodic |= periodic&4;
	subgrid[nsg++] = G;
	interpolation = &grid3d::bad_value;
	return 0;
}

int grid3d::generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf, double dx, double dy, double dz, double (*fn)(double x, double y, double z)) {
	free_F();
	if (dx <= 0 || dy <= 0 || dz <= 0 || xi == xf || yi == yf || zi == zf) {
		N[0] = N[1] = N[2] = 0;
		return -1;
	}
	if (xi > xf)
		swap(xi, xf);
	if (xf - xi < dx)
		dx = xf - xi;
	if (yi > yf)
		swap(yi, yf);
	if (yf - yi < dy)
		dy = yf - yi;
	if (zi > zf)
		swap(zi, zf);
	if (zf - zi < dz)
		dz = zf - zi;
	N[0] = int((xf - xi)/dx + 0.5);
	N[1] = int((yf - yi)/dy + 0.5);
	N[2] = int((zf - zi)/dz + 0.5);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	d[0] = dx; d[1] = dy; d[2] = dz;
	r0[0] = xi; r0[1] = yi; r0[2] = zi;
	if (fn) {
		double x, y, z = zi;
		for (size_t k = 0; k <= N[2]; ++k) {
			y = yi;
			for (size_t j = 0; j <= N[1]; ++j) {
				x = xi;
				for (size_t i = 0; i <= N[0]; ++i) {
					F[k][j][i] = (GRD_data_type)fn(x,y,z);
					x += dx;
				}
				y += dy;
			}
			z += dz;
		}
	}
	for (int i = 0; i != 3; ++i) {
		L[i] = N[i]*d[i];
#ifndef GRD_orthogonal
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return 0;
}


//******************************************************************
/*
read_grd reads a filename file (the file must be a output *.grd file from the
DMol program), it returns a pointer to struct _GRD that contains all the grid
data.
*/
int grid3d::read_grd(const char *filename) {
	char cs[32];
	unsigned int i, j, k;
	int xi[3], order;

	free_F();
	ifstream in(filename);
	if (!in)
		return -1;
	nonortho = 0;
	periodic = 0;
	in.getline(title, 159);
	in.ignore(60, '\n');
#ifdef GRD_orthogonal
	float Ang[3];
#endif
	in >> L[0] >> L[1] >> L[2] >> Ang[0] >> Ang[1] >> Ang[2] >> N[0] >> N[1] >> N[2];
	in >> order >> xi[0] >> cs >> xi[1] >> cs >> xi[2];
	in.ignore(20, '\n');
	if (N[0] < 1 || N[1] < 1 || N[2] < 1) return -1;
	if (order != 1 && order != 3) return -1;
	for (i = 0; i != 3; ++i) {
		d[i] = L[i]/N[i];
		r0[i] = xi[i]*d[i];
	}
	periodic = (xi[0] == 0)|((xi[1] == 0)<<1)|((xi[2] == 0)<<2);

#ifndef GRD_orthogonal
	update_matrices();
#endif

	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (k = 0; k <= N[2]; ++k)
		if (order == 1) {
			for (j = 0; j <= N[1]; ++j)
				for (i = 0; i <= N[0]; ++i)
				{
					in.getline(cs, 31);
					F[k][j][i] = stof(cs);
				}
		} else {
			for (i = 0; i <= N[0]; ++i)
				for (j = 0; j <= N[1]; ++j)
				{
					in.getline(cs, 31);
					F[k][j][i] = stof(cs);
				}
		}
	return (in.good()? 0: -4);
}

/*
internal binary format
*/
int grid3d::read_grd_binary(const char* filename) {
	unsigned int i;
	free_F();
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i != 0x4452475f) // _GRD
		return -1;

	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i > 159)
		return -1;
	in.read(title, i);
	in.read(reinterpret_cast<char*>(N), sizeof N);
	in.read(reinterpret_cast<char*>(L), sizeof L);
	in.read(reinterpret_cast<char*>(r0), sizeof r0);
	in.read(reinterpret_cast<char*>(d), sizeof d);
#ifndef GRD_orthogonal
	in.read(reinterpret_cast<char*>(&nonortho), sizeof(int));
	if (nonortho) {
		in.read(reinterpret_cast<char*>(Ang),3*sizeof(float));
		in.read(reinterpret_cast<char*>(_A), sizeof _A);
		in.read(reinterpret_cast<char*>(A_), sizeof A_);
		multAbf = mult_Abf;
	} else {
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
#else
	in.read(reinterpret_cast<char*>(&i), sizeof(int));
	if (i)
		in.ignore(3*sizeof(float) + 18*sizeof(double));
#endif
	periodic = (r0[0] == 0 && r0[1] == 0 && r0[2] == 0? 7: 0);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	for (unsigned int k = 0; k <= N[2]; ++k)
		for (unsigned int j = 0; j <= N[1]; ++j)
#if GRD_type_size == 8
			for (unsigned int i = 0; i <= N[0]; ++i) {
				float f;
				in.read(reinterpret_cast<char*>(&f), sizeof(float));
				F[k][j][i] = f;
			}
#else
			in.read(reinterpret_cast<char*>(F[k][j]),(N[0] + 1)*sizeof(float));
#endif
	return (in.good()? 0: -4);
}
//******************************************************************

/*
Reads a set of files that contain a slab of res*res scan data points, the data
points are read as unsigned short int (if order is different from 0, the bytes
of the unsigned short are exchanged). The filename must end with a number, and
the fuction read all files with end number greater or equal to filename.
(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)
*/
int grid3d::read_scanfiles(const char *filename, unsigned int res, int order) {
	string nm(filename);
	ifstream in;
	unsigned int i, j, l;
	int m, k = -1;
	unsigned short int n;
	free_F();
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (i = 0; i != 2; ++i) {
		d[i] = 1.0;
		N[i] = res - 1;
		L[i] = float(res - 1);
	}
	d[2] = 1.0;
	l = unsigned(nm.size() - 1);
	while (nm[l] >= '0' && nm[l] <= '9')
		l--;
	m = stoi(nm.substr(++l));
	while (1) {
		nm.replace(l,6,to_string(m++));
		in.open(nm, ios::binary|ifstream::in);
		if (!in) {
			m = 0;
			break;
		}
		if (!((++k)&63)) {
			GRD_data_type ***pt = new (nothrow) GRD_data_type**[k + 64];
			if (!pt) {
				m = -2;
				break;
			}
			if (F) {
				memcpy(pt, F, k*sizeof(void*));
				delete[] F;
			}
			F = pt;
		}
		F[k] = new (nothrow) GRD_data_type*[res];
		if (!F[k]) {
			m = -2;
			break;
		}
		for (j = 0; j != res; ++j)
			F[k][j] = new (nothrow) GRD_data_type[res];
		if (!F[k][N[1]]) {
			m = -2;
			break;
		}

		if (order)
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i) {
					in.read(reinterpret_cast<char*>(&n), sizeof(short int));
					F[k][j][i] = static_cast<unsigned short int>((n>>8)|(n<<8));
				}
		else
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i) {
					in.read(reinterpret_cast<char*>(&n), sizeof(short int));
					F[k][j][i] = n;
				}
		if (in.fail()) {
			for (j = 0; j != res; ++j)
				delete[] F[k][j];
			delete[] F[k];
			m = (--k > 0? 0: -4);//m = -4;
			break;
		}
		in.close();
	}
	N[2] = k;
	L[2] = float(k);
	x_data = 0;
	interpolation = &grid3d::trilinear;
	if (m == -2)
		free_F();
	else {
		j = (k + 1)/2;
		for (i = 0; i != j; ++i)
			swap(F[i], F[k - i]);
	}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return m;
}

//******************************************************************
/*
read_raw_file reads a file that contains integer (8, 16 or 32 bits) or float (or double)
data points. byte is the number of bytes of the integer (1, 2 or 4), or of the float (4
or 8). If the data is big endian, byte must be negative (only for integers). The vector
n[3] contains the number of points in each dimension. The size of file must be
abs(byte)*n[0]*n[1]*n[2]. The function returns zero when succeeds.
*/
int grid3d::read_raw_file(const char *filename, unsigned int *n, int byte, int isfloat) {
	unsigned int i, j, k;
	free_F();
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
	unsigned int ui = 0;
	if (isfloat) {
		if (byte != sizeof(float) && byte != sizeof(double))
			return -1;
	} else if (abs(byte) > 4 || abs(byte) == 3 || !byte)
		return -1;
	if (byte == -1) byte = 1;
	if (n[0] == 0 || n[1] == 0 || n[2] == 0)
		return -1; // bad input parameters
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h) {
		d[h] = 1.0;
		N[h] = n[h] - 1;
		L[h] = float(N[h]);
	}
	if (alloc_F()) {
		free_F();
		return -2;
	}
#ifdef integer_GRD
	if (!isfloat && GRD_type_size == byte)
#else
	if (isfloat && GRD_type_size == byte)
#endif
	{
		byte *= n[0];
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				in.read(reinterpret_cast<char*>(F[k][j]), byte);
	} else if (isfloat) {
		if (byte == 8) {
#if defined (integer_GRD) || GRD_type_size == 4
			double df;
			for (k = 0; k != n[2]; ++k)
				for (j = 0; j != n[1]; ++j)
					for (i = 0; i != n[0]; ++i)
					{
						in.read(reinterpret_cast<char*>(&df), byte);
						F[k][j][i] = GRD_data_type(df);
					}
#endif
		} else {
#if defined (integer_GRD) || GRD_type_size == 8
			float f;
			for (k = 0; k != n[2]; ++k)
				for (j = 0; j != n[1]; ++j)
					for (i = 0; i != n[0]; ++i) {
						in.read(reinterpret_cast<char*>(&f), byte);
						F[k][j][i] = GRD_data_type(f);
					}
#endif
		}
	} else if (byte < 0) {
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				for (i = 0; i != n[0]; ++i) {
					in.read(reinterpret_cast<char*>(&ui), -byte);
					if (byte == -2)
						F[k][j][i] = GRD_data_type((ui>>8)|(ui<<8));
					else if (byte == -4)
						F[k][j][i] = GRD_data_type((ui>>24)|((ui>>8)&0x00f0)|((ui<<8)&0x0f00)|(ui<<24));
					else
						F[k][j][i] = GRD_data_type((ui>>16)|(ui&0x00f0)|(ui<<16));
				}
	} else {
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				for (i = 0; i != n[0]; ++i) {
					in.read(reinterpret_cast<char*>(&ui), byte);
					F[k][j][i] = GRD_data_type(ui);
				}
	}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
/*
Reads a dat file:
http://www.cg.tuwien.ac.at/research/vis/datasets/
*/
int grid3d::read_dat_file(const char *filename) {
	unsigned short int n, nx, ny, nz;
	free_F();
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h)
		d[h] = 1.0;
	in.read(reinterpret_cast<char*>(&nx), sizeof(short int));
	in.read(reinterpret_cast<char*>(&ny), sizeof(short int));
	in.read(reinterpret_cast<char*>(&nz), sizeof(short int));
	N[0] = nx - 1;
	N[1] = ny - 1;
	N[2] = nz - 1;
	for (int h = 0; h != 3; ++h)
		L[h] = float(N[h]);
	if (alloc_F()) {
		free_F();
		return -2;
	}
	while (nz--)
		for (unsigned int j = 0; j < ny; ++j)
			for (unsigned int i = 0; i < nx; ++i) {
				in.read(reinterpret_cast<char*>(&n), sizeof(short int));
				F[nz][j][i] = n;
			}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
/*
Save all point values of the grid in binary format
*/
int grid3d::save_raw_file(const char *filename) {
	if (!F)
		return -1;
	ofstream out(filename, ios::binary);
	if (!out)
		return -1;
	for (unsigned int k = 0; k <= N[2]; ++k)
		for (unsigned int j = 0; j <= N[1]; ++j) {
			if (x_data > 1)
				for (unsigned int i = 0; i <= N[0]; ++i)
					out.write(reinterpret_cast<char*>(F[k][j] + i*x_data), sizeof(GRD_data_type));
			else
				out.write(reinterpret_cast<char*>(F[k][j]), (N[0] + 1)*sizeof(GRD_data_type));
		}
	return (out.good()? 0: -4);
}

