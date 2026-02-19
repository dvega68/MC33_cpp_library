/*
	File: MC33.h
	Programmed by: David Vega - dvega@uc.edu.ve
	version: 5.3
	March 2022
	February 2026
	This library is the C++ version of the library described in the paper:
	Vega, D., Abache, J., Coll, D., A Fast and Memory Saving Marching Cubes 33
	implementation with the correct interior test, Journal of Computer Graphics
	Techniques (JCGT), vol. 8, no. 3, 1–17, 2019.
*/

#ifndef MC33_h_
#define MC33_h_

/********************************** USAGE ************************************/
/*
//1. Header
#include <MC33.h>

//2. Read a grid file.
	grid3d G;
	G.read_dat_file("filename.dat");

//3. create a MC33 object and assign it the grid3d.
	MC33 MC;
	MC.set_grid3d(G);
//4. calculate an isosurface.
	surface S;
	MC.calculate_isosurface(S, isovalue);

*/

/********************************CUSTOMIZING**********************************/
//The following defines can be only changed before compiling the library:
//#define GRD_INTEGER // for dataset with integer type
#define GRD_TYPE_SIZE 4 // 1, 2, 4 or 8 (8 for double, if not defined GRD_INTEGER)
#define MC33_DOUBLE_PRECISION 0 // 1 means double type for MC33 class members, used only with double or size 8 integer grid data
//#define GRD_ORTHOGONAL // If defined, the library only works with orthogonal grids.
//#define MC33_NORMAL_NEG // the front and back surfaces are exchanged.
//#define DEFAULT_SURFACE_COLOR 0xFF80FF40// RGBA 0xAABBGGRR: red 64, green 255, blue 128
/*****************************************************************************/

#define MC33_VERSION_MAJOR 5
#define MC33_VERSION_MINOR 3

#ifdef GRD_INTEGER
#if GRD_TYPE_SIZE == 4
typedef unsigned int GRD_data_type; // variable type of the grid data, by default it is float.
#elif GRD_TYPE_SIZE == 2
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 0
typedef unsigned short int GRD_data_type;
#elif GRD_TYPE_SIZE == 1
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 0
typedef unsigned char GRD_data_type;
#else
#error "Incorrect size of the integer data type. GRD_TYPE_SIZE permitted values: 1, 2 or 4."
#endif
#elif GRD_TYPE_SIZE == 8
#undef MC33_DOUBLE_PRECISION
#define MC33_DOUBLE_PRECISION 1
typedef double GRD_data_type;
#elif MC33_DOUBLE_PRECISION
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 8
typedef double GRD_data_type;
#else
typedef float GRD_data_type;
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 4
#endif

#if MC33_DOUBLE_PRECISION
typedef double MC33_real;
#else
typedef float MC33_real;
#endif

#include <vector>
#include <functional>
class MC33;

/*
The class grid3d contains a 3D matrix (F[][][]) that stores the values of a function
evaluated at points of a 3D regularly spaced grid. N[] is the number of intervals in
each dimension. The grid contains (N[2] + 1)*(N[1] + 1)*(N[0] + 1) points. L[] is the
grid size in each dimension. r0[] are the coordinates of the first grid point. d[]
is the distance between adjacent points in each dimension (can be different for each
dimension), nonortho has a value of 1 when the grid is inclined else 0. _A and A_
are the matrices that transform from inclined to orthogonal coordinates and vice
versa, respectively. If the grid is periodic (is infinitely repeated along each
dimension) the flag periodic must be different from 0.

In this library, if GRD_ORTHOGONAL is defined, then nonortho, _A and A_ can be
removed from this structure, and it only works with orthogonal grids.
*/
class grid3d {
private:
	GRD_data_type ***F;
	int x_data;
	// if x_data is negative, the data of grid are external an they cannot be modified using
	// the member functions of this class. If it is positive, the grid data is managed by the
	// functions of this class. In the subgrids x_data is equal to the inner index step.
	int periodic; // to work with periodical grids.
	double r0[3], d[3]; // grid origin coordinates, distances between grid points
	unsigned int N[3];
	float L[3]; // size of the grid axes
	char title[160];
	grid3d **subgrid;
	unsigned int nsg, maxnsg; // for subgrids
	GRD_data_type (grid3d::*interpolation)(double*) const; // pointer to interpolation function
#ifndef GRD_ORTHOGONAL
	int nonortho;
	float Ang[3]; // angles between grid axes.
	double _A[3][3], A_[3][3];
	void update_matrices(); // set _A and A_ by using Ang
#endif
	int alloc_F(); // allocates memory for the grid data
	void free_F(); // release the allocated memory
	GRD_data_type bad_value(double*) const; // always returns nan
//Interpolation functions:
	GRD_data_type trilinear(double *r) const;
	GRD_data_type tricubic(double *r) const;
public:
//Get pointers to some class members:
	const unsigned int* get_N();
	const float* get_L();
	const double* get_r0();
	const double* get_d();
	const char* get_title();
#ifndef GRD_ORTHOGONAL
	const float* get_Ang();
	const double (*get__A())[3];
	const double (*get_A_())[3];
	int isnotorthogonal(); // returns nonortho value
#endif

//Generates an orthogonal grid from a function fn(x,y,z). xi and xf are the
//limits of the interval along the x axis, yi and yf along the y axis and zi
//and zf along the z axis. dx, dy and dz are the respective step sizes.
	int generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf,
		double dx, double dy, double dz, double (*fn)(double x, double y, double z));
//Get a grid point value
	GRD_data_type get_grid_value(unsigned int i, unsigned int j, unsigned int k);
//Calculate a value at position (x, y, z) by interpolating the grid values.
//This function don't work with subgrids.
	GRD_data_type interpolated_value(double x, double y, double z);
//Select the interpolation type, 1 for trilinear and 3 for tricubic.
	int set_interpolation(int i);

//Modifying the grid parameters:
	/******************************************************************
	set_grid_dimensions sets the new dimensions of a grid data. It overrides the effect of
	the other functions that modify the grid parameters. Nx, Ny and Nz are the number of
	grid points in each dimension. It returns 0 if memory allocation was successful.*/
	int set_grid_dimensions(unsigned int Nx, unsigned int Ny, unsigned int Nz);
	void set_grid_value(unsigned int i, unsigned int j, unsigned int k, GRD_data_type value); // set a grid point value
	void set_ratio_aspect(double rx, double ry, double rz); // modifies d and L
	void set_r0(double x, double y, double z); // modifies r0
#ifndef GRD_ORTHOGONAL
	void set_Ang(float angle_bc, float angle_ca, float angle_ab); // modifies Ang
#endif
	void set_title(const char *s); // copy the c style string s to title
	void delete_grid_data(); // Delete the internal grid data

	/******************************************************************
	set_data_pointer creates internal pointers that point to the external data array. data
	must be stored with the nested inner loop running from i = 0 to Nx - 1 and the outer
	loop from k = 0 to Nz - 1. The data content cannot be modified using class functions.
	It returns 0 if memory allocation of internal pointers was successful.*/
	int set_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data);

// Managing subgrids:
	/******************************************************************
	The grid can contain several subgrids that use the same data. The add_subgrid function
	make a subgrid and stores it in an internal array (the subgrids do not contain data,
	only pointers to the main grid data). The input parameters of the function are the three
	indices of the origin grid point, the number of points in each dimension and the index
	step in each dimension. The function returns 0 if the subgrid was successfully added,
	-1 for wrong input parameters and -2 for error in memory allocation. */
	int add_subgrid(unsigned int Oi, unsigned int Oj, unsigned int Ok,
					unsigned int Ni, unsigned int Nj, unsigned int Nk,
					unsigned int Si, unsigned int Sj, unsigned int Sk);
	void clear_subgrid(); // erase all subgrids
	grid3d *get_subgrid(unsigned int i); // returns a pointer to the subgrid i
	void del_subgrid(unsigned int i); // delete the subgrid i
	unsigned int subgrid_size(); // returns the number of subgrids

//Reading grid data from files:
	/******************************************************************
	The functions returns zero when succeeds. If the file could not be opened or the data
	does not match the format, the return value is -1. If a memory error occurred, the
	return value is -2. A return value of -4 means that the data read may be incomplete. */
	/******************************************************************
	read_grd reads a *.grd file from the DMol3 program.*/
	int read_grd(const char *filename);

	/** read_grd_binary reads a file with an internal binary format */
	int read_grd_binary(const char *filename);
	/******************************************************************
	read_scanfiles reads a set of files that contain a slab of res*res scan data points,
	the data points are read as unsigned short int (if order is different from 0, the
	bytes of the unsigned short are exchanged). The filename must end with a number, and
	the function reads all files with end number greater or equal to filename.
	(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)*/
	int read_scanfiles(const char *filename, unsigned int res, int order);

	/******************************************************************
	read_raw_file reads a file that contains integer (8, 16 or 32 bits) or float (or double)
	data points. byte is the number of bytes of the integer (1, 2 or 4), or of the float (4
	or 8). If the data is big endian, byte must be negative (only for integers). The vector
	n[3] contains the number of points in each dimension. The size of file must be
	abs(byte)*n[0]*n[1]*n[2].*/
	int read_raw_file(const char *filename, unsigned int *n, int byte, int isfloat = 0);

	/******************************************************************
	read_dat_file reads a dat file. The function returns 0 when succeeds.
	http://www.cg.tuwien.ac.at/research/vis/datasets/*/
	int read_dat_file(const char *filename);

	/******************************************************************
	save the grid points values in a raw data file of MC33_real (float or double). The
	function returns 0 when succeeds.*/
	int save_raw_file(const char *filename);

	grid3d();
	// Copy constructor
	grid3d(const grid3d &);
	~grid3d();
friend MC33;
};

template <class T>
struct MC33_v3 {
	T v[3];
};

/* The class surface contains the data of a isosurface. The isovalue is iso, the
number of surface points is nV, and the number of triangles is nT. The vector V
contains the vertex coordinates, the vector T contains the triangle indices, The
vector N contains the normal coordinates (one normal for each vertex), and the
color vector contains the color index of each point.*/
class surface {
private:
	unsigned int nV, nT;
	std::vector<MC33_v3<unsigned int>> T;
	std::vector<MC33_v3<MC33_real>> V;
	std::vector<MC33_v3<float>> N;
	std::vector<int> color;
	MC33_real iso;
	unsigned int sflag;
public:
	union {
		void *p;
		long long ul;
		int i[2];
		short si[4];
		char c[8];
		float f[2];
		double df;
	} user; // user data
	MC33_real get_isovalue(); // returns the isovalue
	unsigned int get_num_vertices(); // gets the number of vertices
	unsigned int get_num_triangles(); // gets the number of triangles
	const unsigned int *getTriangle(unsigned int n); // gets a pointer to indices of triangle n
	const MC33_real *getVertex(unsigned int n); // gets a pointer to coordinates of vertex n
	const float *getNormal(unsigned int n); // gets a pointer to the normal vector n
	void flipNormals(); // reverses the direction of all normals
	void flipTriangles(); // toggle CW / CCW vertex order in triangles
	const unsigned char *getColor(unsigned int n); // gets a pointer to the color of vertex n
	void setColor(unsigned int n, unsigned char *pcolor);

	/******************************************************************
	Saves all the surface *S data (in binary format) to a "filename" file. The
	return value is 0 if the call succeeds, else -1.*/
	int save_bin(const char *filename);

	/******************************************************************
	Saves all the surface *S data (in plain text format) to a "filename" file.
	The return value is 0 if the call succeeds, else -1.*/
	int save_txt(const char *filename);

	/******************************************************************
	Saves the surface *S data (without the color) to Wavefront .obj file.
	The return value is 0 if the call succeeds, else -1.*/
	int save_obj(const char *filename);

	/******************************************************************
	Saves the surface *S data to Polygon File Format .ply file.
	https://paulbourke.net/dataformats/ply/
	The return value is 0 if the call succeeds, else -1.*/
	int save_ply(const char *filename, const char* author = 0, const char* object = 0);

	/******************************************************************
	Reads (from a "filename" file) the surface data stored in binary format.
	The return value is 0 if the call succeeds, else -1.*/
	int read_bin(const char *filename);

	/* Draw the surface, this function can be implemented by the user.*/
	void draw();

	/* Draw some points of the surface, this function can be implemented by the user.*/
	void drawdraft();

	/* Clear all vector data */
	void clear();

	/* Correct the data vector sizes */
	void adjustvectorlenght();

	surface();
friend MC33;
};

/* Marching cubes 33 class.
The function member set_grid3d must be called once before calculate an isosurface.
If the grid3d object is modified or you want to use another grid3d object with the
same MC33 object, this function must be called again.
The function calculate_isosurface fill the surface object with the isosurface data.
*/
class MC33 {
private:
	static int DefaultColor;
	surface *S;
	//Auxiliary grid variables
	unsigned int nx, ny, nz;
	const GRD_data_type ***F;
	MC33_real MC_O[3], MC_D[3], ca, cb;
#ifndef GRD_ORTHOGONAL
	double _A[3][3], A_[3][3];
#endif
	/*Assign memory for the vertex r[3], normal (r + 3)[3]. The return value is
	the new vertex label.*/
	std::function<unsigned int(MC33_real*)> store_point;

	//Other auxiliary variables
	int memoryfault;
	unsigned int di; // for subgrids, index step for inner loop
	// temporary structures that store the indexes of triangle vertices:
	unsigned int **Dx, **Dy, **Ux, **Uy, **Lz;
	MC33_real *v;
	const unsigned short int table[2310]; // Triangle pattern look up table
	//Procedures
	int face_tests(int *, int) const;
	int face_test1(int) const;
	int interior_test(int, int) const;
	unsigned int surfint(unsigned int, unsigned int, unsigned int, MC33_real *);
	void find_case(unsigned int, unsigned int, unsigned int, unsigned int);
	void case_count(unsigned int, unsigned int, unsigned int, unsigned int);
	int init_temp_isosurface();
	void free_temp_D_U();
	void clear_temp_isosurface();
public:
	// Set the color of the next isosurface
	void set_default_surface_color(unsigned char *color);
	// set the grid parameters:
	int set_grid3d(grid3d *G);
	int set_grid3d(grid3d &G);
	// Calculate the isosurface with isovalue iso and store the data in the surface Sf:
	int calculate_isosurface(surface &Sf, MC33_real iso);
	/* Return the size in bytes of an isosurface with out calculate it (nV and nT are
	the number of vertices and triangles):*/
	std::size_t size_of_isosurface(MC33_real iso, unsigned int &nV, unsigned int &nT);
	// Return the size in bytes of an isosurface with out calculate it:
	std::size_t size_of_isosurface(MC33_real iso);
	MC33();
	~MC33();
};

#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c; // grey red 92 green 92 blue 92
#endif

#ifndef compiling_libMC33

/* Define MC33_USE_DRAW_OPEN_GL before include MC33.h (only once in your project), to
compile the following surface::draw and surface::drawdraft functions. */
#ifdef MC33_USE_DRAW_OPEN_GL
#if MC33_DOUBLE_PRECISION
#define GL_MC33_real GL_DOUBLE
#else
#define GL_MC33_real GL_FLOAT
#endif
void surface::draw() {
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_MC33_real, 0, &V[0]);
	glNormalPointer(GL_FLOAT, 0, &N[0]);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, &color[0]);//rgb alpha

	glDrawElements(GL_TRIANGLES, 3*nT, GL_UNSIGNED_INT, &T[0]);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void surface::drawdraft() {
	glDisable(GL_LIGHTING);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_MC33_real, 12*sizeof(MC33_real), &V[0]);
	glColorPointer(3, GL_UNSIGNED_BYTE, 16, &color[0]);//rgb

	glDrawArrays(GL_POINTS, 0, nV>>2);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_LIGHTING);
}
#endif // MC33_USE_DRAW_OPEN_GL

#ifndef GRD_ORTHOGONAL
//c = Ab, A is a 3x3 upper triangular matrix. If t != 0, A is transposed.
template<typename T> void T_multTSA_b(const double (*A)[3], T *b, T *c, int t) {
	if (t) {
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
		c[1] = A[0][1]*b[0] + A[1][1]*b[1];
		c[0] = A[0][0]*b[0];
	} else {
		c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		c[1] = A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][2]*b[2];
	}
}
//Performs the multiplication of the matrix A and the vector b: c = Ab. If t != 0, A is transposed.
template<typename T> void T_multA_b(const double (*A)[3], T *b, T *c, int t) {
	double u,v;
	if (t) {
		u = A[0][0]*b[0] + A[1][0]*b[1] + A[2][0]*b[2];
		v = A[0][1]*b[0] + A[1][1]*b[1] + A[2][1]*b[2];
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
	} else {
		u = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		v = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
	}
	c[0] = u;
	c[1] = v;
}

extern void (*multAbf)(const double (*)[3], MC33_real *, MC33_real *, int);
extern void (*mult_TSAbf)(const double (*)[3], MC33_real *, MC33_real *, int);
extern void (*mult_Abf)(const double (*)[3], MC33_real *, MC33_real *, int);
#endif // GRD_ORTHOGONAL

#endif // compiling_libMC33

#endif // MC33_h_

