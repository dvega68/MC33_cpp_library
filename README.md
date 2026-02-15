### Feel free to use the MC33++ library.

---

#### INFO:

MC33++ library version 5.3

This library is a C++ version based on the MC33 library of the paper:  
Vega, D., Abache, J., Coll, D., [A Fast and Memory-Saving Marching Cubes 33 implementation with the correct interior test](http://jcgt.org/published/0008/03/01), *Journal of Computer Graphics Techniques (JCGT)*, vol. 8, no. 3, 1-18, 2019.

The MC33 library is an open source software. The distribution and use rights are under the terms of the [MIT license](https://opensource.org/licenses/MIT), described in the file "LICENSE.txt".

![FLTK example](https://repository-images.githubusercontent.com/469939412/decb05bb-c9dc-4019-96bc-11f1f6dee8c8 "Screenshot of the FLTK example")

---

#### FILES:

- Makefile (Linux or MinGW/msys GCC makefile)
- MakefileMSVC.mak (NMAKE makefile)
- compileMSVC.bat (batch script to compile with visual c++)
- include/MC33.h (header file)
- source/MC33.cpp (`MC33` class code)
- source/grid3d.cpp (`grid3d` class code)
- source/surface.cpp (`surface` class code)
- source/MC33_LookUpTable.h (Triangulation pattern for each MC33 case)
- source/libMC33++.cpp (source file used to compile the library)
- FLTK_example/TestMC33.cpp (Example of use. FLTK library is required)
- FLTK_example/makefileMinGW-w64.mak (MinGW/msys GCC makefile)
- FLTK_example/makefiledebian.mak (Debian GCC makefile)
- GLUT_example/TestMC33_glut.cpp (Example of use. GLUT or FREEGLUT library is required)
- GLUT_example/makefileMinGW-w64.mak (MinGW/msys GCC makefile)
- GLUT_example/makefiledebian.mak (Debian GCC makefile)

The visualc folder contains the solution and project files to compile the library with Visual Studio 2019 Community.

---

#### CUSTOMIZING:

There are 4 options that can be modified before compiling the library. You can do it by editing the MC33.h or libMC33++.cpp file before compiling the library:

1. To change the data type of the grid (the default value is float) define GRD_TYPE_SIZE and/or GRD_INTEGER (MC33.h). For example:
	```c
	#define GRD_TYPE_SIZE 8 // the data type is double

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 4 // the data type is unsigned int

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 2 // the data type is unsigned short int

	#define GRD_INTEGER
	#define GRD_TYPE_SIZE 1 // the data type is unsigned char
	```

2. If you do not use inclined grids, you can define GRD_ORTHOGONAL (MC33.h):
	```c
	#define GRD_ORTHOGONAL
	```

3. By default the members of MC33 class are float. The member can be changed to double by defining MC33_DOUBLE_PRECISION to 1 (MC33.h). This option is also enabled by defining GRD_TYPE_SIZE to 8. And it is disabled if GRD_TYPE_SIZE is defined to 1 or 2 when GRD_INTEGER is defined. The vertex array type of `surface` class is also modified when using this option.
	```c
	#define MC33_DOUBLE_PRECISION 1 // double type for MC33 class members
	```

4. If you need to exchange the front and back surfaces, define MC33_NORMAL_NEG (libMC33++.cpp):
	```c
	#define MC33_NORMAL_NEG
	```

---

#### INSTALLING:

1. Compile the libMC33++.cpp file as a static C++ library:
	- A GCC makefile is supplied with the library. In a Linux terminal or in msys2 mingw console go to the folder where the Makefile file is, and type: make
	- If you are using Visual Studio, open the visualc/MC33++.sln file with the Visual Studio IDE, select the appropriated configuration (eg. release x64) and build the project. If you don't want to open the Visual Studio IDE, you can run the compileMSVC.bat script from file explorer.

	Once the library is compiled copy the libMC33++.a (or MC33++.lib) from the local lib directory to the compiler lib directory. Also copy the MC33.h file from the local include directory to the compiler include directory.

	Include the header file in your C++ code:
	```c
	#include <MC33.h>
	```
	and put in the linker options of your program makefile: -lMC33++


2. Instead of compiling the library, you can directly include the library code files in your code (see the FLTK or the GLUT example code). Put at the beginning of your C++ code:
	```c
	#include "..Path../source/MC33.cpp"
	#include "..Path../source/grid3d.cpp"
	#include "..Path../source/surface.cpp"
	```

---

#### COMPILING THE EXAMPLES:

In Debian terminal window, go to the FLTK_example or GLUT_example folder and write:
```sh
make -f makefiledebian.mak
```

Or in a msys2 MinGW64 Shell (Windows), write:
```sh
make -f makefileMinGW-w64.mak
```

For the FLTK example in any operating system you also can use the fltk-config script:
```sh
path/fltk-1.X.Y/fltk-config --use-gl --compile TestMC33.cpp
```

The makefiles use the -Ofast optimization option and the fltk-config script uses a lower optimization level.

In the GLUT example, the file containing the grid must be passed to the program on the command line, and no other grid files can be read from the running program. The grid file can be dragged and dropped into the executable in the Windows File Explorer. Examples of usage of the `generate_grid_from_fn` function of the `grid3d` class were included in this code, and are available if the grid file is not specified.

In the FLTK example, a new grid can be read from the running program, but all previous surfaces will be removed from memory. This code has keyboard shortcuts similar to the GLUT example, but using ctrl instead of alt. No usage examples for the `generate_grid_from_fn` function were included, but subgrid management was included.

---

#### USAGE THE LIBRARY IN YOUR CODE:

1. Create a `grid3d` object, and read a data file (use the member functions `read_grd`, `read_grd_binary`, `read_scanfiles`, `read_raw_file` or `read_dat_file`):
	```c
	  grid3d G;
	  G.read_dat_file("filename.dat");

	  // or as a pointer

	  grid3d *Z = new grid3d;
	  Z->read_dat_file("filename.dat");
	```

2. Create an `MC33` object and assign it the `grid3d` object:
	```c
	  MC33 MC;
	  MC.set_grid3d(G); // or MC.set_grid3d(Z);
	```

3. Calculate an isosurface with the MC33 algorithm:
	```c
	  surface S;
	  MC.calculate_isosurface(S, isovalue);
	```

See MC33.h file for the use of other functions.

---

#### GRIDS

The `grid3d` class has functions for building and managing subgrids that use the same data as the main grid. For example, the resolution can be reduced, and the generated isosurfaces will use less memory:
```c
  grid3d G;
  G.read_dat_file("filename.dat");
  const unsigned int *N = G.get_N();
  // build a subgrid with the half of resolution in each dimension:
  G.add_subgrid(0, 0, 0, N[0]/2, N[1]/2, N[2]/2, 2, 2, 2);

  MC33 MC;
  // The first subgrid uses one eighth of the data from the main grid
  MC.set_grid3d(G.get_subgrid(0));
```

By modifying the parameters of `add_subgrid` the grid can also be split.  
The subgrids can be deleted by using `del_subgrid(i)`, where `i` is the subgrid index.

There is another way to create a `grid3d` object. The `grid3d::generate_grid_from_fn` function permits build a grid by using a scalar function `double fn(double x, double y, double z)`.

for example:
```c
// sphere function
double fs(double x, double y, double z) {
  const double radius = 1.0;
  const double cx = 2.0, cy = 2.0, cz = 2.0;
  x -= cx; y -= cy; z -= cz;
  return radius*radius - x*x - y*y - z*z;
}

  .
  .
  grid3d G;
  G.generate_grid_from_fn(0.5, 0.5, 0.5, // coordinates of the grid origin
                          3.5, 3.5, 3.5, // coordinates of the opposite corner
                          0.03, 0.03, 0.03, // steps
                          fs);
  MC33 MC;
  MC.set_grid3d(G);
  .
  .
```

If fn (the last argument of `generate_grid_from_fn`) is NULL, an empty grid will be created but with memory reserved for the data. The data can be filled using the `set_grid_value` function.

If you already have a data array of the same type as the data in the `grid3d` class, you can use the `set_data_pointer` function to set the internal pointers to the grid data. This avoids duplicating the data. When the `grid3d` object is destroyed, the external data will not be modified.

This library contains interpolation functions (trilinear and tricubic type so far). The default interpolation function is the trilinear type. It can be changed to tricubic using the `set_interpolation function`. The `interpolated_value(x, y, z)` function is used to get the interpolated value at the x, y, z position.

```c
  grid3d G;
  G.generate_grid_from_fn(0.5, 0.5, 0.5, 3.5, 3.5, 3.5, 0.03, 0.03, 0.03, fs);
  std::cout.precision(9);
  std::cout << "\nSphere: "  << fs(2.3, 2.3, 2.3); // value at a grid point
  std::cout << " " << fs(2.301, 2.302, 2.301); // point that is not on the grid
  //G.set_interpolation(1);
  std::cout << "\nLinear: " << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
  G.set_interpolation(3);
  std::cout << "\n Cubic: "  << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
  G.set_interpolation(0); // no valid value
  std::cout << "\n   Nan: "  << G.interpolated_value(2.3, 2.3, 2.3);
  std::cout << " " << G.interpolated_value(2.301, 2.302, 2.301);
```

For more information, see the `grid3d` class in the MC33.h file.

---

#### OTHERS:

To display the surface, you can use the `draw()` and `drawdfaft()` functions of the `surface` class by defining `MC33_USE_DRAW_OPEN_GL` (in only one of the project files) before including the MC33.h file in your code. These functions use the OpenGL library. Alternatively, you can implement your own drawing functions in your code after including MC33.h.

To calculate the size (in bytes) of an isosurface, without calculating the isosurface, use:
```c
  size_t size = MC.size_of_isosurface(iso, nV, nT);
```
where iso is the isovalue (a `float` or `double`), nV and nT are unsigned integers that will contain the number of vertices and triangles, respectively. If you do not want to calculate nV and nT, use:
```c
  size_t size = MC.size_of_isosurface(iso);
```
See [this link](https://stackoverflow.com/questions/65066235/estimating-size-of-marching-cubes-output-geometry)

Two new funtions where added to save the surface: `surface::save_obj` and `surface::save_ply`, the first saves the surface data in a Wavefront .obj file, the other saves the data in a "Polygon File Format" (.ply) file.

---

The MC33.h file contains a description of all the functions of this library.

---

#### ACKNOWLEDGEMENT:

Thanks to Miguel D&iacute;az (\*), for testing the library.  
And thanks again to Julien Hess (\*\*) for his helpful suggestions.

(\*) Miguel D&iacute;az e-mail: <mdiaz92@outlook.com>, IVIC doctoral student.  
(\*\*) Julien Hess e-mail: <julien.hess.ch@gmail.com>, consultant of [MATHICSE-Group](https://www.epfl.ch/labs/mathicse/)

---

See [MC33_libraries](https://facyt-quimicomp.neocities.org/MC33_libraries.html) web page.  
Mail to: <dvega@uc.edu.ve>
