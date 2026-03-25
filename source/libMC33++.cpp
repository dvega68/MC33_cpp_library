/*
	File: libMC33++.cpp
	Programmed by: David Vega
	August 2019
	August 2020
	April 2021
	December 2021
	February 2026
	March 2026
	Only this file must be included in a project to compile the MC33 library
*/

/*****************************************************************************
You can change the following lines before compiling the library: */
#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c /* RGBA 0xAABBGGRR: red 92 green 92 blue 92 (grey) */
#endif
#ifndef MC33_NORMAL_NEG
#define MC33_NORMAL_NEG 0 /* If it is 1, the front and back surfaces are exchanged. */
#endif
#ifndef USE_INTERNAL_SIGNBIT
#define USE_INTERNAL_SIGNBIT 1 /* definition of signbf function, see MC33.cpp */
#endif
#ifndef USE_MM_RSQRT_SS
#define USE_MM_RSQRT_SS 1 /* definition of invSqrt function, see MC33.cpp */
#endif
/*****************************************************************************/

#include "../include/MC33.h"

#define compiling_libMC33
#include "grid3d.cpp"
#include "surface.cpp"
#include "MC33.cpp"

