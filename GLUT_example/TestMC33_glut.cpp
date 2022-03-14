/*
	File: TestMC33_glut.cpp
	Programmed by: David Vega - dvega@uc.edu.ve
	June 2021
	December 2021
	This is an open source code. The distribution and use rights are under the terms of the MIT license (https://opensource.org/licenses/MIT)
*/

#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>

#ifdef __APPLE__
/* Defined before OpenGL and GLUT includes to avoid deprecation messages */
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <io.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

/*
#include <MC33.h> // put -lMC33++ in the linker libraries
/*/
//#define GRD_type_size 8
#include "../source/MC33.cpp"
#include "../source/grid3d.cpp"
#include "../source/surface.cpp"
//*/

using namespace std;

int WIDTH = 600, HEIGHT = 450;

list <surface> ls;
list <surface>::iterator cS;
grid3d G;
MC33 MC;

surface *surf;
float bgr, bgg, bgb;
float lightPosition[2][4];
float MGL[16];
float oldM[9], iniV[3];
float tX, tY;
int oldX, oldY;
float scaleSize, scaleFactor;
float cX, cY, cZ;

char surfSizeStr[80], timeStr[32];

string isoString, auxInputString;
bool showMessage, inputIso, saveIso, showHelp, setSurfColor;
float txtF, txtX, txtY;

bool drawPoints;
int mouseButton;
int mouseButtonRotate, mouseButtonMove, mouseButtonScale;

void map_to_sphere(float x, float y, float *NewV) {
	float length = x*x + y*y;
	if (length > 1.0f) {
		float norm = 1.0f/sqrtf(length);
		NewV[0] = x * norm;
		NewV[1] = y * norm;
		NewV[2] = 0.0f;
	} else {
		NewV[0] = x;
		NewV[1] = y;
		NewV[2] = sqrtf(1.0f - length);
	}
}

void reset_view() {
	float *f;
	for (f = MGL + 12; --f != MGL;) // MGL[0] is assigned later
		*f = 0.0f;
	for (f = oldM + 8; --f != oldM;) // oldM[0] and oldM[8] are assigned later
		*f = 0.0f;
	oldM[0] = oldM[4] = oldM[8] = MGL[0] = MGL[5] = MGL[10] = scaleFactor = 1.0f;
	tX = tY = 0.0f;
	MGL[15] = scaleSize;
	drawPoints = false;
	glutPostRedisplay();
}

void reshape(int w, int h) {
	if (w < 120 || h < 90) {
		if (w < 120)
			w = 120;
		if (h < 90)
			h = 90;
		glutReshapeWindow(w, h);
	}
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w < h) {
		float f = (float)h/w;
		txtF = 8.0f/w;
		txtX = 2*txtF - 1.0f;
		txtY = f - txtF*(26.0f/8);
		glOrtho(-1.0f, 1.0f, -f, f, -10.0f, 10.0f);
	} else {
		float f = (float)w/h;
		txtF = 8.0f/h;
		txtX = 2*txtF - f;
		txtY = 1.0f - txtF*(26.0f/8);
		glOrtho(-f, f, -1.0f, 1.0f, -10.0f, 10.0f);
	}
	glMatrixMode(GL_MODELVIEW);
	WIDTH = w;
	HEIGHT = h;
	glutPostRedisplay(); // update the drawing
}

void set_glMaterial(GLenum face) {
	float mat_specular[4] = { 1.0f, 0.7f, 0.2f, 1.0f};
	float shininess = 50.0f;

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, face == GL_FRONT_AND_BACK);
	glColorMaterial(face, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(face, GL_SPECULAR, mat_specular);
	glMaterialf(face, GL_SHININESS, shininess);
}

void init_gl() {
	float light0Ambient[4] = {0.2f, 0.2f, 0.2f, 1.0f};
	float light0Diffuse[4] = {0.8f, 0.8f, 0.8f, 1.0f};

	float light1Ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
	float light1Diffuse[4] = {0.6f, 0.6f, 0.6f, 1.0f};

	// light 0 settings
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0Ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0Diffuse);

	// light 1 settings
	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1Ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1Diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1Diffuse);

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	set_glMaterial(GL_FRONT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(bgr, bgg, bgb, 0.0);
}

void init_param() {
	mouseButtonRotate = GLUT_LEFT_BUTTON;
	mouseButtonMove = GLUT_RIGHT_BUTTON;
	mouseButtonScale = GLUT_MIDDLE_BUTTON;
	tX = tY = bgr = bgg = bgb = cX = cY = cZ = 0.0f;
	scaleFactor = scaleSize = 1.0f;
	lightPosition[0][0] = 1.5f;
	lightPosition[0][1] = -0.5f;
	lightPosition[0][2] = 2.0f;
	lightPosition[0][3] = 0.0f;
	lightPosition[1][0] = -1.5f;
	lightPosition[1][1] = 0.5f;
	lightPosition[1][2] = 2.0f;
	lightPosition[1][3] = 0.0f;
	setSurfColor = inputIso = saveIso = false;
	showHelp = true;
	init_gl();
}

// display a char string in the glut window
void displaytext(const char *s, float x, float y, bool pos = true) {
	if (pos)
		glRasterPos2f(x, y);
	while (*s)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(s++));
}

void info_text() {
	glLoadIdentity();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	if (showMessage || saveIso) {
		glColor4f(0.0f, 0.0f, 1.0f, 0.8f);
		float halfw = (saveIso || setSurfColor? 48.0f: 32.0f)*txtF;
		glEnable(GL_BLEND);
		glRectf(-halfw, -16.0f*txtF, halfw, 16.0f*txtF);
		glDisable(GL_BLEND);
		glColor3f(0.75f, 0.75f, 0.75f);
		displaytext("Press Esc to close this box", txtF - halfw, -15.0f*txtF);
		if (saveIso) {
			if (showMessage) {
				displaytext("Type the surface filename (*.txt or *.sup):", -43*txtF, (20.0f/8)*txtF);
				glColor3f(1.0f, 1.0f, 1.0f);
				displaytext(auxInputString.c_str(), -40*txtF, (-20.0f/8)*txtF);
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
			} else
				displaytext("File already exists, overwrite? (y/N)", -37*txtF, 0.0f);
		} else if (setSurfColor) {
			displaytext("Input r (0-255) g (0-255) b (0-255) colors:", -43*txtF, (20.0f/8)*txtF);
			glColor3f(1.0f, 1.0f, 1.0f);
			displaytext(auxInputString.c_str(), -11*txtF, (-20.0f/8)*txtF);
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
		} else
			displaytext(auxInputString.c_str(), -txtF*auxInputString.size(), 0.0f);
	} else if (WIDTH > 328 && HEIGHT > 225) {
		if (showHelp) {
			glColor4f(0.20f, 0.20f, 0.20f, 0.6f);
			glEnable(GL_BLEND);
			glRectf(-45.0f*txtF, -27.5f*txtF, 45.0f*txtF, 28.5f*txtF);
			glDisable(GL_BLEND);
			glColor3f(1.0f, 1.0f, 1.0f);
			float f = -41*txtF;
			displaytext("KEYBOARD MENU:", -30*txtF, (90.0f/4)*txtF);
			displaytext("I            : Enter an isovalue", f, (70.0f/4)*txtF);
			displaytext("H            : Show/hide this information", f, (55.0f/4)*txtF);
			displaytext("S            : Save the current surface", f, (40.0f/4)*txtF);
			displaytext("C            : Sets the surface color", f, (25.0f/4)*txtF);
			displaytext("R            : Resets the view", f, (10.0f/4)*txtF);
			displaytext("F            : Swap back and front faces", f, -(5.0f/4)*txtF);
			displaytext("L            : back face lighting on/off", f, -(20.0f/4)*txtF);
			displaytext("Escape key   : Hide help and messages", f, (-35.0f/4)*txtF);
			displaytext("Space bar    : Show the next surface", f, (-50.0f/4)*txtF);
			displaytext("Alt + Space  : Show the previous surface", f, (-65.0f/4)*txtF);
			displaytext("Alt + Delete : Delete the current surface", f, (-80.0f/4)*txtF);
			displaytext("Alt + Q      : Exits the program", f, (-95.0f/4)*txtF);
		} else {
			glColor3f(0.75f, 0.75f, 0.75f);
			displaytext("Press h for help", -32*txtF - txtX , txtY);
		}
	}
		
	glColor3f(0.75f, 0.75f, 0.75f);
	displaytext("Isovalue: ", txtX, txtY);
	if (inputIso) {
		glColor3f(1.0f, 1.0f, 1.0f);
		displaytext(auxInputString.c_str(), txtX + 20*txtF, txtY);
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '_');
		glColor3f(0.75f, 0.75f, 0.75f);
	} else
		displaytext(surf? isoString.c_str(): "", 0, 0, false);
	displaytext(surfSizeStr, txtX, (-14.0/8)*txtF - txtY);
	if (8*(int)(strlen(timeStr) + strlen(surfSizeStr)) + 16 <= WIDTH)
		displaytext(timeStr, -txtX - 2*txtF*strlen(timeStr), (-14.0/8)*txtF - txtY);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}

// render function
void draw() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition[0]);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition[1]);
	MGL[12] = tX*scaleSize - cX*MGL[0] - cY*MGL[4] - cZ*MGL[8];
	MGL[13] = tY*scaleSize - cX*MGL[1] - cY*MGL[5] - cZ*MGL[9];
	MGL[14] = -cX*MGL[2] - cY*MGL[6] - cZ*MGL[10];
	//MGL[15] = scaleFactor*scaleSize;
	glLoadMatrixf(MGL);

	if (surf) {
		if (drawPoints)
			surf->drawdraft();
		else
			surf->draw();
	}

	info_text();
	glutSwapBuffers();
}

// mouse motion
void motion(int x, int y) {
	if (mouseButton == mouseButtonRotate) {
		float qm[9]; // the 3 first elements of qm are the end vector to map to sphere
		// the elements 3 to 6 of qm are the rotation quaternion.
		// To simplify the calculation, the angle of rotation is multiplied by 2.
		map_to_sphere(2.0f*x/(WIDTH - 1) - 1.0f,1.0f - 2.0f*y/(HEIGHT - 1), qm);
		qm[4] = iniV[2]*qm[1] - iniV[1]*qm[2]; // qi
		qm[5] = iniV[0]*qm[2] - iniV[2]*qm[0]; // qj
		qm[6] = iniV[1]*qm[0] - iniV[0]*qm[1]; // qk
		float s = qm[4]*qm[4] + qm[5]*qm[5] + qm[6]*qm[6]; // sin^2(angle)
		if (s < 0.00001f*0.00001f)
			return;
		qm[3] = iniV[0]*qm[0] + iniV[1]*qm[1] + iniV[2]*qm[2]; // qr, cos(angle) instead cos(angle/2)
		s += qm[3]*qm[3];
		s = (s > 0.0f ? 2.0f/s : 0.0f);
		// for notation see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
		float qi2 = qm[4]*s, qj2 = qm[5]*s, qk2 = qm[6]*s;
		float qiqr = qi2*qm[3], qjqr = qj2*qm[3], qkqr = qk2*qm[3];
		float qiqj = qm[4]*qj2, qiqk = qm[4]*qk2, qjqk = qm[5]*qk2;
		qi2 *= qm[4];
		qj2 *= qm[5];
		qk2 *= qm[6];
		// now qm is the rotation matrix:
		qm[0] = 1.0f - (qj2 + qk2); qm[1] = qiqj - qkqr; qm[2] = qiqk + qjqr;
		qm[3] = qiqj + qkqr; qm[4] = 1.0f - (qi2 + qk2); qm[5] = qjqk - qiqr;
		qm[6] = qiqk - qjqr; qm[7] = qjqk + qiqr; qm[8] = 1.0f - (qi2 + qj2);
		//Update the OpenGL transformation matrix:
		MGL[0] = qm[0]*oldM[0] + qm[3]*oldM[1] + qm[6]*oldM[2];
		MGL[4] = qm[0]*oldM[3] + qm[3]*oldM[4] + qm[6]*oldM[5];
		MGL[8] = qm[0]*oldM[6] + qm[3]*oldM[7] + qm[6]*oldM[8];
		MGL[1] = qm[1]*oldM[0] + qm[4]*oldM[1] + qm[7]*oldM[2];
		MGL[5] = qm[1]*oldM[3] + qm[4]*oldM[4] + qm[7]*oldM[5];
		MGL[9] = qm[1]*oldM[6] + qm[4]*oldM[7] + qm[7]*oldM[8];
		MGL[2] = qm[2]*oldM[0] + qm[5]*oldM[1] + qm[8]*oldM[2];
		MGL[6] = qm[2]*oldM[3] + qm[5]*oldM[4] + qm[8]*oldM[5];
		MGL[10] = qm[2]*oldM[6] + qm[5]*oldM[7] + qm[8]*oldM[8];
	} else if (mouseButton == mouseButtonMove) {
		float t = 2*(WIDTH > HEIGHT ? scaleFactor/HEIGHT : scaleFactor/WIDTH);
		tX += (x - oldX)*t;
		tY -= (y - oldY)*t;
		oldX = x;
		oldY = y;
	} else if (mouseButton == mouseButtonScale) {
		scaleFactor += (y - oldY)*scaleFactor*0.001f;
		MGL[15] = scaleFactor*scaleSize;
		oldY = y;
	}
	glutPostRedisplay();
}

void show_message(const char *msg) {
	showMessage = true;
	auxInputString = string(msg);
	glutMouseFunc(0);
	glutMotionFunc(0);
	glutPostRedisplay();
}

void mouse(int, int, int, int);
void hide_message() {
	showMessage = false;
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
}

void fill_surface_info() {
	if (ls.empty()) {
		surfSizeStr[0] = 0;
		timeStr[0] = 0;
		surf = 0;
	} else {
		ostringstream info;
		surf = &(*cS);
		int i = distance(ls.begin(), cS);
		info << "Surface: " << i + 1 << "|" << ls.size() << "  Vertices: ";
		info << surf->get_num_vertices() << "  Triangles: " << surf->get_num_triangles();
		strcpy(surfSizeStr, info.str().c_str());
		info.str("");
		info << "Time: " << surf->user.i[0] << " ms";
		strcpy(timeStr, info.str().c_str());
		isoString = to_string(surf->get_isovalue());
		while (isoString.back() == '0') {
			isoString.pop_back();
		}
		if (isoString.back() == '.')
			isoString.pop_back();
	}
	glutPostRedisplay();
}

bool file_exist(const char *filename) {
	return (access(filename, F_OK)? false: true);
}

void delete_iso() {
	if (ls.empty())
		return;
	auto t = cS;
	if (ls.size() > 1) {
		if (cS == ls.begin())
			++cS;
		else
			--cS;
	}
	ls.erase(t);
	fill_surface_info();
}

void calc_iso(double iso) {
	for (auto it = ls.begin(); it != ls.end(); it++)
		if (it->get_isovalue() == iso) {
			cS = it;
			fill_surface_info();
			return;
		}
	surface S;
	int t = clock();
	switch (MC.calculate_isosurface(S, iso)) {
	case -2:
		show_message("No grid data!");
		break;
	case -1:
		show_message("Insufficient memory!");
		break;
	default:
		t = clock() - t;
		if (S.get_num_vertices()) {
			S.user.i[0] = t*1000/CLOCKS_PER_SEC;
			ls.push_back(move(S));
			cS = prev(ls.end());
			fill_surface_info();
		} else
			show_message("Empty surface!");
	}
}

void select_iso(bool next) {
	if (ls.empty())
		return;
	if (next) {
		if (++cS == ls.end())
			cS = ls.begin();
	} else {
		if (cS == ls.begin())
			cS = ls.end();
		--cS;
	}
	fill_surface_info();
}

void save_iso2() {
	const char *c = strrchr(auxInputString.c_str(), '.');
	if (c) {
		if (strcmp(c,".sup"))
			c = 0;
	}
	if (c)
		surf->save_bin(auxInputString.c_str());
	else
		surf->save_txt(auxInputString.c_str());
}

// called when a mouse button is pressed
void mouse(int button, int state, int x, int y) {
	if (state == GLUT_DOWN) { // if any button is pressed
		oldX = x;
		oldY = y;
		if (button == mouseButtonRotate)
			map_to_sphere(2.0f*oldX/(WIDTH - 1) - 1.0f,1.0f - 2.0f*oldY/(HEIGHT - 1), iniV);
		mouseButton = button;
		drawPoints = true;
	} else {
		memcpy(oldM, MGL, 3*sizeof(float));
		memcpy(oldM + 3, MGL + 4, 3*sizeof(float));
		memcpy(oldM + 6, MGL + 8, 3*sizeof(float));
		drawPoints = false;
		glutPostRedisplay();
	}
}

void keyboard(unsigned char c, int, int);

void keyb_string(unsigned char c, int, int) {
	if (c >= ' ' && c <= '}')
		auxInputString += static_cast<char>(c);
	else if (c == '\b' && auxInputString.size())
		auxInputString.pop_back();
	else if (c == '\r' || c == 27) {
		if (inputIso) {
			if (c == '\r') {
				double iso;
				try {
					iso = stod(auxInputString);
				} catch (exception&) {
					auxInputString.clear();
					glutPostRedisplay();
					return;
				}
				isoString = auxInputString;
				calc_iso(iso);
			}
			inputIso = false;
		} else if (saveIso) {
			hide_message();
			if (auxInputString.empty())
				saveIso = false;
			else if (!file_exist(auxInputString.c_str())) {
				saveIso = false;
				save_iso2();
			}
		} else if (setSurfColor) {
			setSurfColor = false;
			hide_message();
			if (auxInputString.size() > 4) {
				stringstream ss(auxInputString);
				int rgb[3];
				unsigned char c[4] = {0, 0, 0, 0xff};
				ss >> rgb[0] >> rgb[1] >> rgb[2];
				for (int i = 0; i != 3; i++)
					c[i] = rgb[i]&0xff;
				MC.set_default_surface_color(c);
				if (ls.size()) {
					int i = surf->get_num_vertices();
					while(i)
						surf->setColor(--i, c);
				}
			}
		}
		glutKeyboardFunc(keyboard);
	}
	glutPostRedisplay();
}

void iso_input() {
	auxInputString.clear();
	glutKeyboardFunc(keyb_string);
	inputIso = true;
	hide_message();
}

void set_color() {
	auxInputString.clear();
	glutKeyboardFunc(keyb_string);
	showMessage = setSurfColor = true;
}

void save_iso() {
	if (ls.empty())
		return;
	auxInputString.clear();
	saveIso = showMessage = true;
	glutKeyboardFunc(keyb_string);
}

void keyboard(unsigned char key, int, int) {
	if (saveIso) {
		if (key == 'Y' || key == 'y')
			save_iso2();
		saveIso = false;
	} else
		switch(key) {
		case 27: // escape key
			showHelp = false;
			hide_message();
			break;
		case ' ':
			select_iso(GLUT_ACTIVE_ALT != glutGetModifiers());
			break;
		case '\b': // Backspace key
		case 127: // Delete key
			if (GLUT_ACTIVE_ALT == glutGetModifiers()) // if the Control key is pressed
				delete_iso();
			break;
		case 'C':
		case 'c':
			set_color();
		break;
		case 'F':
		case 'f': // swap back and front faces
			if (ls.size()) {
				int v;
				glGetIntegerv(GL_FRONT_FACE, &v);
				glFrontFace(v == GL_CW? GL_CCW: GL_CW);
				surf->flipNormals();
				glutPostRedisplay();
			}
			break;
		case 'H':
		case 'h':
			if (showMessage) {
				showHelp = true;
				hide_message();
			} else
				showHelp = !showHelp;
			break;
		case 'I':
		case 'i':
			iso_input();
			break;
		case 'L':
		case 'l': // enable/disable back face lighting
			{
				GLboolean v;
				glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE, &v);
				set_glMaterial(v? GL_FRONT: GL_FRONT_AND_BACK);
			}
			break;
		case 'Q':
		case 'q':
			if (GLUT_ACTIVE_ALT == glutGetModifiers()) // if the Control key is pressed
				exit(0);
			break;
		case 'R':
		case 'r':
			reset_view();
			break;
		case 'S':
		case 's':
			save_iso();
			break;
		}
	glutPostRedisplay();
}

string terminal_input(const char *msg) {
	string is;
	cout << msg;
	cin.clear();
	getline(cin, is);
	return is;
}

bool input_res(int *n, long long size) {
	if (size < 8)
		return false;
	string s = terminal_input("Grid dimensions (Nx Ny Nz): ");
	if (s.size() > 4) {
		stringstream ss(s);
		ss >> n[0] >> n[1] >> n[2];
	} else
		return false;
	long long div = n[0]*n[1]*n[2];
	if (!div)
		return false;
	n[5] = n[3] = size/div;
	if (n[3] == 4)
		s = terminal_input("Data type r (real) or i (integer): ");
	else {
		s[0] = (n[3] < 4? 'i': 'r');
	}
	switch (s[0]) {
	case 'R':
	case 'r':
		n[4] = 1;
		if (n[3] != 4 && n[3] != 8)
			return false;
		break;
	case 'I':
	case 'i':
		n[4] = 0;
		if (n[3] != 1 && n[3] != 2 && n[3] != 4)
			return false;
		if (n[3] > 1) {
			s = terminal_input("big-endian (b) or little-endian (l): ");
			if (s[0] == 'b' || s[0] == 'B')
				n[5] = -n[5];
		}
		break;
	default:
		return false;
	}
	return div*n[3] == size;
}

long long filesize(const char *filename) {
	ifstream is (filename, ifstream::binary);
	is.seekg (0, is.end);
	return is.tellg();
}

void set_grid_spacing(double *r, int setr) {
	if (setr < 0) {
		string s = terminal_input("Grid spacings (dx dy dz): ");
		if (s.size() > 4) {
			stringstream ss(s);
			ss >> r[0] >> r[1] >> r[2];
			for (int i = 0; i != 3; i++)
				if (r[i] <= 0.0)
					r[i] = 1.0;
		}
	}
}


// output the grid info in console terminal:
void display_grid_info() {
	const double *d = G.get_d();
	const double *o = G.get_r0();
	const float *L = G.get_L();
	const unsigned int *N = G.get_N();
	unsigned int nx = N[0] + 1, ny = N[1] + 1, nz = N[2] + 1;
	cout << "\n ***** GRID INFO *****\nNum. data:\n\t" << nx << " x " << ny << " x " << nz;
	cout << "\nDimensions:\n\t" << L[0] << " x " << L[1] << " x " << L[2];
	cout << "\nOrigin:\n\t(" << o[0] << ", " << o[1] << ", " << o[2];
	cout << ")\nstep:\n\t" << d[0] << ", " << d[1] << ", " << d[2];
	cout << "\nmemory used:\n\t" << (sizeof(grid3d) + nz*(ny + 1)*sizeof(void*) +
	nx*ny*nz*sizeof(GRD_data_type))*9.5367431640625e-7 << " MiB\n\n";
}


int open_grid_file(const char *s) {
	int t = -1, setr = -1;
	double r[3] = {1.0, 1.0, 1.0};
	const char *c = strrchr(s, '.');
	if (c) {
		if (!strcmp(c,".grb"))
			setr = t = G.read_grd_binary(s);
		else if (!strcmp(c,".grd"))
			setr = t = G.read_grd(s);
		else if (!strcmp(c,".dat")) {
			set_grid_spacing(r, setr);
			setr = 1;
			t = G.read_dat_file(s);
		}
	}
	if (t) {
		long long size = filesize(s);
		switch ((int)size) {
		case 8192:
			set_grid_spacing(r, setr);
			t = G.read_scanfiles(s, 64, 0);
			break;
		case 131072:
			set_grid_spacing(r, setr);
			t = G.read_scanfiles(s, 256, 1);
			break;
		case 524288:
			set_grid_spacing(r, setr);
			t = G.read_scanfiles(s, 512, 1);
			break;
		default:
			{
				int n[6];
				if (input_res(n, size)) {
					set_grid_spacing(r, setr);
					t = G.read_raw_file(s, reinterpret_cast<unsigned int*>(n), n[5], n[4]);
				}
			}
		}
	}
	if (t) {
		cout << "\nError loading grid file!\n";
	} else {
		surf = 0;
		if (setr)
			G.set_ratio_aspect(r[0], r[1], r[2]);
	}
	return !t;
}

// set the scale and center of the graphic window
void set_scale() {
	const double *r0 = G.get_r0();
	const float *L = G.get_L();
	MC33_real r[3] = {0.5f*L[0], 0.5f*L[1], 0.5f*L[2]};
	if (G.isnotorthogonal())
		mult_Abf(G.get__A(), r, r, 0);
	cX = r[0] + r0[0];
	cY = r[1] + r0[1];
	cZ = r[2] + r0[2];
	int i = L[1] > L[0];
	if (L[2] > L[i])
		i = 2;
	scaleSize = (2.0f/3)*L[i];
}

/*********
Test functions used to generate the grid (implicit surfaces)
double f(double x, double y, double z)
**********/

double f1(double x, double y, double z) {
	return (sin(x*y) + sin(y*z) + sin(x*z))/(1.0 + x*x + y*y + z*z);
}

// https://blogs.ams.org/visualinsight/2016/04/15/barth-sextic/
#define phi ((sqrt(5) + 1)/2)
double Barthsextic(double x, double y, double z) {
	const double phi2 = phi*phi, w2 = 1.0*1.0;
	x *= x; y *= y; z *= z;
	double t = (x + y + z - w2)*w2;
	return 4*(phi2*x - y)*(phi2*y - z)*(phi2*z - x) - (1 + 2*phi)*t*t;
}

double ga(double cx, double cy, double cz, double x, double y, double z) {
	x -= cx; y -= cy; z -= cz;
	return exp(-0.5*(x*x + y*y + z*z));
}

double fourGaussian(double x, double y, double z) {
	return ga(-1.5, -1.5, -1.5, x, y, z) + ga(1.5, 1.5, -1.5, x, y, z) +
		ga(1.5, -1.5, 1.5, x, y, z) + ga(-1.5, 1.5, 1.5, x, y, z);
}

double Cube_f(double x, double y, double z) {
	x = fabs(x); y = fabs(y); z = fabs(z);
	if (x > y) {
		if (x > z)
			z = x;
	} else {
		if (y > z)
			z = y;
	}
	return 1.0 - (1.0/3.0)*z;
}

double Cylinder_f(double x, double y, double z) {
	double t = sqrt(x*x + z*z);
	y = sqrt(2.0)*fabs(y);
	if (t < y)
		t = y;
	return 1.0 - (1.0/3.0)*t;
}

//(9z^2 - 1)(1 - z^2) - 2y(y^2 - 3x^2)(1 - z^2) - (x^2+y^2)^2
// https://en.wikipedia.org/wiki/Implicit_surface
double Genus2(double x, double y, double z) {
	double t = y*y;
	x *= x; z = z*z - 1.0;
	y = (2.0*y*(t - 3.0*x) - 9.0*z - 8.0)*z;
	x += t;
	return y - x*x;
}

// -x^4 + 5*x^2 - y^4 + 5*y^2 - z^4 + 5*z^2 - 11.8
// https://www-sop.inria.fr/galaad/surface/
double TangleCube(double x, double y, double z) {
	x *= x; y *= y; z *= z;
	return x*(5 - x) + y*(5 - y) + z*(5 - z) - 11.8;
}

// https://en.wikipedia.org/wiki/Implicit_surface#Smooth_approximations_of_several_implicit_surfaces
double threeTori(double x, double y, double z) {
	const double R2 = 1.0*1.0, a2 = 0.2*0.2, r = 0.01;
	x *= x; y *= y; z *= z;
	double t = x + y + z + R2 - a2;
	t *= t;
	return r - (t - 4*R2*(x + y))*(t - 4*R2*(x + z))*(t - 4*R2*(y + z));
}

// http://virtualmathmuseum.org/Surface/decocube/decocube.html
const double cc = 1.3;
double DecoCube(double x, double y, double z) {
	const double c2 = 2.0 - cc*cc, ff = 0.02;
	x = x*x - 1.0; y = y*y - 1.0; z = z*z - 1.0;
	double t1 = x + y + c2, t2 = y + z + c2, t3 = z + x + c2;
	t1 *= t1; t2 *= t2; t3 *= t3;
	return ff - (t1 + z*z)*(t2 + x*x)*(t3 + y*y);
}

// https://math.stackexchange.com/q/46222
double f(double x, double y, double z) {
	x *= x; y *= y; z *= z;
	double t1 = (2.92*(x - 1)*x + 1.7*y)*(y - 0.88);
	double t2 = (2.92*(y - 1)*y + 1.7*z)*(z - 0.88);
	double t3 = (2.92*(z - 1)*z + 1.7*x)*(x - 0.88);
	return 0.04 - t1*t1 - t2*t2 - t3*t3;
}
//(2.92(x-1)x^2(x+1)+1.7y^2)^2*(y^2-0.88)^2 +
//(2.92(y-1)y^2(y+1)+1.7z^2)^2*(z^2-0.88)^2 +
//(2.92(z-1)z^2(z+1)+1.7x^2)^2*(x^2-0.88)^2 - 0.04
/******* End of test functions *********/

int main(int argc, char **argv) {

	int i = 1;
	if (argc == 1 || !open_grid_file(argv[1])) {
		cout << "To read a grid file, close this program and type:\n" << argv[0] << " filename\n";
		cout << "\nSelect one of the default grids:\n";
		cout << "\t1. (sin(x*y) + sin(y*z) + sin(x*z))/(1 + x*x + y*y + z*z)\n";
		cout << "\t2. Barth sextic\n";
		cout << "\t3. Four Gaussian functions, isovalue (0, 1)\n";
		cout << "\t4. Cube, isovalue [0, 1)\n";
		cout << "\t5. Cylinder, isovalue [-0.4142, 1)\n";
		cout << "\t6. (9z^2 - 1)(1 - z^2) - 2y(y^2 - 3x^2)(1 - z^2) - (x^2+y^2)^2\n";
		cout << "\t7. -x^4 + 5*x^2 - y^4 + 5*y^2 - z^4 + 5*z^2 - 11.8\n";
		cout << "\t8. Three intersecting tori, isovalue 0\n";
		cout << "\t9. DecoCube, isovalue 0\n";
		cout << "\t0. Leo cube?, isovalue 0\n";
		char c;
		do
			c = cin.get();
		while (c < '0' || c > '9');
		switch (c) {
			case '1':
				G.generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, f1);
				break;
			case '2':
				G.generate_grid_from_fn(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0, 0.02, 0.02, 0.02, Barthsextic);
				break;
			case '3':
				G.generate_grid_from_fn(-4.0, -4.0, -4.0, 4.0, 4.0, 4.0, 0.04, 0.04, 0.04, fourGaussian);
				break;
			case '4':
				G.generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, Cube_f);
				break;
			case '5':
				G.generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.04, 0.04, 0.04, Cylinder_f);
				break;
			case '6':
				G.generate_grid_from_fn(-1.875, -2.2, -1.875, 1.875, 1.55, 1.875, 0.025, 0.025, 0.025, Genus2);
				break;
			case '7':
				G.generate_grid_from_fn(-3.0, -3.0, -3.0, 3.0, 3.0, 3.0, 0.03, 0.03, 0.03, TangleCube);
				break;
			case '8':
				G.generate_grid_from_fn(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5, 0.015, 0.015, 0.015, threeTori);
				break;
			case '9':
				G.generate_grid_from_fn(-1.5, -1.5, -1.5, 1.5, 1.5, 1.5, 0.015, 0.015, 0.015, DecoCube);
				break;
			case '0':
				G.generate_grid_from_fn(-1.2, -1.2, -1.2, 1.2, 1.2, 1.2, 0.006, 0.006, 0.006, f);
				break;
			default:
				i = 0;
		}
	}
	if (i) {
		display_grid_info();
		glutInit(&i, argv);
		glutInitWindowSize(WIDTH, HEIGHT);
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowPosition(50, 50);
		glutCreateWindow("Marching cubes 33");

		init_param();
		MC.set_grid3d(G);
		set_scale();
		reset_view();

		glutDisplayFunc(draw); // set the render function
		glutKeyboardFunc(keyboard); // register the keyboard functions
		glutReshapeFunc(reshape); // register the reshape function
		hide_message();
		glutMainLoop();
	}
	return 0;
}

