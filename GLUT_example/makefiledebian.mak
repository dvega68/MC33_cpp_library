#--------------------------------------------------------
CPP      = g++
SOURCE   = TestMC33_glut.cpp
OBJ      = $(SOURCE:.cpp=.o)
LIBS     = -lglut -lGLU -lGL -ldl -lm -s
CPPINCS  =
BIN      = TestMC33_glut
OPTIM		 = -Ofast -Wall -Wextra -funroll-loops
CPPFLAGS = $(CPPINCS) -std=c++11 $(OPTIM)
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

all:	all-before $(BIN) all-after

$(BIN): $(OBJ)
	#$(SOURCE)
	$(CPP) -o $(BIN) $(OBJ) $(LIBS)

clean: clean-custom
	$(RM) $(OBJ) $(BIN)
#--------------------------------------------------------
