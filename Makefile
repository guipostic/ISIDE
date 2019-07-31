CC      = g++
CFLAGS  = -O3 -fopenmp -march=native -std=c++11 -c -Wall -Wextra
LDFLAGS = -fopenmp
SRC     = $(wildcard src/*.cpp) $(wildcard src/PTools/*.cpp)
OBJ     = $(SRC:.cpp=.o)
EXE     = ISIDE

all: $(SRC) $(EXE)
	
$(EXE): $(OBJ) 
	$(CC) $(LDFLAGS) $(OBJ) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJ) $(EXE)
