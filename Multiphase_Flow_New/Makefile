CC = g++
CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(notdir $(CPP_FILES:.cpp=.o))
LD_FLAGS = -L/usr/local/lib/
CC_FLAGS =-Wall -O3 -std=c++11 -I/usr/local/include/ -static -g

all: main


debug: CC_FLAGS+= -DDEBUG 
debug: $(OBJ_FILES)
	$(CC) $(LD_FLAGS) -o $@ $^

main: $(OBJ_FILES)
	$(CC) $(LD_FLAGS) -o $@ $^ -llapackpp

%.o: %.cpp
	$(CC) $(CC_FLAGS) -c -o $@ $<

$(TARGET): $(OBJ_FILES)
	$(CC) $(OBJ_FILES) -o $(TARGET) $(LDFLAGS)

clean: 
	rm -f $(OBJ_FILES) $(TARGET) core main debug


constants.o: constants.h
grid.o: grid.h constants.h
computation.o: computation.h constants.h
solver.o: solver.h computation.h constants.h
lax_friedrich.o: lax_friedrich.h constants.h solver.h
force.o: force.h computation.h constants.h solver.h
main.o: lax_friedrich.h force.h numerical_method.h
numerical_method.o: numerical_method.h constants.h grid.h constants.h computation.h solver.h







