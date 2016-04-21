CC = g++
CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(notdir $(CPP_FILES:.cpp=.o))
LD_FLAGS = 
CC_FLAGS =-Wall -O3 -std=c++11

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
