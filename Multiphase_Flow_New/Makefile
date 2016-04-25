CC = g++
CPP_FILES = $(wildcard *.cpp)
OBJ_FILES = $(notdir $(CPP_FILES:.cpp=.o))
LD_FLAGS = -L/usr/local/lib/
CC_FLAGS =-Wall -O3 -std=c++11 -I/usr/local/include/

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


FORCE.o: FORCE.h Gleichungssystem.h Konstanten.h
Gleichungssystem.o: Gleichungssystem.h Konstanten.h Zelle.h
Konstanten.o: Konstanten.h
LaxFriedrichMethod.o: LaxFriedrichMethod.h numerische_methode.h Gleichungssystem.h Konstanten.h
Raster.o: Raster.h Konstanten.h
Zelle.o: Zelle.h
main.o: LaxFriedrichMethod.h FORCE.h numerische_methode.h
numerische_methode.o: numerische_methode.h Konstanten.h Raster.h Konstanten.h Gleichungssystem.h







