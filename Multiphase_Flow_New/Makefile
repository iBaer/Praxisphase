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


constants.o: constants.h
Raster.o: Raster.h constants.h
Gleichungssystem.o: Gleichungssystem.h constants.h Zelle.h
LaxFriedrichMethod.o: LaxFriedrichMethod.h numerische_methode.h Gleichungssystem.h constants.h numerische_methode.h
FORCE.o: FORCE.h Gleichungssystem.h constants.h numerische_methode.h
Zelle.o: Zelle.h
main.o: LaxFriedrichMethod.h FORCE.h numerische_methode.h
numerische_methode.o: numerische_methode.h constants.h Raster.h constants.h Gleichungssystem.h







