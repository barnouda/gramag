#CC=icc

#FLAGS=-g
#FLAGS=-DDEBUG
#FLAGS=-DDEBUGMAG
#FLAGS=-O3
#FLAGS=-openmp -DNTHREAD=16
#FLAGS=-parallel
#FLAGS=-fopenmp

PROGS=gralite grafull maglite magfull

BIN=./bin
BUILD=./build
SRC=./src

CC=gcc
FLAGS=-O3 -I./$(SRC)/cubature-1.0
ENDFLAGS=-lm

# makefile tuto
#target: prerequistes
#	command
# $@ = target
# $< = first prerequiste
# $^ = all prerequistes

$(BUILD)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c -o $@ $^

gralite: $(BUILD)/maingralite.o $(BUILD)/interp.o $(BUILD)/loop8.o $(BUILD)/gravilib.o $(BUILD)/primitive8.o $(BUILD)/rhoVp.o $(BUILD)/topprism8.o
	$(CC) $(FLAGS) -o $(BIN)/gralite $^ $(ENDFLAGS)

grafull: $(BUILD)/maingrafull.o $(BUILD)/interp.o $(BUILD)/loop9.o $(BUILD)/gravilib.o $(BUILD)/primitive9.o $(BUILD)/rhoVp.o $(BUILD)/outderivative.o $(BUILD)/topprism9.o 
	$(CC) $(FLAGS) -o $(BIN)/grafull $^ $(ENDFLAGS)

maglite: $(BUILD)/mainmaglite.o $(BUILD)/interp.o $(BUILD)/loop10.o $(BUILD)/gravilib.o $(BUILD)/cubature-1.0/hcubature.o $(BUILD)/maglib.o $(BUILD)/primitive10.o $(BUILD)/topprism10.o
	$(CC) $(FLAGS) -o $(BIN)/maglite $^ $(ENDFLAGS)

magfull: $(BUILD)/mainmagfull.o $(BUILD)/interp.o $(BUILD)/loop11.o $(BUILD)/gravilib.o $(BUILD)/cubature-1.0/hcubature.o $(BUILD)/maglib.o $(BUILD)/primitive11.o $(BUILD)/outderivative.o $(BUILD)/topprism11.o
	$(CC) $(FLAGS) -o $(BIN)/magfull $^ $(ENDFLAGS)
 
all: $(PROGS)

clean:
	/bin/rm -f $(BUILD)/*.o $(BUILD)/cubature-1.0/hcubature.o $(BIN)/*

