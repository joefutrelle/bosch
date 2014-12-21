# depends on libtiff-dev, libsndfile-dev, libevent-pthreads-2.0-5
CC = gcc
CFLAGS = -O3
# -g -DDEBUG
OBJ = bosch.o
EXE = bosch

#LIB_DIR = -L/opt/local/lib
LIBS = -lm -ltiff -lsndfile -lpthread

INC_DIR = -I/opt/local/include

CFLAGS = $(INC_DIR)

all: $(OBJ) $(EXE)

clean:
	rm -f $(OBJ) $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $(EXE) $(INC_DIR) $(LIB_DIR) $(LIBS)
