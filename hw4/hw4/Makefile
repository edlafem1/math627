OBJS := main.o memory.o utilities.o lin_alg.o

EXECUTABLE := power

# choose compiler:
CC := mpiicc
# CC := gcc

# choose flags:
# flags for Intel compiler icc on maya:
CFLAGS := -O3 -std=c99 -Wall -mkl
# flags for Portland Group compiler pgcc on maya:
# CFLAGS := -O3 -c99 -Minform=warn -fastsse
# flags for GNU compiler gcc anywhere:
# CFLAGS := -O3 -std=c99 -Wall -Wno-unused-variable

DEFS := -DBLAS -DPARALLEL
#  -DPARALLEL
# -DBLAS
INCLUDES :=
LDFLAGS := -lm

%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(DEFS) $(INCLUDES) $(OBJS) -o $@ $(LDFLAGS)

clean:
	-rm -f *.o $(EXECUTABLE)

