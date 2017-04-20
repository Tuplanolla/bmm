flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -O0 -g \
	-Weverything \
	-Wno-aggregate-return -Wno-bad-function-cast -Wno-disabled-macro-expansion \
	-Wno-switch
endif
ifeq ($(CONFIG), profile)
flags=-DNDEBUG -O3 -g -save-temps
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -O3 -Wl,-s -w
endif
# TODO Remove this.
flags+=-Wno-used-but-marked-unused
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -Og -g \
	`cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat \
	-Wno-long-long -Wno-traditional -Wno-traditional-conversion \
	-Wno-declaration-after-statement -Wno-unsuffixed-float-constants \
	-Wno-address -Wno-aggregate-return \
	-Wno-switch -Wno-switch-enum -Wno-switch-default \
	-Wno-missing-declarations -Wno-missing-prototypes
endif
ifeq ($(CONFIG), profile)
flags=-D_GNU_SOURCE -DNDEBUG -O3 -g -pg -save-temps
endif
ifeq ($(CONFIG), release)
flags=-D_GNU_SOURCE -DNDEBUG -O3 -s -w
endif
endif

# TODO These are a bit messy.

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 $(flags)
LDLIBS=-lm -lrt
CFLAGSGSL=`pkg-config --cflags cheat gsl`
LDLIBSGSL=`pkg-config --libs cheat gsl`
CFLAGSSDL=`pkg-config --cflags freeglut gl glew gsl sdl2`
LDLIBSSDL=`pkg-config --libs freeglut gl glew gsl sdl2`

build: bmm-dem bmm-filter bmm-glut bmm-sdl

run: build
	GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=42 time -v \
	./bmm-dem | \
	./bmm-filter --mode whitelist --pass 142 --pass 144 --pass 168 --verbose yes | \
	./bmm-sdl

run-to: build
	./bmm-dem | \
	./bmm-filter --mode whitelist --pass 142 --pass 144 --pass 168 | \
	gzip -c > bmm.run.gz

run-from: build
	gunzip -c < bmm.run.gz | \
	./bmm-sdl

start-server: bmm-sdl
	mkfifo bmm.fifo
	./bmm-sdl < bmm.fifo &

run-client: bmm-dem bmm.fifo
	./bmm-dem > bmm.fifo

stop-server: bmm.fifo
	$(RM) bmm.fifo

check-static: build
	cppcheck -I/usr/include --enable=all *.c *.h

check-dynamic: build
	valgrind --tool=memcheck --leak-check=full ./bmm-dem > /dev/null

profile-valgrind: build
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
	./bmm-dem > /dev/null
	callgrind_annotate callgrind.out

profile-gprof: build
	./bmm-dem > /dev/null
	gprof ./bmm-dem

test: tests
	./tests

deep-clean: clean
	$(RM) *.data *.log *.out *.run

clean: shallow-clean
	$(RM) bmm-dem bmm-filter bmm-glut bmm-sdl tests

shallow-clean:
	$(RM) *.gch *.i *.o *.s

bmm-dem: bmm-dem.o \
	dem.o endy.o fp.o geom.o geom2d.o hack.o io.o msg.o \
	opt.o sec.o sig.o size.o str.o tle.o
	$(CC) $(CFLAGS) $(CFLAGSGSL) -o $@ $^ $(LDLIBS) $(LDLIBSGSL)

bmm-filter: bmm-filter.o \
	endy.o filter.o fp.o hack.o io.o msg.o \
	opt.o sec.o sig.o size.o str.o tle.o
	$(CC) $(CFLAGS) $(CFLAGSGSL) -o $@ $^ $(LDLIBS) $(LDLIBSGSL)

bmm-glut: bmm-glut.o \
	dem.o endy.o fp.o geom.o geom2d.o gl.o hack.o io.o msg.o \
	opt.o sec.o sig.o size.o str.o tle.o
	$(CC) $(CFLAGS) $(CFLAGSSDL) -o $@ $^ $(LDLIBS) $(LDLIBSSDL)

bmm-sdl: bmm-sdl.o \
	dem.o endy.o fp.o geom.o geom2d.o gl.o hack.o io.o msg.o \
	opt.o sdl.o sec.o sig.o size.o str.o tle.o
	$(CC) $(CFLAGS) $(CFLAGSSDL) -o $@ $^ $(LDLIBS) $(LDLIBSSDL)

tests: tests.o \
	dem.o endy.o fp.o geom.o geom2d.o hack.o io.o msg.o \
	opt.o sec.o sig.o size.o str.o tle.o
	$(CC) $(CFLAGS) $(CFLAGSGSL) -o $@ $^ $(LDLIBS) $(LDLIBSGSL)

%.o: %.c *.h
	$(CC) $(CFLAGS) $(CFLAGSGSL) -c -o $@ $<
