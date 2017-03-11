flags=

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -O0 -g \
	-Weverything \
	-Wno-bad-function-cast -Wno-disabled-macro-expansion \
	-Wno-aggregate-return -Wno-covered-switch-default -Wno-unused-function
endif
ifeq ($(CONFIG), release)
flags=-DNDEBUG -O3 -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
flags=-D_GNU_SOURCE -DDEBUG -Og -g \
	`cat gcc-$$(./gcc-version | tr . _)-release` \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-declaration-after-statement \
	-Wno-traditional -Wno-traditional-conversion \
	-Wno-unsuffixed-float-constants \
	-Wno-address -Wno-bad-function-cast -Wno-long-long \
	-Wno-aggregate-return -Wno-switch-default -Wno-unused-function
endif
ifeq ($(CONFIG), release)
flags=-D_GNU_SOURCE -DNDEBUG -O3 -s -w
endif
endif

CFLAGS=-D_POSIX_C_SOURCE=200809L -std=c11 $(flags)
LDLIBS=-lm -lrt
CFLAGSGSL=`pkg-config --cflags gsl`
LDLIBSGSL=`pkg-config --libs gsl`
CFLAGSSDL=`pkg-config --cflags sdl`
LDLIBSSDL=`pkg-config --libs sdl`

build: bmm-dem bmm-sdl

run: build
	GSL_RNG_TYPE=mt19937 GSL_RNG_SEED=0 time -v ./bmm-dem

check: build
	cppcheck -I/usr/include --enable=all *.c *.h
	valgrind --leak-check=full --tool=memcheck ./bmm-dem

deep-clean: clean
	$(RM) run-*

clean: shallow-clean
	$(RM) bmm-dem bmm-sdl

shallow-clean:
	$(RM) *.gch *.o

bmm-dem: bmm-dem.o clopts.o errors.o strs.o
	$(CC) $(CFLAGS) $(CFLAGSGSL) -o $@ $^ $(LDLIBS) $(LDLIBSGSL)

bmm-sdl: bmm-sdl.o
	$(CC) $(CFLAGS) $(CFLAGSSDL) -o $@ $^ $(LDLIBS) $(LDLIBSSDL)
