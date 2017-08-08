CFLAGS+=-D_POSIX_C_SOURCE=200809L -std=c11
LDLIBS+=-lm -lrt

ifeq ($(CC), clang)
ifeq ($(CONFIG), debug)
CFLAGS+=-D_GNU_SOURCE -DDEBUG -O0 -ftrapv -g \
	-Weverything \
	-Wno-aggregate-return -Wno-bad-function-cast -Wno-disabled-macro-expansion \
	-Wno-missing-prototypes -Wno-padded -Wno-unused-parameter -Wno-attributes -Wno-unused-function -Wno-shadow \
	-Wno-switch -Wno-used-but-marked-unused \
	-Wno-dollar-in-identifier-extension
# TODO The third to last line should be removed later.
endif
ifeq ($(CONFIG), profile)
CFLAGS+=-D_GNU_SOURCE -DNDEBUG -O3 -g -save-temps
endif
ifeq ($(CONFIG), release)
CFLAGS+=-D_GNU_SOURCE -DNDEBUG -O3 -Wl,-s -w
endif
endif

ifeq ($(CC), gcc)
ifeq ($(CONFIG), debug)
CFLAGS+=-D_GNU_SOURCE -DDEBUG -Og -ftrapv -g \
	$$(cat gcc-$$(./gcc-version | tr . _)-release) \
	-Wno-error -Wno-fatal-errors -Wno-system-headers \
	-Wno-c++-compat -Wno-c90-c99-compat \
	-Wno-long-long -Wno-traditional -Wno-traditional-conversion \
	-Wno-declaration-after-statement -Wno-unsuffixed-float-constants \
	-Wno-address -Wno-aggregate-return \
	-Wno-switch -Wno-switch-enum -Wno-switch-default \
	-Wno-missing-prototypes -Wno-padded -Wno-unused-parameter -Wno-attributes -Wno-unused-function -Wno-shadow \
	-Wno-missing-declarations -Wno-missing-prototypes
# TODO The second to last line should be removed later.
endif
ifeq ($(CONFIG), profile)
CFLAGS+=-D_GNU_SOURCE -DNDEBUG -O3 -g -pg -save-temps
endif
ifeq ($(CONFIG), release)
CFLAGS+=-D_GNU_SOURCE -DNDEBUG -O3 -s -w
endif
endif

build: bmm-dem bmm-filter bmm-glut bmm-nc bmm-sdl

run: bmm-dem bmm-filter bmm-sdl
	./bmm-dem --script shear --trap no --verbose yes | \
	./bmm-filter --mode whitelist --pass opts --pass istep --pass parts --pass neigh --verbose yes | \
	./bmm-sdl

run-store: bmm-dem bmm-filter
	./bmm-dem --script mix | \
	./bmm-filter --mode whitelist --pass npart --pass parts --pass neigh | \
	gzip -c > bmm.run.gz

run-load: bmm-sdl
	gunzip -c < bmm.run.gz | \
	./bmm-sdl

run-client: bmm-dem bmm.fifo
	./bmm-dem > bmm.fifo

start-server: bmm-sdl
	mkfifo bmm.fifo
	./bmm-sdl < bmm.fifo &

stop-server: bmm.fifo
	$(RM) bmm.fifo

check-static: build
	cppcheck -I/usr/include --enable=all *.c *.h | grep '^\['

check-dynamic: build
	valgrind --tool=memcheck --leak-check=full ./bmm-dem --script couple > /dev/null

profile-sample: build
	./bmm-dem > /dev/null
	gprof ./bmm-dem

profile-emulate: build
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out \
	./bmm-dem > /dev/null
	callgrind_annotate callgrind.out

test: tests
	./tests

deep-clean: clean
	$(RM) *.data *.log *.out *.run

clean: shallow-clean
	$(RM) bmm-dem bmm-filter bmm-glut bmm-nc bmm-sdl tests

shallow-clean:
	$(RM) *.gch *.i *.o *.s

bmm-dem: CFLAGS+=$$(pkg-config --cflags gsl)
bmm-dem: LDLIBS+=$$(pkg-config --libs gsl)
bmm-dem: bmm-dem.o \
	common.o dem.o endy.o fp.o geom.o geom2d.o hack.o kernel.o io.o ival.o msg.o \
	neigh.o opt.o random.o sec.o sig.o str.o tle.o wrap.o

bmm-filter: bmm-filter.o \
	common.o endy.o filter.o fp.o hack.o kernel.o io.o msg.o \
	opt.o sec.o sig.o str.o tle.o wrap.o

bmm-glut: CFLAGS+=$$(pkg-config --cflags freeglut gl glew gsl)
bmm-glut: LDLIBS+=$$(pkg-config --libs freeglut gl glew gsl)
bmm-glut: bmm-glut.o \
	common.o dem.o endy.o fp.o geom.o geom2d.o gl2.o glut.o hack.o kernel.o io.o ival.o msg.o \
	neigh.o opt.o random.o sec.o sig.o str.o tle.o wrap.o

bmm-nc: CFLAGS+=$$(pkg-config --cflags netcdf)
bmm-nc: LDLIBS+=$$(pkg-config --libs netcdf)
bmm-nc: bmm-nc.o \
	common.o endy.o fp.o hack.o kernel.o io.o msg.o \
	nc.o opt.o sec.o sig.o str.o tle.o wrap.o

bmm-sdl: CFLAGS+=$$(pkg-config --cflags freeglut gl gsl sdl2)
bmm-sdl: LDLIBS+=$$(pkg-config --libs freeglut gl gsl sdl2)
bmm-sdl: bmm-sdl.o \
	common.o dem.o endy.o fp.o geom.o geom2d.o gl.o hack.o kernel.o io.o ival.o msg.o \
	neigh.o opt.o sdl.o random.o sec.o sig.o str.o tle.o wrap.o

tests: CFLAGS+=$$(pkg-config --cflags cheat gsl)
tests: LDLIBS+=$$(pkg-config --libs cheat gsl)
tests: tests.o \
	common.o endy.o fp.o geom.o geom2d.o hack.o kernel.o io.o ival.o msg.o \
	neigh.o opt.o random.o sec.o sig.o str.o tle.o wrap.o
