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
	# Run with the arguments in this order.
	# ./bmm-dem --stat yes --ds 0.5 --pfac 1.0 --script fin --trap no --verbose yes
	./bmm-dem --script shear --trap no --verbose yes | \
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

# The rest is automatically generated by `gcc -MM *.c`.

bmm-dem.o: bmm-dem.c dem.h conf.h cpp.h geom2d.h common.h alias.h ext.h \
 common_mono.h common_poly.h common_ord.h common_num.h common_int.h \
 common_sint.h common_uint.h common_fp.h wrap.h fp.h ival.h io.h msg.h \
 endy.h msg_.h geom.h opt.h str.h tle.h tle_.h
bmm-filter.o: bmm-filter.c conf.h cpp.h ext.h filter.h io.h msg.h endy.h \
 msg_.h opt.h str.h tle.h tle_.h
bmm-glut.o: bmm-glut.c glut.h ext.h cpp.h io.h opt.h str.h tle.h tle_.h
bmm-nc.o: bmm-nc.c ext.h cpp.h nc.h io.h opt.h str.h tle.h tle_.h
bmm-sdl.o: bmm-sdl.c sdl.h dem.h conf.h cpp.h geom2d.h common.h alias.h \
 ext.h common_mono.h common_poly.h common_ord.h common_num.h common_int.h \
 common_sint.h common_uint.h common_fp.h wrap.h fp.h ival.h io.h msg.h \
 endy.h msg_.h opt.h str.h tle.h tle_.h
common.o: common.c alias.h common.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h common_mono.c common_poly.c \
 common_ord.c common_num.c common_int.c common_sint.c common_uint.c \
 common_fp.c
common_fp.o: common_fp.c ext.h cpp.h
common_int.o: common_int.c ext.h cpp.h
common_mono.o: common_mono.c
common_num.o: common_num.c ext.h cpp.h
common_ord.o: common_ord.c ext.h cpp.h
common_poly.o: common_poly.c ext.h cpp.h
common_sint.o: common_sint.c ext.h cpp.h
common_uint.o: common_uint.c ext.h cpp.h
concat.o: concat.c
dem.o: dem.c common.h alias.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h conf.h dem.h geom2d.h fp.h ival.h io.h msg.h endy.h \
 msg_.h geom.h kernel.h neigh.h random.h sig.h tle.h tle_.h
endy.o: endy.c endy.h ext.h cpp.h
fact.o: fact.c
filter.o: filter.c conf.h cpp.h filter.h ext.h io.h fp.h common.h alias.h \
 common_mono.h common_poly.h common_ord.h common_num.h common_int.h \
 common_sint.h common_uint.h common_fp.h wrap.h msg.h endy.h msg_.h sig.h \
 tle.h tle_.h
fp.o: fp.c fp.h common.h alias.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h
geom.o: geom.c geom.h common.h alias.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h fp.h
geom2d.o: geom2d.c geom2d.h common.h alias.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h fp.h ival.h
gl.o: gl.c fp.h common.h alias.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h gl.h
gl2.o: gl2.c gl2.h ext.h cpp.h io.h tle.h tle_.h
glut.o: glut.c gl2.h ext.h cpp.h glut.h io.h conf.h dem.h geom2d.h \
 common.h alias.h common_mono.h common_poly.h common_ord.h common_num.h \
 common_int.h common_sint.h common_uint.h common_fp.h wrap.h fp.h ival.h \
 msg.h endy.h msg_.h sig.h tle.h tle_.h
hack.o: hack.c hack.h
io.o: io.c io.h ext.h cpp.h tle.h tle_.h
ival.o: ival.c ival.h common.h alias.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h fp.h
kernel.o: kernel.c kernel.h common.h alias.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h fp.h
maskbits.o: maskbits.c
msg.o: msg.c common.h alias.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h conf.h endy.h io.h msg.h msg_.h tle.h tle_.h
nc-ex.o: nc-ex.c
nc.o: nc.c conf.h cpp.h dem.h geom2d.h common.h alias.h ext.h \
 common_mono.h common_poly.h common_ord.h common_num.h common_int.h \
 common_sint.h common_uint.h common_fp.h wrap.h fp.h ival.h io.h msg.h \
 endy.h msg_.h nc.h sig.h tle.h tle_.h
neigh.o: neigh.c neigh.h common.h alias.h ext.h cpp.h common_mono.h \
 common_poly.h common_ord.h common_num.h common_int.h common_sint.h \
 common_uint.h common_fp.h wrap.h
opt.o: opt.c opt.h ext.h cpp.h tle.h tle_.h
pow.o: pow.c
random.o: random.c random.h ext.h cpp.h
sdl.o: sdl.c common.h alias.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h dem.h conf.h geom2d.h fp.h ival.h io.h msg.h endy.h \
 msg_.h gl.h sdl.h tle.h tle_.h
sec.o: sec.c sec.h ext.h cpp.h
sig.o: sig.c sig.h
splice.o: splice.c
str.o: str.c str.h ext.h cpp.h tle.h tle_.h
tests.o: tests.c alias.h common.h ext.h cpp.h common_mono.h common_poly.h \
 common_ord.h common_num.h common_int.h common_sint.h common_uint.h \
 common_fp.h wrap.h endy.h fp.h geom2d.h ival.h neigh.h msg.h io.h msg_.h
tle.o: tle.c ext.h cpp.h hack.h sec.h tle.h tle_.h
wrap.o: wrap.c ext.h cpp.h wrap.h alias.h
