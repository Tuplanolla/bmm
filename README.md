# Brittle Matter Matters

This repository contains a simplified implementation
of discrete element method as applied to brittle surface contact.
While the code seems to produce correct results and run reasonably fast,
it is primarily intended for pedagogical purposes.
Clarity and conciseness of the presentation are favored
over numerical stability and performance tricks.
If you would prefer something with the opposite design goals, see

* [University of Oslo][uio]'s [ESyS-Particle][esys-particle] under OSL-3.0,
* [Yade Project][yade-dem]'s [Yade][yade] under GPL-2.0 or
* [CFDEM][cfdem]'s [LIGGGHTS][liggghts] under GPL-2.0.

## Overview

Brittle Matter Matters was written by Sampsa "Tuplanolla" Kiiskinen
to support his master's thesis on complex material physics.
The project was written between 2017-03-09 and 2017-08-31.

## License

Brittle Matter Matters is free software and as such
licensed under the GNU General Public License version 3 or later.
The full license can be found in the `LICENSE` file that
resides in the same directory as this file.
In short, copies and derivative works are permitted
as long as they use a compatible license.

## The Plan

### The Goal

Find $\\mu_{macro}(p, \\mu_{micro})$.

### Project Structure

The project is divided into several programs with different purposes.
This makes it easier to distribute over several machines
with different specializations.
One can, for example, run the simulation on a headless server
while viewing the visualization on a single-user workstation.

The main program `bmm-dem` does all the heavy lifting.
In addition to running the discrete element method and
producing two output streams of results,
it parses the command line options, reads the initialization files and
handles signals in a synchronous fashion.
Aside from these effects the program resembles a pure function.

The utility program `bmm-filter` removes certain messages from a stream.
If, for example, `bmm-dem` produces real-time progress reports,
they can be stripped away before the results are saved into a file.

The analysis programs `bmm-sdl` and `bmm-gp`
draw visualizations in different ways.
The former is a real-time visualizer built on top of SDL and
the latter is a batch visualizer
that produces data files and Gnuplot scripts for them.

New programs may pop up unexpectedly.

### Streams

The programs can be roughly classified
into producers, consumers and transformers.
Producers only write to their output streams and
consumers only read from their input streams
while transformers do both.
As an example, this classification makes `bmm-dem` a producer,
`bmm-sdl` a consumer and `bmm-filter` a transformer.

Since transformers participate in all parts of inter-process communication,
the overall streaming architecture is best illustrated from their perspective.
They read from the input stream `stdin` and
write to both the primary output stream `stdout` and
the secondary output stream `stderr`.

The primary output stream `stdout` is a high-bandwidth channel
for transporting a sequence of events that follow a strict messaging protocol.
It carries all the essential information of the simulation and
is thus intended to be piped into another process for visualization or
saved into a file for later analysis.
The input stream `stdin` is essentially dual to `stdout` and
as such another strictly regulated high-bandwidth channel.

The secondary output stream `stderr` is quite different
from `stdout` and `stdin` as it is a low-bandwidth channel
for unstructured diagnostics and error messages.
It is meant to be read by the user as is and
does not participate in inter-process communication.

If multiple threads are used,
some care must be taken to ensure the coherence of the streams.
For `stdout` and `stdin` this can be accomplished
with a simple lock or a size-restricted atomic use policy.
For `stderr` coherence is less imporant,
so line-buffering should be sufficient.

### Messaging Protocol

To achieve maximal throughput,
information is passed around in a compact binary format,
which is specifically designed for this use case.
Conceptually the format consists of streams,
which are finite sequences of messages.
Within each stream the messages are completely independent of each other,
so sequencing them works by mere concatenation.
Messages themselves come with a header and a body,
where the header has a fixed format,
but the structure of the body is determined by the header.
Just like streams,
messages are also formed from their constituents by concatenation.

At the lowest level of abstraction messages are defined
in terms of bits and bytes (octets) instead of other abstract concepts.
The following table contains all the possible bit patterns and their meanings.
In the table digits and repeated lowercase letters denote bits and
repeated uppercase letters denote bytes.
Free variables that carry no information should be zero by default.
Anything not mentioned in the table is assumed to be invalid and
sending or attempting to interpret such a thing is a protocol violation.

| Bit Pattern | Meaning
|:------------|:--------
| `0bbbdddd` | The message has low priority.
| `1bbbdddd` | The message has high priority.
| `a000dddd` | Integers are in little-endian (stream order).
| `a001dddd` | Integers are in middle-endian with 2-byte blocks swapped.
| `a010dddd` | Integers are in middle-endian with 4-byte blocks are swapped.
| `a011dddd` | Integers are in middle-endian with 2-byte and 4-byte blocks swapped.
| `a100dddd` | Integers are in middle-endian with 8-byte blocks swapped.
| `a101dddd` | Integers are in middle-endian with 2-byte and 8-byte blocks swapped.
| `a110dddd` | Integers are in middle-endian with 4-byte and 8-byte blocks swapped.
| `a111dddd` | Integers are in big-endian (network order).
| `abbb0ddd` | Message body has a fixed size of `d` bytes.
| `abbb0000` | Message body is empty (0 B).
| `abbb0001` | Message body is very small (1 B).
| `abbb0010` | Message body is small (2 B).
| `abbb0011` | Message body is medium small (3 B).
| `abbb0100` | Message body is medium (4 B).
| `abbb0101` | Message body is medium large (5 B).
| `abbb0110` | Message body is large (6 B).
| `abbb0111` | Message body is very large (7 B).
| `abbb1ddd` | Message body has a variable size.
| `abbb10dd` | Message body size fits into `d` bytes.
| `abbb1000 E` | Message body has a size of `e` bytes (256 B).
| `abbb1001 E E` | Message body has a size of `e` bytes (64 kiB).
| `abbb1010 E E E E` | Message body has a size of `e` bytes (4 GiB).
| `abbb1011 E E E E E E E E` | Message body has a size of `e` bytes (16 EiB).
| `abbb11dd` | Message body is terminated by a literal that spans `d` bytes.
| `abbb1100 E` | Message body is terminated by a literal `e` (1 B).
| `abbb1101 E E` | Message body is terminated by a literal `e` (2 B).
| `abbb1110 E E E E` | Message body is terminated by a literal `e` (4 B).
| `abbb1111 E E E E E E E E` | Message body is terminated by a literal `e` (8 B).

The purposes of some of the patterns overlap intentionally,
so that one can neglect to implement the more complex parts
(higher in bits to indicate) if the simpler ones are sufficient.

The first octet is called the flags and then comes the (possibly empty) prefix.
Together they form the header.
Afterwards comes the body, which begins with a (possibly empty) number.
The size in the prefix covers the whole body.

    | Header         | Body
    | Flags | Prefix | Number | Data

To summarize the table informally,
each message is prefixed by two nibbles,
the first of which contains user-set flags,
the second of which contains derived flags.
Additionally, for those messages whose bodies may vary in size,
there is an extra prefix to signal the size of the message.
The next byte contains the message number and
after that comes the message body.

In GNU C it would not go as follows.

    struct bmm_msg {
      uint8_t userset;
      uint8_t derived;
      union {
        union {};
        uint8_t octets0;
        uint16_t octets1;
        uint32_t octets2;
        uint64_t octets3;
      } size;
      enum bmm_msg_num num;
      unsigned char body[];
    };

For words of length $n$ bytes there are $n!$ byte ordering schemes.
To enumerate all byte orderings for words of $m = 8$ bytes,
as is the maximal message length,
one would need $b = \\lceil\\log_2 m!\\rceil = 16$ bits.

It is better to assume that byte ordering schemes are built from pair swaps,
so the permutations form a binary tree.
Now $d = \\lceil\\log_2 n\\rceil$ is the tree depth, starting from zero.
There are $2^{d + 1} - 1$ orderings and
enumerating them only takes $b = 4$ bits.

If all the swaps per level are constant,
that would result in $n$ choices and an enumeration of $b = 3$ bits.

One more cut could be made by requiring all the swaps are constant,
giving $2$ options and an enumeration of $b = 1$ bits.
These are known as little and big endian.

#### Index Files

To make browsing past simulations easier and more efficient,
one could generate index files from the recorded streams.
These would consist of fixed chunks that contain
the message number and flags as usual,
but instead of the message body there would simply be an offset number
to the message within the record.

Look at this familiar table.

| Bit Pattern | Meaning
|:------------|:--------
| `abbb0ddd` | Message offset has the endianness `b` and the range `d`.
| `abbb0000` | Message is at offset `g` (0 B).
| `abbb0001 G` | Message is at offset `g` (256 B).
| `abbb0010 G G` | Message is at offset `g` (64 kiB).
| `abbb0011 G G G` | Message is at offset `g` (16 MiB).
| `abbb0100 G G G G` | Message is at offset `g` (4 GiB).
| `abbb0101 G G G G G` | Message is at offset `g` (1 TiB).
| `abbb0110 G G G G G G` | Message is at offset `g` (256 TiB).
| `abbb0111 G G G G G G G` | Message is at offset `g` (64 PiB).
| `abbbdddd` | Protocol violation.

The range `d` has to stay the same throughout the index.

### Streaming

Compressing movies is easy.

    $ ./bmm-dem | gzip -c > bmm.run.gz
    $ gunzip -c < bmm.run.gz | ./bmm-sdl

Sending data through the network works fine.

    $ nc -l 9001 | ./bmm-sdl
    $ ./bmm-dem | nc 127.0.0.1 9001

### Analysis Software

ParaView specializes in

* Structural Analysis
* Fluid Dynamics
* Astrophysics
* Climate Science
* LiDAR/Point-Cloud

so it is not an option.
Also

* AtomEye
* PyMol
* Raster3d
* RasMol

are chemist-oriented.
We shall review

* VMD
* OVITO

next.
Also take a look at GLE.

### File Format Choices

OVITO is a bit picky about file formats;
table as of version 2.8.2.

| Path | Format | Binary | Import | Export
|:-----|:-------|:-------|:-------|:-------
| `ovito/src/plugins/particles/lammps` | LAMMPS Dump | Yes | Yes | Yes
| `ovito/src/plugins/particles/parcas` | PARCAS | Yes | Yes | No
| `ovito/src/plugins/netcdf` | NetCDF | Yes | Yes | No
| `ovito/src/plugins/particles/gsd` | GSD/HOOMD | Yes | Yes | No
| `ovito/src/plugins/particles/lammps` | LAMMPS Data | No | Yes | Yes
| `ovito/src/plugins/particles/xyz` | XYZ | No | Yes | Yes
| `ovito/src/plugins/particles/vasp` | POSCAR/XDATCAR | No | Yes | Yes
| `ovito/src/plugins/particles/imd` | IMD | No | Yes | Yes
| `ovito/src/plugins/particles/fhi_aims` | FHI-Aims | No | Yes | Yes
| `ovito/src/plugins/particles/cfg` | CFG | No | Yes | No
| `ovito/src/plugins/particles/pdb` | PDB | No | Yes | No
| `ovito/src/plugins/crystalanalysis` | Crystal Analysis | No | No | Yes
| `ovito/src/core/dataset/importexport` | Calculation Results File | No | No | Yes
| `ovito/src/plugins/povray` | POV-Ray Scene | No | No | Yes

Of the documented importable binary formats,
NetCDF and GSD/HOOMD have formal specifications,
while LAMMPS and POSCAR do not.

VMD actually accepts lots of file formats,
but the documentation only mentions a few;
table as of version 1.9.3.

| Path | Format | Binary | Import | Export
|:-----|:-------|:-------|:-------|:-------
| `vmd/plugins/molfile_plugin/src/binposplugin.c` | AMBER "BINPOS" Trajectory Reader (.binpos) | Yes | Yes | Yes
| `vmd/plugins/molfile_plugin/src/dcdplugin.c` | CHARMM, NAMD, X-PLOR "DCD" Reader/Writer (.dcd) | Yes | Yes | Yes
| `vmd/plugins/molfile_plugin/src/gromacsplugin.C` | Gromacs TRR/XTC Reader (.trr, .xtc) | Yes | Yes | Yes
| `vmd/plugins/molfile_plugin/src/netcdfplugin.c` | AMBER NetCDF Trajectory Reader (.nc) | Yes | Yes | No
| `vmd/plugins/molfile_plugin/src/gromacsplugin.C` | Gromacs TNG Reader (.tng) | Yes | Yes | No
| `vmd/plugins/molfile_plugin/src/h5mdplugin.c` | H5MD Plugin (.h5) | Yes | Yes | No
| `vmd/plugins/molfile_plugin/src/netcdfplugin.c` | MMTK NetCDF Trajectory Reader (.nc) | Yes | Yes | No
| `vmd/plugins/molfile_plugin/src/crdplugin.c` | AMBER "CRD" Trajectory Reader (.crd, .crdbox) | No | Yes | Yes
| `vmd/plugins/molfile_plugin/src/lammpsplugin.c` | LAMMPS Trajectory Reader (.lammpstrj) | No | Yes | Yes
| `vmd/plugins/molfile_plugin/src/xyzplugin.c` | XYZ Trajectory Files (.xyz) | No | Yes | Yes
| `vmd/plugins/molfile_plugin/src/cpmdplugin.c` | CPMD (CPMD Trajectory) Reader (.cpmd) | No | Yes | No
| `vmd/plugins/molfile_plugin/src/dlpolyplugin.c` | DLPOLY HISTORY File Reader (.dlpolyhist) | No | Yes | No
| `vmd/plugins/molfile_plugin/src/carplugin.c` | VASP Trajectories of Ionic Steps (.xml, .outcar, .xdatcar) | No | Yes | No
| `vmd/plugins/molfile_plugin/src/vtfplugin.c` | VTF Trajectory Files (.vtf) | No | Yes | No
| `vmd/plugins/molfile_plugin/src/xsfplugin.C` | XCrySDen, Quantum Espresso XSF/AXSF Trajectory Files (.axsf, .xsf) | No | Yes | No
| `vmd/plugins/molfile_plugin/src/graspplugin.C` | GRASP Surface File Reader (.grasp, .srf) | Yes | No | Yes
| `vmd/plugins/molfile_plugin/src/msmsplugin.C` | MSMS Surface File Reader (.face, .vert) | No | No | Yes
| `vmd/plugins/molfile_plugin/src/raster3dplugin.C` | Raster3D Scene Reader (.r3d) | No | No | Yes
| `vmd/plugins/molfile_plugin/src/stlplugin.C` | STL Solid Model Triangulated Geometry Files (.stl) | No | No | Yes

Of the documented importable binary trajectory formats,
NetCDF (AMBER) is the only one that OVITO supports as well.
Textual formats that overlap are LAMMPS, XYZ and POSCAR/XDATCAR.

Notably XYZ is the simplest, but does not support custom radii.
However OVITO has a nonstandard extension
that allows its XYZ importer to support pretty much any properties.

Bonds are based on distance search or protein databank properties.
It is impossible to manually bond particles in either program.
This could be worked around by exporting phantom particles
that reside inside real particles and
visualizing the bonds as displacements of said phantom particles.
However real and phantom particles would be indistinguishable and
bloat the data file.

### Program Options

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction of runs.

The following incomplete table lists the options for `bmm-dem`.

| Key | Value | Meaning
|:----|:------|:--------
| `--ncellx` | Nonnegative Integer below `BMM_MCELL` | Number of horizontal cells.
| `--ncelly` | Nonnegative Integer below `BMM_MCELL` | Number of vertical cells.
| `--nbin` | Nonnegative Integer below `BMM_MBIN` | Number of histogram bins.
| `--npart` | Nonnegative Integer below `BMM_MPART` | Number of particles.
| `--nstep` | Nonnegative Integer below `BMM_MSTEP` | Number of simulation steps.

The following table lists the options for `bmm-filter`.

| Key | Value | Meaning
|:----|:------|:--------
| `--mode` | `blacklist` or `whitelist` | Pass or stop all messages.
| `--pass` | Message Name | Pass a certain message.
| `--stop` | Message Name | Stop a certain message.
| `--verbose` | Truth Value | Print statistics at the end.

The following incomplete table lists the options for `bmm-sdl`.

| Key | Value | Meaning
|:----|:------|:--------
| `--width` | Positive Integer | Default window width.
| `--height` | Positive Integer | Default window height.
| `--fps` | Positive Integer below `1000` | Visual frame rate.
| `--ms` | Small Power of Two | Multisample anti-aliasing factor.

### Building a Pipeline

The following stutters or chokes the simulation.

    ./bmm-dem --npart 3 --nstep 20 | \
    ./bmm-sdl --fps 10

The following does not stutter, but chokes the simulation.
However lots of system time may be used.

    stdbuf -o 0 ./bmm-dem --npart 3 --nstep 20 | \
    stdbuf -i 0 ./bmm-sdl --fps 10

The following stutters, but does not choke the simulation.
However lots of memory may be used.

    ./bmm-dem --npart 3 --nstep 20 | \
    sponge | \
    ./bmm-sdl --fps 10

The following does not stutter or choke the simulation.
However lots of memory and system time may be used.

    stdbuf -o 0 ./bmm-dem --npart 3 --nstep 20 | \
    stdbuf -i 0 -o 0 sponge | \
    stdbuf -i 0 ./bmm-sdl --fps 10

These things need better explanations.

### Implementation Details

Here be notes.

#### Error Messages

Since printing error messages taints library procedures and
`errno` and `strerror` cannot carry dynamic error information,
it could be useful to have another mechanism.
The proposal in the `tlerr` translation unit could work.

#### Allocation

Note that `T xs[N][N]` as indexed with `xs[i][j]`
has distinctly different memory access characteristics
from `T xs[N * N]` as indexed with `xs[i + j * n]`.
The same applies to `T xs[N * N]` and `T* xs[N * N]`.
Profile the cache behavior of each one to find the optimal solution.

#### Double Buffering

Some numerical operations need a copy of the universe.
This is provided via double buffering,
but operations that rely on it are obliged to copy the entire universe over
to guarantee no stale data is left behind when the buffers are swapped.

#### Critical Failure

In case a message writer dies in the middle of a message,
the receiver has no way to detect and correct for this.
This is by design.

Failure also occurs if a consumer is fed a message it does not care about.
This could be mitigated with the following reception system,
but that might not be worth the effort (both programming and computational).

    typedef bool (* bmm_recv_proc)(unsigned char*, void*);
    struct bmm_recv_list {
      size_t n;
      bmm_recv_proc procs[BMM_PROC_MAX];
    };
    struct bmm_recv_mask {
      struct bmm_recv_list lists[BMM_MMSG];
    };
    bool bmm_recv_reg(struct bmm_recv_mask*, enum bmm_msg, bmm_recv_proc);
    bool bmm_recv_unreg(struct bmm_recv_mask*, enum bmm_msg, bmm_recv_proc);
    bool bmm_recv_unreg_for(struct bmm_recv_mask*, enum bmm_msg);
    bool bmm_recv_unreg_all(struct bmm_recv_mask*);
    bool bmm_recv(struct bmm_recv_mask const*);

#### Interesting Idea

Particle speed distribution could be calculated and
then partitioned into slow and fast regions.
Only the slow region would participate in the linked cell algorithm.
This did not work due to nonlocality:
transfer of kinetic energy quickly changed the distribution.

#### Pair Symmetry

Currently links are stored in a particle list
each of which with a list of indices representing links *to*.
Another data structure option would be a freestanding index pair list,
where the first index would always be strictly less than the second.

Forces between pairs could be worked out asymmetrically,
which would cut the Moore neighborhood sizes and iterations into half.
However this would have some downsides too.
First, finding whether a particle is linked to any other particle
would require $n g$ instead of $g$ operations
for $n$ particles and $g$ amortized group size.
Second, interaction calculations would require almost twice as much code.

Observe static structure time complexities.
Assume 0-based indexing.
Lists must store a "free list" as an unsorted array to manage allocation.
Trees must use "bit pairs" to represent pointers.

| Data Structure | Insertion | Deletion  | Access    | Search    | Space
|:---------------|:----------|:----------|:----------|:----------|:------
| Sorted Array   | $n - i$   | $n - i$   | $1$       | $\\log n$ | $k$
| Unsorted Array | $1$       | $1$       | $1$       | $i$       | $k$
| Sorted List    | $i$       | $i$       | $i$       | $2 i$     | $3 k$
| Unsorted List  | $2$       | $i$       | $i$       | $2 i$     | $3 k$
| Sorted Tree    | $\\log n$ | $\\log n$ | $\\log n$ | $\\log n$ | $3 k$
| Unsorted Tree  | $\\log n$ | $n$?      | $n$?      | $n$?      | $3 k$

#### Representation of Time

From the user perspective,
it is easiest to specify timespan $t_1 - t_0$ and time step $dt$.
For simplicity set $t_0 = 0$
so that the timespan coincides with the end time $t_1$.

We have two ways to derive the current time $t$ at step $i$.
Either $t = \\sum_{k = 1}^i dt$ or $t = i dt$.
On paper the two are the same,
but in a numerical simulation they slowly diverge.

In terms of total step number $n = t_1 / dt$
the latter would be $t = (i / n) t_1$
(give or take some $\\pm 1$).

Haskell notation follows.
It is unreasonable to expect `sum (replicate n dx) == n * dx`,
but `sum (replicate n dx) <= n * dx + epsilon` should hold for small `epsilon`.

#### Simulation Stages

First particles are placed randomly without overlapping
until enough consecutive failures happen.
Then gravity is applied and the simulation is run
until energy reaches a low point.
The process is repeated until failures happen immediately.
This introduces a statistical size bias, which has to be taken into account.
The time step may be varied (from large to small here)
as long as it produces a stable simulation.

Beams are formed by building the convex hull of the particle point set and
adding and flipping vertices until we run out.
This is basically Delaunay triangulation with length constraints.

#### Unit Stuff

Wikipedia says the following.

> For homogeneous isotropic materials simple relations exist
> between elastic constants (Young's modulus $Y$,
> shear modulus $G$, bulk modulus $K$, and Poisson's ratio $\\nu$)
> that allow calculating them all as long as two are known:
> $Y = 2 G (1 + \\nu) = 3 K (1 - 2 \\nu)$.

#### Types

Ignoring `const` and other qualifiers, the following basic types exist in C11.

| Category                        | Sign-Agnostic    | Signed                | Unsigned
|:--------------------------------|:-----------------|:----------------------|:---------
| Bytes and Characters            | `char`           | `signed char`         | `unsigned char`
| Pointers                        | `void *`         |                       |  
| Function Pointers               | `void (*)(void)` |                       |  
| Integers                        |                  | `short int`           | `unsigned short int`
| Integers                        |                  | `int`                 | `unsigned int`
| Integers                        |                  | `long int`            | `unsigned long int`
| Integers                        |                  | `long long int`       | `unsigned long long int`
| Floating-Point Numbers          |                  | `float`               |  
| Floating-Point Numbers          |                  | `double`              |  
| Floating-Point Numbers          |                  | `long double`         |  
| Special-Purpose Integers        |                  |                       | `size_t`
| Special-Purpose Integers        |                  | `ptrdiff_t`           |  
| Enumerations                    | `bool`           |                       |  
| Extended Integers               |                  | `int8_t`              | `uint8_t`
| Extended Integers               |                  | `int_least8_t`        | `uint_least8_t`
| Extended Integers               |                  | `int_fast8_t`         | `uint_fast8_t`
| Extended Integers               |                  | `int16_t`             | `uint16_t`
| Extended Integers               |                  | `int_least16_t`       | `uint_least16_t`
| Extended Integers               |                  | `int_fast16_t`        | `uint_fast16_t`
| Extended Integers               |                  | `int32_t`             | `uint32_t`
| Extended Integers               |                  | `int_least32_t`       | `uint_least32_t`
| Extended Integers               |                  | `int_fast32_t`        | `uint_fast32_t`
| Extended Integers               |                  | `int64_t`             | `uint64_t`
| Extended Integers               |                  | `int_least64_t`       | `uint_least64_t`
| Extended Integers               |                  | `int_fast64_t`        | `uint_fast64_t`
| Extended Integers               |                  | `intptr_t`            | `uintptr_t`
| Extended Integers               |                  | `intmax_t`            | `uintmax_t`
| Extended Floating-Point Numbers |                  | `complex float`       |  
| Extended Floating-Point Numbers |                  | `complex double`      |  
| Extended Floating-Point Numbers |                  | `complex long double` |  

The following classes or traites are employed.

    |------------------- poly -------------------|  Polymorphic
    |-------------------- eq --------------------|  Equivalence with =
    |-------------- ord ---------------|            Order with < and >
    |-------------- bnd ---------------|            Bounded with ^ and v
       |------------------ num ------------------|  Numeric with +, *, 0 and 1
       |-------- int --------|                      Integral with ++ and --
       |-- sint --|                                 Signed integral with -
                  |-- uint --|                      Unsigned integral
                             |-------- fp -------|  Floating-point with -, / and %
                             |-- sfp --|            Scalar floating-point
                                       |-- vfp --|  Vector floating-point with ||
    |------------------- mono -------------------|  Monomorphic

### Do These Things

See if better sedimentation algorithms exist.

    http://www.ime.unicamp.br/~martinez/packmol/

Try Poisson disc sedimentation.

    https://www.jasondavies.com/poisson-disc/

Refactor integrators and data structures for particle dynamics
into their own little module.

Standardize messages to make the development of consumers easier.
Concretely: combine NPART and PARTS for example.

Finish the NetCDF adapter.

Finish the Gnuplot adapter.

Finish the realtime visualizer.

Try other integrators.

Consider `int` instead of `size_t` in critical parts
since the undefinedness of overflows may allow some better optimizations.

Annotate with `__attribute__ ((__flatten__, __hot__))`.

#### Programming Conventions

All names consist of tokens that are eight or fewer characters long.
Procedures are prefixed with the `namespace##_` token.
Higher-order procedures that take closures
are suffixed with the `##_cls` token.
Procedures that need to emphasize mutation
are suffixed with the `##_mut` token.
In-parameters do not have prefixes or suffixes.
Out-parameters that are written only use the prefix `o##`.
Out-parameters that are read and written use the prefix `io##`.

#### Number Systems

These may be implemented for various number types.
Note that in curried languages "more constant" parameters come first,
so for example `x / y == div y x`; here the opposite is the case.

    /// The call `from(x)`
    /// returns the value of `x`.
    /// This is analogous to the unary operator `(type)`.

    /// The call `to(x)`
    /// returns the value of `x`.
    /// This is analogous to the unary operator `(type)`.

    /// The call `zero()`
    /// returns zero.
    /// This is analogous to the constant `0`.

    /// The call `add(x, y)`
    /// returns the sum of `x` and `y`.
    /// This is analogous to the binary operator `+`.

    /// The call `neg(x)`
    /// returns the negation of `x`.
    /// This is analogous to the unary operator `-`.

    /// The call `sub(x, y)`
    /// returns the difference of `x` and `y`.
    /// This is analogous to the binary operator `-`.

    /// The call `succ(x)`
    /// returns the successor of `x`.

    /// The call `pred(x)`
    /// returns the predecessor of `x`.

    /// The call `one()`
    /// returns one.
    /// This is analogous to the constant `1`.

    /// The call `mul(x, y)`
    /// returns the product of `x` and `y`.
    /// This is analogous to the binary operator `*`.

    /// The call `ipow(x, n)`
    /// returns the value of `x` raised to the integer power of `n`.

    /// The call `two()`
    /// returns two.
    /// This is analogous to the constant `2`.

    /// The call `recip(x)`
    /// returns the reciprocal of `x`.

    /// The call `div(x, y)`
    /// returns the division of `x` and `y`.
    /// This is analogous to the binary operator `/`.

    /// The call `quott(x, y)`
    /// returns the truncated quotient of `x` and `y`.
    /// This is analogous to the binary operator `/`.

    /// The call `remt(x, y)`
    /// returns the truncated remainder of `x` and `y`.
    /// This is analogous to the binary operator `%`.

    /// The call `quote(x, y)`
    /// returns the Euclidean quotient of `x` and `y`.

    /// The call `reme(x, y)`
    /// returns the Euclidean remainder of `x` and `y`.

    /// The call `pow(x, y)`
    /// returns the value of `x` raised to the power of `y`.

    /// The call `rt(x, y)`
    /// returns the value of `x` dropped to the root of `y`.

More exotic things are polymorphic over two types.

    /// The call `dist(x, y)`
    /// returns the distance between `x` and `y`.

    /// The call `smul(x, a)`
    /// returns the scalar product of `x` and `a`.

    /// The call `norm(x)`
    /// returns the norm of `x`.

    /// The call `conj(x)`
    /// returns the conjugate of `x`.

    /// The call `imul(x, y)`
    /// returns the inner product of `x` and `y`.

Allocated versions follow.

    /// The call `from(oy, x)`
    /// stores into `oy` the value of `x`.
    /// This is analogous to the unary operator `(type)`.

    /// The call `to(oy, x)`
    /// stores into `oy` the value of `x`.
    /// This is analogous to the unary operator `(type)`.

    /// The call `zero(ox)`
    /// stores into `ox` zero.
    /// This is analogous to the constant `0`.

    /// The call `add(oz, x, y)`
    /// stores into `oz` the sum of `x` and `y`.
    /// This is analogous to the binary operator `+`.

    /// The call `neg(oy, x)`
    /// stores into `oy` the negation of `x`.
    /// This is analogous to the unary operator `-`.

    /// The call `sub(oz, x, y)`
    /// stores into `oz` the difference of `x` and `y`.
    /// This is analogous to the binary operator `-`.

    /// The call `succ(oy, x)`
    /// stores into `oy` the successor of `x`.

    /// The call `pred(oy, x)`
    /// stores into `oy` the predecessor of `x`.

    /// The call `one(ox)`
    /// stores into `ox` one.
    /// This is analogous to the constant `1`.

    /// The call `mul(oz, x, y)`
    /// stores into `oz` the product of `x` and `y`.
    /// This is analogous to the binary operator `*`.

    /// The call `ipow(oy, x, n)`
    /// stores into `oy` `x` raised to the integer power of `n`.

    /// The call `two(ox)`
    /// stores into `ox` two.
    /// This is analogous to the constant `2`.

    /// The call `recip(oy, x)`
    /// stores into `oy` the reciprocal of `x`.

    /// The call `div(oz, x, y)`
    /// stores into `oz` the division of `x` and `y`.
    /// This is analogous to the binary operator `/`.

    /// The call `quott(oz, x, y)`
    /// stores into `oz` the truncated quotient of `x` and `y`.
    /// This is analogous to the binary operator `/`.

    /// The call `remt(oz, x, y)`
    /// stores into `oz` the truncated remainder of `x` and `y`.
    /// This is analogous to the binary operator `%`.

    /// The call `quote(oz, x, y)`
    /// stores into `oz` the Euclidean quotient of `x` and `y`.

    /// The call `reme(oz, x, y)`
    /// stores into `oz` the Euclidean remainder of `x` and `y`.

Polymorphism mixes up the following.

    /// The call `dist(oa, x, y)`
    /// stores into `oa` the distance between `x` and `y`.

    /// The call `smul(oy, x, a)`
    /// stores into `oy` the scalar product of `x` and `a`.

    /// The call `norm(oa, x)`
    /// stores into `oa` the norm of `x`.

    /// The call `conj(oy, x)`
    /// stores into `oy` the conjugate of `x`.

    /// The call `imul(oa, x, y)`
    /// stores into `oa` the inner product of `x` and `y`.

Mutating versions follow.

    /// The call `add_mut(iox, y)`
    /// stores into `iox` the sum of `iox` and `y`.
    /// This is analogous to the binary operator `+=`.

    /// The call `neg_mut(iox)`
    /// stores into `iox` the negation of `iox`.

    /// The call `sub_mut(iox, y)`
    /// stores into `iox` the difference of `iox` and `y`.
    /// This is analogous to the binary operator `-=`.

    /// The call `succ_mut(iox)`
    /// stores into `iox` the successor of `iox`.
    /// This is analogous to the unary operator `++`.

    /// The call `pred_mut(iox)`
    /// stores into `iox` the predecessor of `iox`.
    /// This is analogous to the unary operator `--`.

    /// The call `mul_mut(iox, y)`
    /// stores into `iox` the product of `iox` and `y`.
    /// This is analogous to the binary operator `*=`.

    /// The call `recip_mut(iox)`
    /// stores into `iox` the reciprocal of `iox`.

    /// The call `div_mut(iox, y)`
    /// stores into `iox` the division of `iox` and `y`.
    /// This is analogous to the binary operator `/`.

    /// The call `quott_mut(iox, y)`
    /// stores into `iox` the truncated quotient of `iox` and `y`.
    /// This is analogous to the binary operator `/=`.

    /// The call `remt_mut(iox, y)`
    /// stores into `iox` the truncated remainder of `iox` and `y`.
    /// This is analogous to the binary operator `%=`.

    /// The call `quote_mut(iox, y)`
    /// stores into `iox` the Euclidean quotient of `iox` and `y`.

    /// The call `reme_mut(iox, y)`
    /// stores into `iox` the Euclidean remainder of `iox` and `y`.

The other cases are, again, tricky.

    /// The call `smul_mut(iox, a)`
    /// stores into `iox` the scalar product of `iox` and `a`.

    /// The call `conj_mut(iox)`
    /// stores into `iox` the conjugate of `iox`.

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
