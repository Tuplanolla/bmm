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

### Do This

TODO See if better sedimentation algorithms exist.

    http://www.ime.unicamp.br/~martinez/packmol/

### Program Options

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction of runs.

The following incomplete table lists the options for `bmm-dem`.

| Key | Value | Meaning
|:----|:------|:--------
| `--ncellx` | Nonnegative Integer below `BMM_NCELL` | Number of horizontal cells.
| `--ncelly` | Nonnegative Integer below `BMM_NCELL` | Number of vertical cells.
| `--nbin` | Nonnegative Integer below `BMM_NBIN` | Number of histogram bins.
| `--npart` | Nonnegative Integer below `BMM_NPART` | Number of particles.
| `--nstep` | Nonnegative Integer below `BMM_NSTEP` | Number of simulation steps.

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
      struct bmm_recv_list lists[BMM_NMSG];
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

#### Notational Conventions

Heed these.

| Symbol        | Context  | Meaning
|:--------------|:---------|:--------
| $a$           | Physics  | Acceleration
| $e$           | Physics  | Energy
| $f$           | Physics  | Force
| $i$, $j$, $k$ | Software | Storage Index
| $j$           | Physics  | Jerk
| $j$           | Physics  | Moment of Inertia
| $k$           | Software | Storage Capacity
| $l$           | Software | Unique Label
| $m$           | Physics  | Mass
| $n$           | Software | Number
| $r$           | Physics  | Radius
| $s$           | Physics  | Jounce
| $t$           | Physics  | Time
| $v$           | Physics  | Velocity
| $x$           | Physics  | Position
| $y$           | Physics  | Young Modulus

| Symbol      | Meaning
|:------------|:--------
| $\\alpha$   | Angular Acceleration
| $\\epsilon$ | Coefficient of Restitution
| $\\phi$     | Angle
| $\\mu$      | Coefficient of Friction
| $\\nu$      | Poisson Ratio
| $\\xi$      | Compression
| $\\tau$     | Torque
| $\\omega$   | Angular Velocity

#### Partitioning Moore Neighborhoods

Insert a hyperplane into the lattice in such a way that it

* it partitions the neighboring points into two,
* its normal vector's components are positive and increasing and
* it maximizes the minimum distance between each lattice point and the plane.

This construction should be unique.

It seems that approximately right results can be obtained with

* $(1)$ for $1$ dimensions,
* $(1, 2)$ for $2$ dimensions,
* $(1, 2, 4)$ for $3$ dimensions and
* maybe $(1, 2, 4, 8)$ for $4$ dimensions.

The minimum distance between each lattice point (displacement $a$) and
the plane (with normal vector $n$) would then be

$$
|r| = \\frac{|a \\cdot n|}{|n|}
= \\frac 1{\\sqrt{\\sum_{k = 1}^d 2^{2 (k - 1)}}}.
$$

The plane for the normal $n$ has the equation $\\sum_{k = 1}^d n_k x_k = 0$.
Given the lattice indices $i$ from $1$ to $3^d$ and
excluding the origin $(1 + 3^d) / 2$,
the lattice points $x$ follow the lexicographic ordering relation $L$
as $x_j = L_d(i_j) - 2$.
Now the upper partition satisfies $\\sum_{k = 1}^d n_k x_k > 0$ or,
with indices, $\\sum_{k = 1}^d n_k (L_d(i_j)_k - 2) > 0$ or,
considering the rest of the imposed constraints, $i_j > (1 + 3^d) / 2$.

#### Spherical Coordinates

Here is a sane formulation of the spherical coordinate system pictured below.

                x
                 d

    d - 1       ^
    -----       |     +
     | |  x     |  r /|
     | |   i   ,|.  / |
    i = 1   ,-' | `/. |  t
         ,-'    | /\ `|.  d
      ,-'       |/  | | `-.
    -:          +-----+    :-
      `-.               ,-'
         `-.         ,-'
            `-.   ,-'
               `-'

Let $0 \\le \\theta_n < \\twopi$ for some $1 < n \\le d$ and
let $0 \\le \\theta_k < \\twopi / 2$ for all the rest $1 < k \\le d$
satisfying $k \\ne n$.
The components obey the following relations.

$$
\\begin{array}{cccccccc}
x_1       & = & r & (\\sin \\theta_1) & \\cos \\theta_2 & \\cos \\theta_3 & \\dotsb & \\cos \\theta_{d - 1} & \\cos \\theta_d \\\\
x_2       & = & r &                   & \\sin \\theta_2 & \\cos \\theta_3 & \\dotsb & \\cos \\theta_{d - 1} & \\cos \\theta_d \\\\
x_3       & = & r &                   &                 & \\sin \\theta_3 & \\dotsb & \\cos \\theta_{d - 1} & \\cos \\theta_d \\\\
\\vdots   &   &   &                   &                 &                 & \\ddots &                       &                 \\\\
x_{d - 1} & = & r &                   &                 &                 &         & \\sin \\theta_{d - 1} & \\cos \\theta_d \\\\
x_d       & = & r &                   &                 &                 &         &                       & \\sin \\theta_d
\\end{array}
$$

Parenthesized is a suggestive extra factor.
Define $\\theta_1 = \\twopi / 4$ to obtain

$$
x_n = r \\sin \\theta_n \\prod_{k = n}^{d - 1} \\cos \\theta_k
$$

for all $1 \\le n \\le d$.

Flip the angles with $\\theta' = \\twopi / 4 - \\theta$ and
reverse the components with $n' = d - n + 1$.
Also shift the indices of the angles from $1 < n \\le d$ to $1 \\le n < d$.
The components now obey the following relations.

$$
\\begin{array}{cccccccc}
x_1       & = & r & \\cos \\theta_1 &                 &         &                       &                       &                   \\\\
x_2       & = & r & \\sin \\theta_1 & \\cos \\theta_2 &         &                       &                       &                   \\\\
\\vdots   &   &   &                 &                 & \\ddots &                       &                       &                   \\\\
x_{d - 2} & = & r & \\sin \\theta_1 & \\sin \\theta_2 & \\dotsb & \\cos \\theta_{d - 2} &                       &                   \\\\
x_{d - 1} & = & r & \\sin \\theta_1 & \\sin \\theta_2 & \\dotsb & \\sin \\theta_{d - 2} & \\cos \\theta_{d - 1} &                   \\\\
x_d       & = & r & \\sin \\theta_1 & \\sin \\theta_2 & \\dotsb & \\sin \\theta_{d - 2} & \\sin \\theta_{d - 1} & (\\cos \\theta_d)
\\end{array}
$$

This time around define $\\theta_d = 0$ to obtain

$$
x_n = r \\cos \\theta_n \\prod_{k = 1}^{n - 1} \\sin \\theta_k
$$

for all $1 \\le n \\le d$.

Choose the favorable formulation based on the direction of folding.

Acknowledge Tero Harjupatana and Wikipedia for providing good ideas.

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

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
