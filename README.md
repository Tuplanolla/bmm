# Brittle Matter Matters

This repository contains an implementation
of a force-based discrete element method
as applied to shearing brittle surfaces.
While the code seems to produce correct results and run reasonably fast,
it was originally intended to accomplish a single task.
If you would prefer something that is actually useful, see

* [University of Oslo][uio]'s [ESyS-Particle][esys-particle] under OSL-3.0,
* [Yade Project][yade-dem]'s [Yade][yade] under GPL-2.0 or
* [CFDEM][cfdem]'s [LIGGGHTS][liggghts] under GPL-2.0.

## Overview

Brittle Matter Matters was written by Sampsa "Tuplanolla" Kiiskinen
to support his master's thesis on complex material physics.
The project was written between 2017-03-09 and 2018-01-11.

## License

Brittle Matter Matters is free software and as such
licensed under the GNU General Public License version 3 or later.
The full license can be found in the `LICENSE` file that
resides in the same directory as this file.
In short, copies and derivative works are permitted
as long as they use a compatible license.

## Notes and Other Nonsense

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
Of the two

* VMD
* OVITO

the former is too limited and
the latter does not support bonds without hacks.

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

### Critical Failure

In case a message writer dies in the middle of a message,
the receiver has no way to detect and correct for this.
This is by design, because speeed.

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

### Representation of Time

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

### Programming Conventions

All names consist of tokens that are eight or fewer characters long.
Procedures are prefixed with the `namespace##_` token.
Higher-order procedures that take closures
are suffixed with the `##_cls` token.
Procedures that need to emphasize mutation
are suffixed with the `##_mut` token.
In-parameters do not have prefixes or suffixes.
Out-parameters that are written only use the prefix `o##`.
Out-parameters that are read and written use the prefix `io##`.

### Things I Never Got Around Doing

* Refactor integrators and data structures for particle dynamics
  into their own little module.
* Standardize messages to make the development of consumers easier.
* Finish the NetCDF adapter.
* Finish the Gnuplot adapter.
* Finish the realtime visualizer.
* Annotate with `__attribute__ ((__flatten__, __hot__))`.
* Consider `int` instead of `size_t` in critical parts since the undefinedness
  of signed overflows may allow some better optimizations.

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
