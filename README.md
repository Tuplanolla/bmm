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
Afterwards comes the body, which begins with a (possibly empty) type.
The size in the prefix covers the whole body.

    | Header         | Body
    | Flags | Prefix | Type | Data

To summarize the table informally,
each message is prefixed by two nibbles,
the first of which contains user-set flags,
the second of which contains derived flags.
Additionally, for those messages whose bodies may vary in size,
there is an extra prefix to signal the size of the message.
The next byte contains the message type and after that comes the message body.

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
      enum bmm_msg_type type;
      unsigned char body[];
    };

For words of length $n$ bytes there are $n!$ byte ordering schemes.
To enumerate all byte orderings for words of $m = 8$ bytes,
as is the maximal message length,
one would need $b = \lceil\log_2 m!\rceil = 16$ bits.

It is better to assume that byte ordering schemes are built from pair swaps,
so the permutations form a binary tree.
Now $d = \lceil\log_2 n\rceil$ is the tree depth, starting from zero.
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
the message type and flags as usual,
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

TODO Try producing high-quality images and movies with OVITO.

    http://www.ovito.org/manual/usage.import.html
    http://www.ovito.org/manual/usage.export.html

It is a bit picky about file formats according to the documentation and code.

| Path | Format | Import | Export | Binary
|:-----|:-------|:-------|:-------|:-------
| `ovito/src/plugins/particles/lammps` | LAMMPS Dump | Yes | Yes | Yes
| `ovito/src/plugins/particles/lammps` | LAMMPS Data | Yes | Yes | No
| `ovito/src/plugins/particles/xyz` | XYZ | Yes | Yes | No
| `ovito/src/plugins/particles/vasp` | POSCAR / XDATCAR | Yes | Yes | No
| `ovito/src/plugins/particles/imd` | IMD | Yes | Yes | No
| `ovito/src/plugins/particles/fhi_aims` | FHI-Aims | Yes | Yes | No
| `ovito/src/plugins/particles/parcas` | PARCAS | Yes | No | Yes
| `ovito/src/plugins/netcdf` | NetCDF | Yes | No | Yes
| `ovito/src/plugins/particles/gsd` | GSD/HOOMD | Yes | No | Yes
| `ovito/src/plugins/particles/cfg` | CFG | Yes | No | No
| `ovito/src/plugins/particles/pdb` | PDB | Yes | No | No
| `ovito/src/plugins/crystalanalysis` | Crystal Analysis | No | Yes | No
| `ovito/src/core/dataset/importexport` | Calculation Results File | No | Yes | No
| `ovito/src/plugins/povray` | POV-Ray Scene | No | Yes | No

### Program Options

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction of runs.

The following table lists the options for `bmm-dem`.

| Option | Meaning
|:-------|:--------
| `--ncellx` | Number of horizontal cells (at most `BMM_CELL_MAX`).
| `--ncelly` | Number of vertical cells (at most `BMM_CELL_MAX`).
| `--nbin` | Number of histogram bins (at most `BMM_BIN_MAX`).
| `--npart` | Number of particles (at most `BMM_PART_MAX`).
| `--nstep` | Number of simulation steps (at most `BMM_STEP_MAX`).

The following table lists the options for `bmm-sdl`.

| Option | Meaning
|:-------|:--------
| `--width` | Default window width.
| `--height` | Default window height.
| `--fps` | Visual frame rate (between `1` and `1000`).
| `--ms` | Multisample anti-aliasing factor (a small power of two).

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
      struct bmm_recv_list lists[BMM_MSG_MAX];
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

Forces between pairs could be worked out asymmetrically,
which would cut the Moore neighborhood sizes and iterations into half.
However this would prevent one from not having to build neighbor lists and
short-circuiting interaction calculations for fixed particles.

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
