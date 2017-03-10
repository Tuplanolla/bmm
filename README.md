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

## Draft Stuff

### The Plan

The project is divided into several programs.
The main program `bmm` does all the heavy lifting.
It reads initialization files, runs the discrete element method,
handles signals in a synchronous fashion and produces two output streams.
Aside from these effects the program resembles a pure function.

The primary high-bandwidth output stream `stdout`
consists of a sequence of events that follow a strict messaging protocol.
It is used to pass on all the necessary information of the system and
is thus intended to be piped into another process for visualization or
saved into a file for later analysis.
The secondary low-bandwidth stream `stderr`
contains diagnostics and error messages in an unstructured format.
It is meant to be read as is by the user.

If multiple threads are used, some care must be taken to ensure data coherence.
For `stdout` this can be accomplished
with a simple lock or a size-restricted atomic writing policy.
For `stderr` data coherence is less imporant,
so line-buffering should be sufficient.

The utility program `bmm-filter` removes certain messages from a stream.
If, for example, `bmm` produces real-time progress reports,
they can be stripped out before the results are saved into a file.

The analysis programs `bmm-sdl` and `bmm-gp`
draw visualizations in different ways.
The former is a real-time visualizer built on top of SDL and
the latter is a batch visualizer
that produces data files and Gnuplot scripts for them.
Other programs may pop up unexpectedly.

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction of runs.

### Messaging Protocol

To attain maximal throughput,
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

At the lowest level of abstraction
messages are defined in terms of bits and bytes (that are assumed to be octets)
instead of other abstract concepts.
The following table contains all the possible bit patterns and their meanings.
In the table repeated lowercase letters denote bit wildcards and
repeated uppercase letters denote byte wildcards.
Anything not mentioned in the table is assumed to be invalid and
sending or attempting to interpret such a pattern is a protocol violation.

| Bit Pattern | Meaning
|:------------|:--------
| `0xxxxxxx Y` | Integers are in big-endian (network order).
| `1xxxxxxx Y` | Integers are in little-endian.
| `x0xxxxxx Y` | Floating-point numbers are in big-endian (network order).
| `x1xxxxxx Y` | Floating-point numbers are in little-endian.
| `xxuuxxxx Y` | Reserved for other options (`u` is free).
| `xxxx0uxx Y` | Message body has a fixed size (`u` is free).
| `xxxx10xx Y W` | Message body is terminated by a literal `w`.
| `xxxx1100 Y W` | Message body has a size of `w` bytes.
| `xxxx1101 Y W W` | Message body has a size of `w` bytes.
| `xxxx1110 Y W W W W` | Message body has a size of `w` bytes.
| `xxxx1111 Y W W W W W W W W` | Message body has a size of `w` bytes.
| `X 0000vvvv` | Message number `v`.

To summarize the table informally,
each message is prefixed by two bytes,
the first of which sets various flags and
the second of which contains the message number.
Additionally, those messages whose bodies may vary in size
reserve the first four bytes of their bodies
to signal the size of the rest of the body.

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
