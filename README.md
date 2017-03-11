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

The main program `bmm` does all the heavy lifting.
In addition to running the discrete element method and
producing two output streams of results,
it parses the command line options, reads the initialization files and
handles signals in a synchronous fashion.
Aside from these effects the program resembles a pure function.

The utility program `bmm-filter` removes certain messages from a stream.
If, for example, `bmm` produces real-time progress reports,
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
As an example, this classification makes `bmm` a producer,
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
| `0aaaaaaaD` | Integers are in big-endian (network order).
| `1aaaaaaaD` | Integers are in little-endian.
| `a0aaaaaaD` | Floating-point numbers are in big-endian (network order).
| `a1aaaaaaD` | Floating-point numbers are in little-endian.
| `aabbaaaaD` | Reserved for other options (`b` is free).
| `aaaa0baaD` | Message body has a fixed size (`b` is free).
| `aaaa10aaDH` | Message body is terminated by a literal `h`.
| `aaaa1100DH` | Message body has a size of `h` bytes (256 B).
| `aaaa1101DHH` | Message body has a size of `h` bytes (64 kiB).
| `aaaa1110DHHHH` | Message body has a size of `h` bytes (4 GiB).
| `aaaa1111DHHHHHHHH` | Message body has a size of `h` bytes (16 EiB).
| `Aeeee0000` | Initialization options (`e` is free).
| `A0000eeee` | Runtime information (`e` is free).

To summarize the table informally,
each message is prefixed by two bytes,
the first of which sets various flags and
the second of which contains the message number.
Additionally, those messages whose bodies may vary in size
reserve the first few bytes of their bodies
to signal the size of the rest of the body.

### Program Options

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction of runs.

The following table lists all the basic options of `bmm`.

| Option | Meaning
|:-------|:--------
| `--npart` | Number of particles.
| `--nstep` | Number of simulation steps.
| `--magic` | Magic happens.

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
