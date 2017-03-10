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

The project is divided into several executables.
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

The program `bmm-sdl` is a real-time visualizer built on top of SDL.
The program `bmm-gp` is a batch visualizer
that produces data files and Gnuplot scripts for them.
Other executables may pop up unexpectedly.

Since passing a large number of arguments to the programs and
keeping track of them is such a hassle,
only long options are used and their values are replicated in the output
for the sake of easy reproduction.

[cfdem]: http://www.cfdem.com/
[liggghts]: https://github.com/CFDEMproject/LIGGGHTS-PUBLIC
[uio]: https://www.uio.no/english/services/it/research/hpc/abel/help/software/ESyS-Particle.html
[esys-particle]: https://launchpad.net/esys-particle
[yade-dem]: https://yade-dem.org/
[yade]: https://github.com/yade/trunk
