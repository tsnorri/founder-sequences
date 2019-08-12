# founder-sequences

Given a minimum segment length and *m* sequences of length *n* drawn from an alphabet of size σ, create a segmentation in *O(mn log σ)* time and use various matching strategies to join the segment texts to generate founder sequences. Please see [releases](https://github.com/tsnorri/founder-sequences/releases) for pre-built binaries.

## Academic use

If you use the software in an academic setting we kindly ask you to cite the following paper:

```TeX
@article{NorriCKM19,
  author    = {Tuukka Norri and
               Bastien Cazaux and
               Dmitry Kosolobov and
               Veli M{\"{a}}kinen},
  title     = {Linear time minimum segmentation enables scalable founder reconstruction},
  journal   = {Algorithms for Molecular Biology},
  volume    = {14},
  number    = {1},
  pages     = {12:1--12:15},
  year      = {2019},
  url       = {https://doi.org/10.1186/s13015-019-0147-6},
  doi       = {10.1186/s13015-019-0147-6},
  timestamp = {Mon, 29 Jul 2019 15:58:48 +0200}
}
```

## Build/Runtime Requirements

On Linux the following libraries are required:

* [zlib](https://zlib.net/)
* libpthread
* [libbsd](https://libbsd.freedesktop.org/)

## Build Requirements

- A recent version of [Clang](https://clang.llvm.org/) and libclang. C++17 support is required. Building the tools has been tested with Clang 6.0.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [CMake](https://cmake.org)
- [Boost](https://www.boost.org)
- On Linux, libdispatch dependencies including [Ninja](https://ninja-build.org) and [SystemTap](https://sourceware.org/systemtap/) development files are required.

For installing libdispatch dependencies, the package list in the [Building and installing for Linux](https://github.com/apple/swift-corelibs-libdispatch/blob/master/INSTALL.md#building-and-installing-for-linux) section in the libdispatch installation guide can be helpful. On Linux, libdispatch itself is built as part of our build process, so it does not need to be installed. On macOS, the operating system libraries are used instead.

## Building

### Short version

1. `git clone --recursive https://github.com/tsnorri/founder-sequences.git`
2. `cd founder-sequences`
3. `cp linux-static.local.mk local.mk`
4. Edit local.mk.
5. `make -j12`

### Long version

1. Clone the repository with `git clone --recursive https://github.com/tsnorri/founder-sequences.git`.
2. Change the working directory with `cd founder-sequences`.
3. Create the file `local.mk`. `linux-static.local.mk` is provided as an example and may be copied with `cp linux-static.local.mk local.mk`
4. Edit `local.mk` in the repository root to override build variables. Useful variables include `CC`, `CXX`, and `GENGETOPT` for C and C++ compilers and gengetopt respectively. `BOOST_ROOT` is used to determine the location of Boost headers and libraries. `BOOST_LIBS` and `LIBDISPATCH_LIBS` are passed to the linker. See `common.mk` for additional variables.
5. Run make with a suitable numer of parallel jobs, e.g. `make -j12`

Useful make targets include:

<dl>
<dt>all</dt>
<dd>Build everything</dd>
<dt>clean</dt>
<dd>Remove build products except for dependencies (in the <code>lib</code> folder).</dd>
<dt>clean-all</dt>
<dd>Remove all build products.</dd>
</dl>

## Running

The package contains `founder_sequences` as well as some auxiliary tools.

### founder\_sequences

Takes a text file that contains a list of sequence file paths as its input. A FASTA file with short (less than 1 kb) lines may be used instead. It generates a segmentation with substrings not shorter than the value given with `--segment-length-bound`. It then proceeds to join the segments with the joining method specified with `--segment-joining` and writes the founder sequences to the path given with `--output-founders` one sequence per line. In addition, the segments may be written to a separate file with `--output-segments`.

#### Example

    founder_sequences --input=input-list.txt --segment-length-bound=10 --output-segments=segments.txt --output-founders=founders.txt

`input-list.txt` should contain the paths of the sequence files, one path per line. The sequence files should contain one sequence in each file without the terminating newline. The segment length bound specifies the minimum segment length.

### remove\_identity\_columns

Reads the aligned texts file paths given from a given list. Outputs the reduced texts to files created in the current directory. The identity columns will be listed as a sequence of zeros and ones (indicates identity) to the standard output.

### insert\_identity\_columns

Given a set of founder sequences, a reference sequence and a list of identity columns, outputs the founder sequences with the identity columns included.

### match\_founder\_sequences

Matches sequences to founder sequences and outputs statistics. Uses a greedy algorithm to find the longest match in the set of founders.
