# founder-sequences

Create an optimal segmentation with a given minimum segment length from *m* sequences of length *m* in *O(mn)* time and use bipartite matching to join the segment texts.

**Note:** [git-remote-hg](https://github.com/felipec/git-remote-hg) is required to clone this repository.
Please use a command similar to the following ones to clone:

    GIT_ALLOW_PROTOCOL=hg:https git clone --recursive https://github.com/tsnorri/founder-sequences.git
    GIT_ALLOW_PROTOCOL=hg:https:ssh git clone --recursive git@github.com:tsnorri/founder-sequences.git

## Build Requirements

- A recent version of Clang. C++17 support is required. Building the tools has been tested with Clang 5.0.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [CMake](http://cmake.org)
- [Boost](http://www.boost.org)

## Building

### Short version

1. `GIT_ALLOW_PROTOCOL=hg:https git clone --recursive https://github.com/tsnorri/founder-sequences.git`
2. `cd founder-sequences`
3. `cp linux-static.local.mk local.mk`
4. Edit local.mk.
5. `make -j12`

### Long version

1. Clone the repository with `GIT_ALLOW_PROTOCOL=hg:https git clone --recursive https://github.com/tsnorri/founder-sequences.git`. Since one of the submodules uses Mercurial for version control, `git-remote-hg` needs to be installed and git needs the permission to use fetch the submodule.
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

The tool takes a text file that contains a list of sequence file paths as its input. A FASTA file with short (less than 1 kb) lines may be used instead. It reads the sequences into memory and generates the optimal segmentation.
