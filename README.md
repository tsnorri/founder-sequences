# founder-sequences

Given a minimum segment length and *m* sequences of length *n* drawn from an alphabet of size σ, create a segmentation in *O(mn log σ)* time and use various matching strategies to join the segment texts to generate founder sequences. Please see [releases](https://github.com/tsnorri/founder-sequences/releases) for pre-built binaries.

## Academic use

If you use the software in an academic setting we kindly ask you to cite the following paper:

    @InProceedings{norri_et_al:LIPIcs:2018:9317,
      author ={Tuukka Norri and Bastien Cazaux and Dmitry Kosolobov and Veli M{\"a}kinen},
      title ={{Minimum Segmentation for Pan-genomic Founder Reconstruction in Linear Time}},
      booktitle ={18th International Workshop on Algorithms in  Bioinformatics (WABI 2018)},
      pages ={15:1--15:15},
      series ={Leibniz International Proceedings in Informatics (LIPIcs)},
      ISBN ={978-3-95977-082-8},
      ISSN ={1868-8969},
      year ={2018},
      volume ={113},
      editor ={Laxmi Parida and Esko Ukkonen},
      publisher ={Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik},
      address ={Dagstuhl, Germany},
      URL ={http://drops.dagstuhl.de/opus/volltexte/2018/9317},
      URN ={urn:nbn:de:0030-drops-93175},
      doi ={10.4230/LIPIcs.WABI.2018.15},
      annote ={Keywords: Pan-genome indexing, founder reconstruction, dynamic programming, positional Burrows-Wheeler transform, range minimum query}
    }

## Build Requirements

- A recent version of Clang. C++17 support is required. Building the tools has been tested with Clang 6.0.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [CMake](http://cmake.org)
- [Boost](http://www.boost.org)

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
