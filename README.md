# founder-sequences

Create an optimal segmentation with a given minimum segment length from *m* sequences of length *n* in *O(mn)* time and use various matching strategies to join the segment texts.

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

- A recent version of Clang. C++17 support is required. Building the tools has been tested with Clang 5.0.
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

The tool takes a text file that contains a list of sequence file paths as its input. A FASTA file with short (less than 1 kb) lines may be used instead. It reads the sequences into memory and generates the optimal segmentation.

### Example

    founder_sequences --input=input-list.txt --segment-length-bound=10 --output-segments=segments.txt --output-founders=founders.txt

`input-list.txt` should contain the paths of the sequence files, one path per line. The sequence files should contain one sequence in each file without the terminating newline. The segment length bound specifies the minimum segment length.
