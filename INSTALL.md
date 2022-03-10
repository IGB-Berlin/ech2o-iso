# Installation: compiling the code

The following page contains information on how to compile the latest
version of EcH<sub>2</sub>O-iso, whose executable is `ech2o`.

The current version of EcH<sub>2</sub>O-iso does not have a configure
script. The Makefile has been generated for the gnu c++ compiler and
does not check for dependencies. MINGW or CYGWIN is necessary to compile
in Windows.

EcH<sub>2</sub>O-iso uses PCRASTER as data pre- and post-processor.
Please install PCRASTER free of charge from
[here](https://pcraster.geo.uu.nl/downloads/latest-release/).

## 1\. Dependencies

  - Change directory to your workspace and clone the latest version of
    the source files from the git repository:

<!-- end list -->

    $ git clone http://bitbucket.igb-berlin.de:7990/scm/~ech2o/ech2o_iso.git

  - Install the armadillo development files, either compiling and
    installing from source or from the package manager of your Linux
    distribution.

  - Precompiled versions of the libcsf dependency for Linux, Windows and
    Mac are included in the `lib` folder. The compilation was carried
    assuming little endian 64 bit architectures.
    
    If the linker complains, the library may need to be compiled for
    your system. Please, clone the source code from

<!-- end list -->

    $ git clone https://bitbucket.org/maneta/libcsf.git

or download from

    $  https://sourceforge.net/p/pcraster/rasterformat/ci/master/tree/

and compile from source. Then replace the old libcsf64 library in the
`lib` directory with the newly compiled library. Make sure you change
the name of the new library so it has the same name as the old one.

## 2\. Making `ech2o_iso`

  - Change to the `Release` folder within the source folder
  - If compiling for Windows, edit the objects.mk file and substitute
    item `-lcsf64` for `-llibcsf64` so that `make` will link against the
    correct static library. Save and close the editor
  - If compiling for Mac, edit the objects.mk file and substitute item
    `-lcsf64` for `-lcsfosx` so that `make` will link against the
    correct static library. Save and close the editor
  - from the command line type `make` to make the source.

## 3\. Making `asc2c`

  - Open a command line terminal
  - Change directory to your workspace and clone the latest version of
    the source files from the git repository:

<!-- end list -->

    $ git clone https://bitbucket.org/maneta/asc2c.git

  - Change directory into the source folder and type `make` to make the
    code.

### 4\. Contact

If you need assistance compiling the source, contact
<marco.maneta@umontana.edu>

If you find this documentation to be incomplete, please file a ticket in
the appropriate issue tracker:

  - asc2c compilation issues:
    <https://bitbucket.org/maneta/asc2c/issues>
