.. |ech2o| replace:: EcH\ :sub:`2`\ O

Installation: compiling the code
================================

The following page contains information on how to compile the latest version of |ech2o|-iso, whose executable is ``ech2o_iso``.  
    
The current version of |ech2o|-iso does not have a configure script. The Makefile has been generated for the gnu c++ compiler and does not check for dependencies. MINGW or CYGWIN is necessary to compile in Windows. 

|ech2o|-iso uses PCRASTER as data pre- and post-processor. Please install PCRASTER free of charge from `here <http://pcraster.geo.uu.nl/downloads/latest-release/>`_.
Note that PCRaster is not natively supported on Mac architectures, but can (theoritically) be built from source, see `here <http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/pcraster/build.html>`_.


1. Dependencies
^^^^^^^^^^^^^^^^

* Change directory to your workspace and clone the latest version of the source files from the git repository:

::

  $ git clone https://bitbucket.org/sylka/ech2o_iso.git

* Install the armadillo development files (version 7), either compiling and installing from source or from the package manager of your Linux distribution.

* Precompiled versions of the libcsf dependency for Linux, Windows and Mac are included in the ``lib`` folder. The compilation was carried assuming little endian 64 bit architectures.

  If the linker complains, the library may need to be compiled for your system. Please, clone the source code from 
    
::
   
   $ git clone https://bitbucket.org/maneta/libcsf.git
   

or download from
   
::
   
   $  https://sourceforge.net/p/pcraster/rasterformat/ci/master/tree/
   
and compile from source. Then replace the old libcsf64 library in the ``lib`` directory with the newly compiled library. Make sure you change the name of the new library so it has the same name as the old one. 
   

2. Making ``ech2o_iso``
^^^^^^^^^^^^^^^^^^^^^^^

*  Change to the ``Release-*`` folder within the source folder, where ``*`` is your OS type: Linux, Mac or Windows.

* If compiling for Windows, edit the objects.mk file and substitute item ``-lcsf64`` for ``-llibcsf64`` so that ``make`` will link against the correct static library. Save and close the editor

* If compiling for Mac, edit the objects.mk file and substitute item ``-lcsf64`` for ``-lcsfosx`` so that ``make`` will link against the correct static library. Save and close the editor

* from the command line type ``make`` to make the source.

3. Making ``asc2c``
^^^^^^^^^^^^^^^^^^^^

* Open a command line terminal 
 
* Change directory to your workspace and clone the latest version of the source files from the git repository:

::

   $ git clone https://bitbucket.org/maneta/asc2c.git

* Change directory into the source folder and type ``make`` to make the code. 


4. Contact
----------

If you need assistance compiling the source, contact marco.maneta@umontana.edu, or sylvain.kuppel@abdn.ac.uk.

If you find this documentation to be incomplete, please file a ticket in the appropriate issue tracker:

* ech2o_iso compilation issues:  https://bitbucket.org/sylka/ech2o_iso/issues
* asc2c compilation issues:  https://bitbucket.org/maneta/asc2c/issues
  