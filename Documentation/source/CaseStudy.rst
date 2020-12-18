.. |ech2o| replace:: EcH\ :sub:`2`\ O

Case study
==========

The folder named *CaseStudy*, which is distributed with the the model package,
includes a sample dataset. This dataset is used in this case study to
illustrate the process of creating a database and running for a mountain
watershed. We will use the GIS functionality provided by PCRaster to
assist us in the construction of the database. It is a good idea to have
the PCRaster documentation open in a browser tab and peruse it to learn
more about the commands we will be using in these examples.

Configuration files
-------------------

The configuration files present the main communication interface with
the model. They are plain text files with pairs of
keys and values. The values indicate options, paths to folders, or name
of files that contain information needed by |ech2o|-iso. 
There are two configuration files:

* The *main configuration file* (defaut name: config.ini), called in the execution command of |ech2o|-iso.  The list of keywords in the current version of the main configuration file (v1.7) is shown `here <http://ech2o-iso.readthedocs.io/en/latest/Keywords.html>`_.
* The *tracking configuration file* (default name: configTrck.ini), whose location is defined in main configuration file and read *only if* water isotopes and/or age tracking is activated (keyword ``Tracking`` set to 1). The list of keywords in the current version of the tracking configuration file (v1.0) is shown `here <http://ech2o-iso.readthedocs.io/en/latest/KeywordsTrck.html>`_.
  
The easiest way to set the configuration files is to generate templates that can subsequently be
edited. To generate a configuration file template, navigate to the
``CaseStudy`` folder and use the following command:

::

    ech2o_iso -g config.ini

where the ``-g`` option indicates that we wish to generate a main configuration
file with the name ``config.ini``. 
A tracking configuration file with the default name ``configTrck.ini``
is also generated.

In the configuration files we provide information about the location of the database, the names of the files and also
control the length of the run, the size of the time step, the outputs we
want the model to produce and select a number of options.

Open the file with any text editor. In the ``Folder`` section of the file
make sure the paths to the ``Spatial`` and ``Climate`` folders of the case
study are correct. In these files we will be storing spatial and climate
information. Also make sure the folder where |ech2o|-iso will write the
results (``Results`` folder) exists and the path is correct.

The maps to be read by |ech2o|-iso will be in the PCRaster cross-system format so
make sure ``MapTypes = csf``. Also we will be using tables to initialize the vegetation
state variables so make sure ``Species_State_Variable_Input_Method = tables``.

Turn on the reinfiltration and channel switches (1). We will use the
aerodynamic resistance option as used in the Penman-Monteith equation
(option 0) and the soil resistance as used by Ivanov (option 1).

We will simulate 1 year using daily time steps. We have to provide that
information in seconds. The simulation will start at second 0 and will
end at second 31536000. The simulation time step will be 86400. The
climate input time step and the report intervals will be the same as the
simulation time step (86400 seconds, i.e., daily).

The next few sections is where the maps in the database are associated
with parameters in the model. Give the appropriate file name for each
parameter. A description of each parameter can be found in appendix a.
At this point we should have generated all the needed files.

The report map section is a series of boolean switches (0-1) that turn
on or off the reporting (writing to the results folder) of maps with the
state variables. Turn on (*= 1*) the variables that you like reported. Mind
that writing maps to the disk is an expensive processes in terms of
computer time and space.

Populating the database
-----------------------

Creating a base map and importing the elevation model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we will create the spatial section of the database. Open the
command line and navigate to the ``Spatial`` folder of the case study. The
maps created in this section should be placed in this folder.

The first recommended step in the preparation of the database to run is
to prepare a base map holding information on the geometry of the domain
grid (dimension, resolution, etc). This map can be generated when
importing the digital elevation model (DEM) basemap as explained below.

The easiest way to generate the base maps is to obtain a DEM in ArcInfo
ASCII raster format. needs maps in planar coordinates, with lat-long
coordinates in meters, such as the UTM projection. If the map is
obtained in other projection using degrees a reprojection of the map is
necessary using ArcGIS or any other external tool.

Move to the example folder provided with the package, open the file
named with a text editor and check the metadata header with information
on the geometry of the raster image.

Within the PCRaster environment, type::

    mapattr base.map

to start the interface and crate a new blank base map named ``base.map``. Introduce
the number of rows and columns as indicated in the metadata of the ascii
raster image. Choose the *’scalar’* datatype and the *’small real’* cell
representation. If the projection is UTM you may want to indicate a *’y
increases from bottom to top’* projection. Provide the coordinates for
the x upper left corner and for the y upper left corner and the cell
resolution.

Please, note that the ArcInfo standard provides information for the
lower left corner. You can calculate the value of the upper left y
coordinate by adding to the lower left coordinate the result of
multiplying the number of rows by the resolution.

Once this information is provided, press ’q’ and answer ’y’ to write the
newly created map to the drive. Display the map to check it has the
correct dimensions::

    aguila base.map

This base map will be used to import all other maps and to ensure all
the maps in the database have the exact same geometry. To import the
ArcInfo DEM map into the CSF PCRaster format, type::

    asc2map -a --clone base.map dem.asc DEM.map

This command indicates 1) that we are importing an ascii file named ``dem.asc`` into
the PCRaster format with name ``DEM.map``, 2) that the imported file has Arcinfo
ascii grid format, and 3) that we are cloning the geometry of our
base.map.

Display the map to check it has been correctly imported::

    aguila DEM.map

To display it in 3D you can type::

    aguila -3 DEM.map

These maps will form the core of the database from 
other necessary maps can be derived.

Delineating the drainage network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The drainage network is derived from the DEM using a steepest-descent
algorithm on the 8 neighbor window around each cell. From a PCRaster
environment type the command::

    pcrcalc ldd.map = lddcreate(DEM.map, 1e9,1e9,1e9,1e9)

.. attention:: NOTE TO LINUX USERS
     Please, note that if you are following this tutorial in a linux computer
     you need to place the arguments to ``pcrcalc`` between quotes like

      $ pcrcalc ’ldd.map = lddcreate(DEM.map, 1e9,1e9,1e9,1e9)’

This command instructs PCRaster to calculate the local drainage
direction (ldd) for each cell using the dem (``DEM.map``) and save the drainage
network on a map called ``ldd.map``. The large numbers included as the final four
arguments to the *lddcreate* function are options to remove pits and
core areas (see PCRaster documentation on lddcreate for more details).
Display the results with aguila to visually inspect the drainage
network. You may have to zoom in to see the details of the network.

Pits and outlets are coded with the value 5 in the resulting map. These
cells flow nowhere and are considered flow sinks. There is at least one
sink in each basin (the outlet). Mostly we will want to have a
continuous flow network towards the outlet (unless we are working on a
karst area or similar), so if we see internal flow sinks it may be due
to errors in the DEM that to some extent can be corrected with some of
the functions in PCRaster (see PCRaster documentation for this).

A map of the channels and the width of the channel is provided in the
folder ``Spatial``. Inspect it using aguila and observe that cells with a channel
have a positive number indicating the width of the channel in meters and
cells without a channel should have attribute 0.

The resistance presented by the channel to flow is given by Manning’s
:math:`n` coefficient. Values for Manning’s :math:`n` coefficient needs
ot be provided for each cell where the channel width is larger than 0. A
map of Manning’s :math:`n` values in :math:`sm^{-\frac{1}{3}}` for the
example channel network is provided (``chanmannningn.map``).

The parameter controlling the seepage from the subsurface system to the
channel lets us fine-tune subsurface-channel interactions. A good
starting value for this parameter is 0.02 for the entire channel system.
The larger the value, the more resistance to flow into the channel. We
can produce this map using::

    pcrcalc chanparam.map = chanwidth.map/chanwidth.map * 0.02;

Defining soil and surface properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we will create a set of maps that provide information on
the soil and surface properties. Some of these properties can be derived
from the DEM and for others we will use some simplifying assumptions
about the spatial distribution of the properties.

The slope of the terrain can be obtained directly from the DEM using the
following command

::

    pcrcalc slope.map = slope(DEM.map)

This command will create a map named with the slope (rise over run) of
the basin.

Now we will create a map with constant value 1 that will help us create
maps of soil properties with a spatially uniform distribution.

::

    pcrcalc unit.map = DEM.map/DEM.map

This operation divides the DEM map by itself to produce a map called
holding 1 everywhere in the basin.

Now we use to construct maps of spatially uniform properties

::

    pcrcalc albedo.map = unit.map * 0.3
    pcrcalc emissivity.map = unit.map * 0.98
    pcrcalc soilheatcap.map = unit.map * 2.205e6
    pcrcalc soilthermalK.map = unit.map * 0.2
    pcrcalc dampdepth.map = unit.map * 2
    pcrcalc temp_damp.map = unit.map * 10
    pcrcalc snowmeltCoeff.map = unit.map * 4.1e-8
    pcrcalc randrough.map = unit.map * 0.05
    pcrcalc psi_ae.map = unit.map * 0.2
    pcrcalc BClambda.map = unit.map *  5.3
    pcrcalc KvKh.map = unit.map * 0.4
    pcrcalc theta_r.map = unit.map * 0.05
    pcrcalc Wc.map = unit.map * 0.7
    pcrcalc Wp.map = unit.map * 9

This will create maps of uniform albedo, surface emissivity, soil heat
capacity, soil thermal conductivity, soil depth at which heat exchanges
are negligible, initial soil temperature, snowmelt coefficient, terrain
rugosity, soil air entry pressure, Brooks and Corey :math:`\lambda`
parameter, vertical to horizontal hydraulic conductivity anisotropy
ratio, residual soil moisture, and two soil parameter, Wc, Wp,
respectively, with values equal to the multiplying scalar in the right
side of the expression.

To introduce some spatial variability in the simulation, we will assume
that some geomorphologic sorting of soil particles distributes some key
hydrologic properties throughout the basin. For instance, finer
particles may get washed out of steep upslope areas and be deposited
when water slows down in flatter areas downslope. This may produce
deeper, more porous soils at the valley bottom with lower hydraulic
conductivity than soils located higher up in the hillslopes.

To simulate such a geomorphologic driven distribution we can use a
topographic index that is larger for flatter cells with large
contributing areas (such as valley bottoms) and smaller for steep cells
near the water divide.

::

    pcrcalc topind.map = ln(accuflux(ldd.map,10000)/slope.map)

This expression uses the function *accuflux* to accumulate the area of
the cells (10,000 :math:`m^{2}`) following the drainage direction and
divides it by the map of slopes that we created earlier. The function
*ln* takes the logarithm of the result of the quotient to equalize the
distribution of values, which is highly skewed due to the exponential
distribution of the accumulated areas.

We will assume that this map describes the spatial distribution of soil
depth, porosity and effective hydraulic conductivity. With the help of
some scaling functions we produce the resulting fields for these soil
properties:

::

    pcrcalc depth_soil.map = topind.map 
    	/areaaverage(topind.map,nominal(unit.map))
    pcrcalc Keff.map = 1 / (depth_soil.map * 36000)
    pcrcalc poros.map = 1 / (1 + exp(0.01 * topind.map))

We will set initial conditions for the soil assuming the basin starts
free of snow, with 50% of the pores saturated with water and with a
temperature of 10\ :math:`^{\circ}C` throughout the basin:

::

    pcrcalc swe.map = unit.map * 0 
    pcrcalc Soil_moisture_1.map = poros.map * 0.5
    pcrcalc Soil_moisture_2.map = poros.map * 0.5
    pcrcalc Soil_moisture_3.map = poros.map * 0.5
    pcrcalc soiltemp.map = unit.map * 10
    pcrcalc streamflow.map = unit.map * 0

We will also assume that the first hydraulic layer of the soil is 10 cm
deep (0.1 m). We will also assume that the second hydraulic layer is 10
cm deep. will calculate the depth of the 3rd layer such that the sum of
the three layers equals the soil depth at the pixel.
::

    pcrcalc depth_soil1.map = unit.map * 0.1
    pcrcalc depth_soil2.map = unit.map * 0.1

|ech2o|-iso assumes an exponential root profile: :math:`root(z)=exp(-K_{root}z)`.
Here we chose a value of 10 m\ :sup:`-1` for :math:`K_{root}`,
which, given the depths provided above, results in having ~63%
and ~23% of the roots in the first and second layers, respectively.
It is thus specified in the ``SpeciesParams.tab`` file
(column 37, preceding the ``vegtype`` switch).

Finally, for simplicity we further assume that the bedrock at depth of
the soil is impervious (leakance=0). This parameter varies between 0
(no flow boundary) and 1 (free drainage).
::

    pcrcalc leakance.map = unit.map * 0.0

We will see later we will spin-up the model to equilibrate the initial
conditions for the characteristics and climate of the basin.

Defining vegetation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the sake of simplicity we will assume that there is only one type of
forest homogeneously covering 60% of the basin (proportion of area
covered in each forest patch is specified in file ``SpecsProp.tab``).

The parameters that define the vegetation in the forest is provided in
Table 1.

.. csv-table:: **Sample parameter configuration file for one tree species**
   :header: "Parameter", "Species ID"
   :widths: 20, 20

   SpeciesID , 1
   NPP/GPPRatio , 0.47
   gsmax(ms-1) , 0.009
   CanopyQuantumEffic(gCJ-1) , 0.0000018
   MaxForestAge(yrs) , 290
   OptimalTemp(C) , 18
   MaxTemp(C) , 30
   MinTemp(C) , 0
   FoliageAllocCoef\_a , 2.235
   FoliageAllocCoef\_b , 0.006
   StemAllocCoef\_a , 3.3
   StemAllocCoef\_b , 6.00E-07
   gs\_light\_coeff , 300
   gs\_vpd\_coeff , 0.0019
   gs_psi\_d , 5
   gs_psi\_c , 2
   WiltingPnt , 0.05
   SpecificLeafArea(m2g-1) , 0.003
   SpecificRootArea(m2kg-1) , 0.022
   Crown2StemDRat , 0.25
   TreeShapeParam , 0.4
   WoodDens(gCm-2) , 220000
   Fhdmax , 15
   Fhdmin , 5
   LeafTurnoverRate(s-1) , 8.56E-09
   MaxLeafTurnoverWaterStress(s-1) , 0.000000018
   LeafTurnoverWaterStressParam , 0.2
   MaxLeafTurnoverTempStress(s-1) , 0.000000018
   LeafTurnoverTempStressParam , 0.2
   ColdStressParam(degC) , 1
   RootTurnoverRate(s-1) , 5.34E-09
   MaxCanStorageParam(m) , 0.0000624
   albedo , 0.1
   emissivity , 0.95
   KBeers , 0.55
   CanopyWatEffic(gCm-1) , 800
   Kroot(m-1) , 10
   vegtype , 0
   DeadGrassLeafTurnoverRate(s-1) , 0
   DeadGrassLeafTurnoverTempAdjustment(degC) , 0

The parameters are listed in the order they should appear in the
vegetation confguration file. Make sure you include in the first line
of the header the number of species in the file and the number of
information items per species (2 39). For convenience, the information
in Table 1 is properly formatted in a parameter file
named ``SpeciesParams.tab``, which is is provided in the folder of the case study.

The next step is providing information about the variables for
vegetation. There are two ways to provide this information, through
tables that provide constant variable values for each initial forest
patch and through maps that provide spatially variable values.

The easiest way is to provide first the information using tables and
spin-up the model to provide maps with the distributed variables. Then
restart the model using the maps as initial forest conditions. If we are
using tables we need to provide a map with the spatial distribution of
the types of forest or *patches*. This spatial distribution is done
using integer ID numbers for each patch. In this example we will assume
that only one type of forest exist covering the entire area with ID 1.
We can create the patch map using the *unit.map*:

::

    pcrcalc patches.map = unit.map 

The vegetation variables needed to run the model are the proportion of
canopy coverage, the stem density, the leaf area index, the age, the
total basal area, the species height and the root density of each
species for each path. Each of these variables is contained in an
individual file with the same format.

As mentioned earlier, we will assume that the canopy of vegetation type
1 (the only type) cover 60% of the basin. The canopy coverage file for
this example would be

::

    1 2
    1 0.6

Where the first element in the first line indicate the number of patches
(1), the second element is the number of covers in the patch (1
vegetation type + bare soil = 2). The second line indicates the patch ID
for which this line is providing information (matching the appropriate
ID in the file). the following numbers are the proportion of the patch
covered by canopy for each vegetation type (only 1 in the case). The
proportion of bare soil is calculated internally from this information.
The information for each species must be entered in the same order that
was provided in the table of vegetation parameters including 0.0 if
there is no coverage of a specific species or vegetation type in a given
patch. A file (``SpecsProp.tab``) with this information is included for convenience in
the ``CaseStudy/Spatial`` folder.

The same data structure is used in the files containing information for
the other mandatory vegetation variables, for which files are
conveniently provided: Stem density, leaf area index, stand age, total
stand basal area, effective height and root density.

Climate inputs
~~~~~~~~~~~~~~

Navigate to the ``Climate`` folder of the case study. The maps generated
in this section need to be placed in this folder. A climate zones map
provides the information to spatially distribute the climate time series
and should be created first. In this example we will partition our basin
in ten climate zones following the elevation contours. The easiest way
to create to do that is to reclassify the DEM in ten uniform elevation
zones with unique integer IDs using a classification table (see PCRaster
documentation for formats):

::

    <,1430] 1
    <1430,1580] 2
    <1580,1730] 3
    <1730,1880] 4
    <1880,2030] 5
    <2030,2180] 6
    <2180,2330] 7
    <2330,2480] 8
    <2480,2630] 9
    <2630,> 10

Assuming the previous table is saved under the name we can reclassify
our climate zones map with the following command:

::

    pcrcalc ClimZones.map = lookupnominal(ElevZonesClass.txt,
              ..\Spatial\DEM.map)

In the folder there are a set of example climate files providing climate
information for the 10 time zones at a daily time step during 365 days
starting at the beginning of the water year (October 1st). The
temperature field has been generated assuming a sinusoidal cycle of
temperature. A standard environmental lapse rate has been used to
distribute precipitation with elevation. Precipitation for the
bottom-most zone has been generated randomly to simulate a typical
semiarid climate with precipitation falling during fall and winter. A
linear model has been used to simulate an increase of precipitation with
altitude. The other climate variables have been generated using
polynomial functions to simulate seasonality and are considered to be
spatially uniform. You can import the files to a spreadsheet program
like MSExcel and plot them to inspect the type of climate we are
simulating.

In order to make these files usable for |ech2o|-iso we need to import them
into binary format with the utility provided with |ech2o|-iso. This utility
takes two arguments: the name of the properly formatted ascii file with
the climate information and the desired name for the binary file to be
written.

The following commands will import the climate files and generate the
necessary binary files having the same name as the original text files
but with a *.bin* extension.

::

    asc2c Tavg.txt Tavg.bin
    asc2c Tmin.txt Tmin.bin
    asc2c Tmax.txt Tmax.bin
    asc2c Precip.txt Precip.bin
    asc2c Sdown.txt Sdown.bin
    asc2c Ldown.txt Ldown.bin
    asc2c RH.txt RH.bin
    asc2c windspeed.txt windspeed.bin

If the tracking of deuterium and/or oxygen 18 is activated, the corresponding binary files must be generated

::

    asc2c d2H.txt d2H.bin
    asc2c d18O.txt d18O.bin

Besides, we will introduce some random variability in the precipitation field
using the isohyet map assuming no autocorrelation structure or
directionality of the field. The random fluctuations are produced using
a uniform distribution ranging with a range 0.5-1.5 to simulate
precipitation fluctuations ranging from half to one and half times the
incident local precipitation

::

    pcrcalc isohyet.map = uniform(boolean(..\Spatial\unit.map))+0.5

Reporting results
-----------------

The report time series section is another series of boolean switches
that turn on or off the reporting (writing to the results folder) of
time series for the specified variables. The spatial pixels for which
the time series is produced are these indicated by the map. This map
should contain the value zero or no data everywhere except for the
pixels for which a time series is desired. These pixels should be marked
with an integer ID that will be used to identify the time series in the
resulting output file containing the time series information.

One way to crate the map is by making use of a text file containing
information on the location of the pixels to be monitored. This file
should have one line per pixel to be monitored (probe). Each line
contains the x and y coordinate of the pixel and the pixel
identification number. The file located in the folder of this case study
is an example of such a file with three probes.

Once this file is created you can import it to create a map using the
PCRaster tool. Navigate back to the and transform the information
contained in into a map using

::

    col2map --clone base.map probes.txt Tsmask.map

Running the program
-------------------

**FIRST, MAKE SURE THE FOLDER WHERE ECH2O-iso IS TOLD TO WRITE THE RESULTS
EXISTS. ECH2O-iso WILL NOT CREATE THE FOLDER IF IT DOES NOT EXIST AND WILL
TERMINATE THE RUN WHEN IT ATTEMPTS TO WRITE TO THE NON-EXISTING FOLDER.**

Open the configuration file in a text editor and replace the default
input file names for the soil moisture keys with the correct filenames

::

    Soil_moisture_1 = Soil_moisture_1.map 
    Soil_moisture_2 = Soil_moisture_2.map 
    Soil_moisture_3 = Soil_moisture_3.map 

Once the database is complete and the configuration file correctly set
we are ready to run |ech2o|-iso. This is simply done by navigating to the
folder containing the |ech2o|-iso configuration files and running the following
command:

::

    ech2o_iso config.ini


Where ``config.ini`` stands for the name of the configuration file. Note that this
file and ``configTrck.ini`` can be named in any other way to differentiate different
projects or runs.

After hitting enter you will see the splash screen with the version
number and a report on the pre-processing steps (whether it was able to
successfully read the files and create the components of the model run).

The screen reports information on the water mass (in :math:` m^{3} ` for
the different components of the basin for each time step and information
on the mass balance error (in % of the total input). The mass balance
error should be a very small number (typically :math:`<` 1.0e-10%). If
the number is large or steadily increases as the simulation progresses
it is an indication of some problem in the inputs.

Once the model has finished running you can inspect the results using or
to display the timeseries files or the maps in the results folder.

For instance if you have reported discharge you can display it by typing

::

    aguila OutletDisch.tab

Spatial time series can also be displayed. For instance if you have
reported snow water equivalent maps series for one year at daily
timesteps (365 time steps) you can inspect them with the command

::

    aguila SWE00000.001+365

You can even drape them to the DEM. Assuming you are in the folder:

::

    aguila -3 ..\spatial\DEM.map + SWE00000.001+365

Also a file called is created in the root folder (where file is
located). This file contains summary information on the water balance of
the basin in total volumes of water (:math:`m^{3}`).

Spinning-up the model
---------------------

As you see, the model diverges from the initial conditions provided and
will finish with a very different spatial distribution of the state
variables. Some of the variables will show a declining trend, others
will show an increasing trend rather than a cycle. This indicates that
the state of the model is not in equilibrium with the provided boundary
conditions.

The process of running the model in order to allow it time to achieve a
state of equilibrium is called ’spin-up’. The easiest way is to run the
model long enough or multiple times in a loop with the same forcing
until the different state variables show no significant change with
respect to the previous run.

Probably the simplest way is to create a batch file that will run the
model multiple times using the state-variables from the previous run as
initial conditions for the following run. The first step is configure an
initial run as explained in the previous example, using tables to
initialize the vegetation parameters (``Species_State_Variable_Input_Method = table``) and make sure the required state
variables needed to initialize the model will be reported:

-  Report\_Leaf\_Area\_Index = 1

-  Report\_Veget\_frac = 1

-  Report\_Stem\_Density = 1

-  Report\_Stand\_Age = 1

-  Report\_Root\_Mass = 1

-  Report\_Tree\_Height = 1

-  Report\_Basal\_Area = 1

-  Report\_SWE = 1

-  Report\_Soil\_Water\_Content = 1

-  Report\_Soil\_Temperature = 1

The next step is tell the model that next time the model starts, the
vegetation parameters will not be read from a table but that the
parameters will be given as maps. For this set ``Species_State_Variable_Input_Method = maps``.

When this variable is set to ``maps``, the model expects to find a set of maps
in the spatial information folder with the following names. There has to
be one map per species:

-  lai\_n.map

-  p\_n.map

-  ntr\_n.map

-  age\_n.map

-  root_n.map

-  hgt_n.map

-  bas_n.map

To initiate leaf area index, species proportion in each cell, tree
density, age of stands, root mass, tree height and basal area,
respectively. The value for *n* is the species id within [0, ``NumSpecies`` -1], 
where ``NumSpecies`` is the number of species being simulated.

To run the model in a loop we create a batch file that runs the model,
takes the final state of the basin (as per the reported state
variables), copies them with the right name in the spatial folder for
initialization and runs the model again. Assuming you have set up the
model as in the example of the next section, create a batch file with
the name in the same folder where is located. Type the following
contents into a new file ``spinup.bat``:

::

    @echo off
    set count = 1

    :loop

    echo Running iteration %COUNT%

    start /w ech2o_iso config.ini

    ping -w 1000 1.1.1.1 

    echo finishing and copying files after iteration %COUNT%

    copy /Y .\Results/root0_00.365 .\Spatial/root_0.map
    copy /Y .\Results/p0_00000.365 .\Spatial/p_0.map
    copy /Y .\Results/ntr0_000.365 .\Spatial/ntr_0.map
    copy /Y .\Results/lai0_000.365 .\Spatial/lai_0.map
    copy /Y .\Results/hgt0_000.365 .\Spatial/hgt_0.map
    copy /Y .\Results/bas0_000.365 .\Spatial/bas_0.map
    copy /Y .\Results/age0_000.365 .\Spatial/age_0.map

    copy /Y .\Results/SWE00000.365 .\Spatial/SWE.map
    copy /Y .\Results/SWC1_000.365 .\Spatial/Soil_moisture_1.map
    copy /Y .\Results/SWC2_000.365 .\Spatial/Soil_moisture_2.map
    copy /Y .\Results/SWC3_000.365 .\Spatial/Soil_moisture_3.map
    copy /Y .\Results/Ts000000.365 .\Spatial/soiltemp.map
    copy /Y .\Results/Q0000000.365 .\Spatial/streamflow.map

    type .\Results\lai_0.tab >> .\Results\laiaccum.txt
    type .\Results\NPP_0.tab >> .\Results\NPPaccum.txt
    type .\Results\SoilMoistureAv.tab >> .\Results\SWCaccum.txt
    
    set /A COUNT=%COUNT%+1

    goto loop 
  

Run the batch file by typing ``spinup.bat``. This file will spinup the model until you
stop it pressing . Let the model spin for a period of 5 or 10 years.

For convenience, a linux version of this routine is also given in the file ``spinup.sh``

If you are reporting time series of leaf area index, net primary
production and soil moisture, the batch file will append the results in
a file that contains the time series for the entire spinup period.
Plotting this file in Excel will let us evaluate if the state variables
are equilibrated at the end of the spinup period.


