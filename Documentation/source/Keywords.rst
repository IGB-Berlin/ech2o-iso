Configuration (main) file keywords
==================================

**Path definitions**

+--------------------+--------------+-----------------------------------------------------------+
| Keyword            | Type         | Description                                               |
+====================+==============+===========================================================+
| Maps\_Folder       | System Path  | Path to folder with land surface information              |
+--------------------+--------------+-----------------------------------------------------------+
| Clim\_Maps\_Folder | System Path  | Path to folder with Climate information                   |
+--------------------+--------------+-----------------------------------------------------------+
| Output\_Folder     | System Path  | Path to folder where outputs will be written              |
+--------------------+--------------+-----------------------------------------------------------+

**Tracking switch**

+--------------------+--------------+------------------------------------------------------------+
| Keyword            | Type         | Description                                                |
+====================+==============+============================================================+
| Tracking           | option       | | Switch 1/0 to turn water tracking (isotopes and/or ages) |
|                    |              | | on (1) or off (0)                                        |
+--------------------+--------------+------------------------------------------------------------+
| TrackingConfig     | System Path  | Location and name of the tracking configuration file.      |
+--------------------+--------------+------------------------------------------------------------+

**Options**

+-----------------------------------------+---------------------------------------------------------------+
| Keyword                                 | Description                                                   |
+=========================================+===============================================================+
| MapTypes                                | Format of maps, in this version it is *csf* (PCRaster)        |
+-----------------------------------------+---------------------------------------------------------------+
| Species\_State\_Variable\_Input\_Method | | Specifies the input format of the vegetation state          |
|                                         | | variables. Options are *table* or *maps*                    |
+-----------------------------------------+---------------------------------------------------------------+
| Vegetation\_dynamics                    | | Switches between different approaches to vegetation         |
|                                         | | dynamics: 0 assumes constant LAI (equal to initial value),  |
|                                         | | 1 turns on vegetation allocation and growth module to       |
|                                         | | calculate LAI, and 2 allows for externally forcing LAI.     |
+-----------------------------------------+---------------------------------------------------------------+
| TimeSeries\_LAI                         | | if Vegetation\_dynamics = 2, template name of binary files  |
|                                         | | giving LAI dynamics, one for each species.                  |
|                                         | | Files name format: template + '_' + species # + '.bin'      |
+-----------------------------------------+---------------------------------------------------------------+
| Reinfiltration                          | Boolean switch to turn reinfiltration on (1) or off (0)       |
+-----------------------------------------+---------------------------------------------------------------+
| Aerodyn\_resist\_opt                    | | Switches between different aerodynamic resistance           |
|                                         | | formulations. 0: Penman; 1: Thom and Oliver (1977);         |
+-----------------------------------------+---------------------------------------------------------------+
| Soil\_resistance\_opt                   | | Switches between different soil resistance                  |
|                                         | | formulations. 0: No resistance; 1: Passerat de Silans       |
|                                         | | et al. (1989); 2: Sellers et al. (1992); 3: Sakaguchi and   |
|                                         | | Zeng (2009)                                                 |
+-----------------------------------------+---------------------------------------------------------------+

**Time controls**

+---------------------+---------+---------+--------------------------------------------------------------+
| Keyword             | Type    | Unit    | Description                                                  |
+=====================+=========+=========+==============================================================+
| Simul\_start        | Integer | Seconds | | Time of simulation start. In the current version this      |
|                     |         |         | | value must be 0                                            |
+---------------------+---------+---------+--------------------------------------------------------------+
| Simul\_end          | Integer | Seconds | | Time when simulation ends in seconds. This value           |
|                     |         |         | | indicates the total simulated time                         |
+---------------------+---------+---------+--------------------------------------------------------------+
| Simul\_tstep        | Integer | Seconds | | Size of the integration time step                          |
+---------------------+---------+---------+--------------------------------------------------------------+
| Clim\_input\_tstep  | Integer | Seconds | | Time step of climate forcing. Typically it is the same as  |
|                     |         |         | | *Simul\_tstep* but can be larger (i.e. climate inputs are  |
|                     |         |         | | daily but we are using an hourly integration time step).   |
|                     |         |         | | *Clim\_input\_tstep* cannot be smaller than *Simul\_tstep* |
+---------------------+---------+---------+--------------------------------------------------------------+
| Report\_interval    | Integer | Seconds | | Intervals between time series outputs. *Report\_interval*  |  
|                     |         |         | | cannot be smaller than *Simul\_tstep* and typically it is  |
|                     |         |         | | equal to *Simul\_tstep*                                    |
+---------------------+---------+---------+--------------------------------------------------------------+
| ReportMap\_interval | Integer | Seconds | | Intervals between maps outputs. *ReportMap\_interval*      |
|                     |         |         | | cannot be smaller than *Simul\_tstep*                      |
+---------------------+---------+---------+--------------------------------------------------------------+

**Climate information** (maps and binary files must be placed in ``Clim_Maps_Folder``)

+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| Keyword                     | Type                | Unit                 | Description                                           |
+=============================+=====================+======================+=======================================================+
| Snow\_rain\_temp\_threshold | scalar              | :math:`^{\circ}C`    | Air temperature threshold for snow/rain transition    |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| ClimateZones                | Map file name       | :math:`[-]`          | Map identifying the climate zones                     |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| Isohyet\_map                | Map file name       | :math:`[-]`          | | This map allows to redistribute rainfall within a   |
|                             |                     |                      | | climate zone. It is a map with multiplication       |
|                             |                     |                      | | factors for rain in a given pixel. A map containing |
|                             |                     |                      | | 1 over the domain has the effect of overriding this |
|                             |                     |                      | | factors for rain in a given pixel. A map containing |
|                             |                     |                      | | input (does not modify the precipitation input).    |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| Precipitation               | Binary climate file | :math:`ms^{-1}`      | Precipitation input                                   |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| AirTemperature              | Binary climate file | :math:`^{\circ}C`    | Average air temperature                               |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| MaxAirTemp                  | Binary climate file | :math:`^{\circ}C`    | Maximum air temperature                               |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| MinAirTemp                  | Binary climate file | :math:`^{\circ}C`    | Maximum air temperature                               |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| RelativeHumidity            | Binary climate file | :math:`kPa kPa^{-1}` | Relative humidity of the Air                          |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| WindSpeed                   | Binary climate file | :math:`ms^{-1}`      | Wind speed                                            |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| IncomingLongWave            | Binary climate file | :math:`Wm^{-2}`      | Incoming long wave radiation                          |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+
| IncomingShortWave           | Binary climate file | :math:`Wm^{-2}`      | Incoming solar radiation                              |
+-----------------------------+---------------------+----------------------+-------------------------------------------------------+

**Drainage network** (files must be located in ``Maps_Folder``)

+------------------------------+---------------+-------------------+----------------------------------------------------------+
| Keyword                      | Type          | Unit              | Description                                              |
+==============================+===============+===================+==========================================================+
| local\_drain\_direc          | Map file name | :math:`[-]`       | D8 steepest descent ldd                                  |
+------------------------------+---------------+-------------------+----------------------------------------------------------+
| channel\_width               | Map file name | :math:`m`         | | mask with width of channel network. Pixels with no     |
|                              |               |                   | | channel must be 0 or nodata. Positive numbers indicate |
|                              |               |                   | | the width of the channel in the pixel                  |
+------------------------------+---------------+-------------------+----------------------------------------------------------+
| channel\_gw\_transfer\_param | Map file name | :math:`m^{-1}`    | | Coefficient controlling transfers of water from the    |
|                              |               |                   | | subsurface system to the channel                       |
+------------------------------+---------------+-------------------+----------------------------------------------------------+
| mannings\_n                  | Map file name | :math:`sm^{-1/3}` | Manning's n roughness coefficient for channel            |
+------------------------------+---------------+-------------------+----------------------------------------------------------+

**Initial conditions for soil states** (files must be located in ``Maps_Folder``)

+-------------------------+---------------+--------------------+------------------------------------------------------+
| Keyword                 | Type          | Unit               | Description                                          |
+=========================+===============+====================+======================================================+
| Streamflow              | Map file name | :math:`m^3 s^{-1}` | Streamflow                                           |
+-------------------------+---------------+--------------------+------------------------------------------------------+
| snow\_water\_equivalent | Map file name | :math:`m`          | Snow water equivalent                                |
+-------------------------+---------------+--------------------+------------------------------------------------------+
| Soil\_moisture\_1       | Map file name | :math:`m^3 m^{-3}` | Volumetric soil water content for topmost soil layer |
+-------------------------+---------------+--------------------+------------------------------------------------------+
| Soil\_moisture\_2       | Map file name | :math:`m^3 m^{-3}` | Volumetric soil water content for layer 2            |
+-------------------------+---------------+--------------------+------------------------------------------------------+
| Soil\_moisture\_3       | Map file name | :math:`m^3 m^{-3}` | Volumetric soil water content of bottommost layer    |
+-------------------------+---------------+--------------------+------------------------------------------------------+
| Soil\_temperature       | Map file name | :math:`^{\circ}C`  | Soil temperature at boundary of thermal layer        |
+-------------------------+---------------+--------------------+------------------------------------------------------+


**Soil parameters** (files must be located in ``Maps_Folder``)

+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Keyword                        | Type          | Unit                    | Description                                            |
+================================+===============+=========================+========================================================+
| DEM                            | Map file name | :math:`m`               | Digital elevation model                                |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Slope                          | Map file name | :math:`mm^{-1}`         | Local terrain slope. Rise over run                     |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Horiz\_Hydraulic\_Conductivity | Map file name | :math:`ms^{-1}`         | Effective soil hydraulic conductivity                  |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Vert\_Horz\_Anis\_ratio        | Map file name | :math:`[-]`             | Ratio of vertical to horizontal hydraulic conductivity |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Terrain\_Random\_Roughness     | Map file name | :math:`m`               | Local surface roughness                                |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Porosity                       | Map file name | :math:`[-]`             | Soil porosity                                          |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Air\_entry\_pressure           | Map file name | :math:`m`               | Soil air entry pressure                                |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Brooks\_Corey\_lambda          | Map file name | :math:`[-]`             | Pore size distribution                                 |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Residual\_soil\_moisture       | Map file name | :math:`m^{3}m^{-3}`     | Minimum allowed volumetric soil water content          |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Soil\_depth                    | Map file name | :math:`m`               | Soil depth                                             |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Depth\_soil\_layer\_1          | Map file name | :math:`m`               | Depth of topmost soil layer                            |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Depth\_soil\_layer\_2          | Map file name | :math:`m`               | Depth of second soil layer                             |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Veget\_water\_use\_param1      | Map file name | :math:`m`               | | Vegetation water use parameter as per Landsberg and  |
|                                |               |                         | | Waring (1997)                                        |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Veget\_water\_use\_param2      | Map file name | :math:`m`               | | Vegetation water use parameter as per Landsberg and  |
|                                |               |                         | | Waring (1997)                                        |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Root\_profile\_coeff           | Map file name | :math:`m^{-1}`          | | Coefficient for the exponentiall-decreasing root     |
|                                |               |                         | | profile                                              |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Albedo                         | Map file name | :math:`[-]`             | Surface albedo                                         |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Surface\_emissivity            | Map file name | :math:`[-]`             | Surface emissivity/absorptivity                        |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Dry\_Soil\_Heat\_Capacity      | Map file name | :math:`Jm^{-3}K^{-1}`   | Heat capacity of soil solid particles                  |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Dry\_Soil\_Therm\_Cond         | Map file name | :math:`Wm^{-1}K^{-1}`   | Thermal conductivity of soil solid particles           |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Damping\_depth                 | Map file name | :math:`m`               | Depth of bottom of second soil thermal layer           |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Temp\_at\_damp\_depth          | Map file name | :math:`^{\circ}C`       | Soil temperature at damping depth                      |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Snow\_Melt\_Coeff              | Map file name | :math:`m^{\circ}C^{-1}` | Snowmelt coefficient factor                            |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
| Soil\_bedrock\_leakance        | Map file name | :math:`[-]`             | | Factor between 0 and 1 defining the vertical         |
|                                |               |                         | | hydraulic conductivity at the soil-bedrock interface |
|                                |               |                         | | (in proportion of soil Kv)                           |
+--------------------------------+---------------+-------------------------+--------------------------------------------------------+
  

**Forest parameters** (files must be located in ``Maps_Folder``)

+---------------------+-----------------+-------------+---------------------------------------------------------------------------+
| Keyword             | Type            | Unit        | Description                                                               |
+=====================+=================+=============+===========================================================================+
| ForestPatches       | Map file name   | integers    | Map identifying forest categories (patches)                               |
+---------------------+-----------------+-------------+---------------------------------------------------------------------------+
| Number\_of\_Species | Integer         | :math:`[-]` | Number of vegetation types included in the simulation                     |
+---------------------+-----------------+-------------+---------------------------------------------------------------------------+
| Species\_Parameters | Parameter table | :math:`[-]` | Table containing parameter information for each simulated vegetation type |
+---------------------+-----------------+-------------+---------------------------------------------------------------------------+


**Vegetation tables** (needed only if ``Species_State_Variable_Input_Method=tables``)

+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Keyword                     | Type           | Unit                 | Description                                                 |
+=============================+================+======================+=============================================================+
| Species\_Proportion\_Table  | Variable table | :math:`m^{2} m^{-2}` | | Table with initial proportion of covered area             |
|                             |                |                      | | (canopy cover) for each vegetation type with respect to   |
|                             |                |                      | | cell area                                                 |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Species\_StemDensity\_Table | Variable table | :math:`trees.m^{-2}` | Table with initial tree density for each vegetation type    |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Species\_LAI\_Table         | Variable table | :math:`m^{2} m^{-2}` | Table with initial leaf area index for each vegetation type |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+ 
| Species\_AGE\_Table         | Variable table | :math:`years`        | Table with initial average age each vegetation type         |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Species\_BasalArea\_Table   | Variable table | :math:`m^{2}`        | Table with initial total basal area per vegetation type     |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Species\_Height\_table      | Variable table | :math:`m`            | Table with initial effective height per vegetation type     |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+
| Species\_RootMass\_table    | Variable table | :math:`g m^{-3}`     | | Table with initial root mass per volume of soil for each  |
|                             |                |                      | | vegetation type                                           |
+-----------------------------+----------------+----------------------+-------------------------------------------------------------+

**Map report switches**

+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Keyword                               | Unit                   | Description                                         | File root |
+=======================================+========================+=====================================================+===========+
| Report\_Long\_Rad\_Down               | :math:`W m^{-2}`       | | Downwelling long wave (infrared) radiation at the | LDown     |
|                                       |                        | | top of the canopy (climate input)                 |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Short\_Rad\_Down              | :math:`W m^{-2}`       | | Incoming shortwave (visible) radiation at the top | Sdown     |
|                                       |                        | | of canopy (climate input)                         |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Precip                        | :math:`m s^{-1}`       | Precipitation (climate input)                       | Pp        |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Rel\_Humidity                 | :math:`Pa^{1} Pa^{-1}` | Relative humidity in the atmosphere (climate input) | RH        |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Wind\_Speed                   | :math:`m s^{-1}`       | Horizontal wind speed (climate input)               | WndSp     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_AvgAir\_Temperature           | :math:`^{\circ}C`      | Average air temperature (climate input)             | Tp        |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_MinAir\_Temperature           | :math:`^{\circ}C`      | Minimum air temperature (climate input)             | TpMin     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_MaxAir\_Temperature           | :math:`^{\circ}C`      | Maximum air temperature (climate input)             | TpMax     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_SWE                           | :math:`m`              | Snow water equivalent                               | SWE       |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Infilt\_Cap                   | :math:`m s^{-1}`       | Infiltration _Capacity                              | IfCap     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Streamflow                    | :math:`m^{3}s^{-1}`    | Channel discharge                                   | Q         |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Water\_Content\_Average | :math:`m^{3}m^{-3}`    | | Average volumetric water content for entire soil  | SWCav     |
|                                       |                        | | profile                                           |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Water\_Content\_Up      | :math:`m^{3}m^{-3}`    | | Average volumetric water content for the two      | SWCup     |
|                                       |                        | | upper soil layers                                 |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Water\_Content\_L1      | :math:`m^{3}m^{-3}`    | | Volumetric water content for topmost              | SWC1      |
|                                       |                        | | soil layer                                        |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Water\_Content\_L2      | :math:`m^{3}m^{-3}`    | | Volumetric water content for second               | SWC2      |
|                                       |                        | | soil layer                                        |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Water\_Content\_L3      | :math:`m^{3}m^{-3}`    | | Volumetric water content for bottommost           | SWC3      |
|                                       |                        | | soil layer                                        |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_WaterTableDepth               | :math:`m`              | | Depth the equivalent water table using the        | WTD       |
|                                       |                        | | average soil moisture                             |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Sat\_Deficit            | :math:`m`              | Meters of water needed to saturate soil             | SatDef    |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Ground\_Water                 | :math:`m`              | | Meters of water above field capacity in the third | GW        |
|                                       |                        | | hydrologic layer                                  |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Net\_Rad                | :math:`Wm^{-2}`        | Soil net radiation integrated over the grid cell    | NRs       |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_LE                      | :math:`Wm^{-2}`        | Latent heat for surface layer                       | LEs       |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Sens\_Heat                    | :math:`Wm^{-2}`        | Sensible heat for surface layer                     | SensH     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Grnd\_Heat                    | :math:`Wm^{-2}`        | Ground heat                                         | GrndH     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Snow\_Heat                    | :math:`Wm^{-2}`        | Turbulent heat exchange with snowpack               | SnowH     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Soil\_Temperature             | :math:`^{\circ}C`      | | Soil temperature at the bottom of first thermal   | Ts        |
|                                       |                        | | layer                                             |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Skin\_Temperature             | :math:`^{\circ}C`      | Soil skin temperature                               | Tskin     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Total\_ET                     | :math:`m s^{-1}`       | Total evapotranspiration                            | Evap      |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Transpiration\_sum            | :math:`m s^{-1}`       | | Transpiration integrated over the grid cell using | EvapT     |
|                                       |                        | | vegetation fractions                              |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Einterception\_sum            | :math:`m s^{-1}`       | | Evaporation of intercepted water integrated over  | EvapI     |
|                                       |                        | | the grid cell using vegetation fractions          |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Esoil\_sum                    | :math:`m s^{-1}`       | | Soil evaporation integrated over subcanopy and    | EvapS     |
|                                       |                        | | bare soil fractions                               |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Net\_Rad\_sum                 | :math:`Wm^{-2}`        | | Top-of-canopy net radiation integrated over the   | NRtot     |
|                                       |                        | | grid cell (including bare soil fraction)          |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Veget\_frac                   | :math:`m^{2} m^{-2}`   | | Fraction of cell covered by canopy of vegetation  | p\_*n*    |
|                                       |                        | | type *n*                                          |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Stem\_Density                 | :math:`stems m^{-2}`   | Density of individuals of vegetation type *n*       | ntr\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Leaf\_Area\_Index             | :math:`m^{2} m^{-2}`   | Leaf area index of vegetation type *n*              | lai\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Stand\_Age                    | :math:`years`          | Age of stand of vegetation type *n*                 | age\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_Conductance           | :math:`m s^{-1}`       | Canopy conductance for vegetation type *n*          | gc\_*n*   |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_GPP                           | :math:`gC m^{-2}`      | | Gross primary production for vegetation type *n*  | gpp\_*n*  |
|                                       |                        | | during the time step                              |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_NPP                           | :math:`gC^{-1} m^{-2}` | | Net primary production for vegetation type *n*    | npp\_*n*  |
|                                       |                        | | during the time step                              |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Basal\_Area                   | :math:`m^{2}`          | Total basal area of vegetation type *n*             | bas\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Tree\_Height                  | :math:`m`              | Height of stand of vegetation type *n*              | hgt\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Root\_Mass                    | :math:`g m^{-3}`       | Root mass per volume of soil vegetation type *n*    | root\_*n* |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_Temp                  | :math:`^{\circ}C`      | Canopy temperature of vegetation type *n*           | Tc\_*n*   |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_NetR                  | :math:`W m^{-2}`       | Net radiation above the vegetation type *n*         | NRc\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_LE\_E                 | :math:`W m^{-2}`       | | Latent heat for evaporation of canopy             | LEEi\_*n* |
|                                       |                        | | interception for vegetation type *n*              |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_LE\_T                 | :math:`W m^{-2}`       | Transpiration latent heat for vegetation type *n*   | LETr\_*n* |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_Sens\_Heat            | :math:`W m^{-2}`       | Sensible heat, canopy layer of vegetation type *n*  | Hc\_*n*   |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Canopy\_Water\_Stor           | :math:`m`              | Intercepted water storage of vegetation type *n*    | Cs\_*n*   |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_species\_ET                   | :math:`m s^{-1}`       | Evapotranspiration within the vegetation type *n*   | ETc\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Transpiration                 | :math:`m s^{-1}`       | Transpiration from vegetation type *n*              | Trp\_*n*  |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Einterception                 | :math:`m s^{-1}`       | | Evaporation of intercepted water for the          | Ei\_*n*   |
|                                       |                        | | vegetation type *n*                               |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Esoil                         | :math:`m s^{-1}`       | Soil evaporation under the vegetation type *n*      | Es\_*n*   |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_GW\_to\_Channnel              | :math:`m`              |  Quantity of groundwater seeping in stream water    | GWChn     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Surface\_to\_Channel          | :math:`m`              | | Quantity of surface runoff contributing to        | SrfChn    |
|                                       |                        | | stream water                                      |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Infiltration                  | :math:`m`              | | Meters of water (re)infiltrated water in the      | Inf       |
|                                       |                        | | first hydrological layer                          |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Return\_Flow\_Surface         | :math:`m`              | | Meters of water exfiltrated from the first        | RSrf      |
|                                       |                        | | hydrological layer                                |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Overland\_Inflow              | :math:`m`              | Surface run-on (excluding channel inflow)           | LSrfi     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Stream\_Inflow                | :math:`m`              | Incoming stream water                               | LChni     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Groundwater\_Inflow           | :math:`m`              | Lateral groundwater inflow                          | LGWi      |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Overland\_Outflow             | :math:`m`              | Surface run-off (excluding channel outflow)         | LSrfo     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Groundwater\_Outflow          | :math:`m`              | Lateral groundwater outflow                         | LGWo      |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_GW\_to\_Channnel\_acc         | :math:`m`              | | Cumulated quantity of groundwater seeping in      | GWChnA    |
|                                       |                        | | stream water                                      |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Surface\_to\_Channel\_acc     | :math:`m`              | | Cumulated quantity of surface runoff contributing | SrfChnA   |
|                                       |                        | |  to stream water                                  |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Infiltration\_acc             | :math:`m`              | | Cumulated meters of water (re)infiltrated water   | InfA      |
|                                       |                        | | in the first hydrological layer                   |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Return\_Flow\_Surface\_acc    | :math:`m`              | | Cumulated meters of water exfiltrated from the    | RSrfA     |
|                                       |                        | | first hydrological layer                          |           |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Overland\_Inflow\_acc         | :math:`m`              | Cumulated surface run-on (excluding channel inflow) | LSrfiA    |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Stream\_Inflow\_acc           | :math:`m`              | Cumulated lncoming stream water                     | LChniA    |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Groundwater\_Inflow\_acc      | :math:`m`              | Cumulated lateral groundwater inflow                | LGWiA     |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Overland\_Outflow\_acc        | :math:`m`              | Cumulated surface run-off (excludes discharge)      | LSrfoA    |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+
| Report\_Groundwater\_Outflow\_acc     | :math:`m`              | Cumulated lateral groundwater outflow               | LGWo      |
+---------------------------------------+------------------------+-----------------------------------------------------+-----------+


**Map mask for time series locations**

+----------+---------------+---------------------------------------------------------------------+
| Keyword  | Type          | Description                                                         |
+==========+===============+=====================================================================+
| TS\_mask | Map file name | | Map identifying cells for which state variables will be reported. |
|          |               | | Map should be zero everywhere except for target cells.            |
|          |               | | A maximum of 32 cells can be reported.                            |
+----------+---------------+---------------------------------------------------------------------+

**Time series report switches**

Written outputs file are time series tables at cells identified in ``TS_mask``.

+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Keyword                           | Unit                   | Description                                         | File name                |
+===================================+========================+=====================================================+==========================+
| Ts\_OutletDischarge               | :math:`m^{3} s^{-1}`   | | Discharge at cells with *ldd* value = 5 (outlets  | OutletDisch.tab          |
|                                   |                        | |  and sinks)                                       |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Long\_Rad\_Down               | :math:`W m^{-2}`       | | Downwelling long wave (infrared) radiation at the | LDown.tab                |
|                                   |                        | | top of the canopy (climate input)                 |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Short\_Rad\_Down              | :math:`W m^{-2}`       | | Incoming shortwave (visible) radiation at the top | Sdown.tab                |
|                                   |                        | | of canopy (climate input)                         |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Precip                        | :math:`m s^{-1}`       | Precipitation (climate input)                       | Precip.tab               |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Rel\_Humidity                 | :math:`Pa^{1} Pa^{-1}` | Relative humidity in the atmosphere (climate input) | RelHumid.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Wind\_Speed                   | :math:`m s^{-1}`       | Horizontal wind speed (climate input)               | WindSpeed.tab            |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_AvgAir\_Temperature           | :math:`^{\circ}C`      | Average air temperature (climate input)             | AvgTemp.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_MinAir\_Temperature           | :math:`^{\circ}C`      | Minimum air temperature (climate input)             | MinTemp.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_MaxAir\_Temperature           | :math:`^{\circ}C`      | Maximum air temperature (climate input)             | MaxTemp.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_SWE                           | :math:`m`              | Snow water equivalent                               | SWE.tab                  |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Infilt\_Cap                   | :math:`m s^{-1}`       | Infiltration _Capacity                              | InfiltCap.tab            |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Streamflow                    | :math:`m^{3}s^{-1}`    | Channel discharge                                   | Streamflow.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Ponding                       | :math:`m^{3}m^{-3}`    | | Surface water height                              | Ponding.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Water\_Content\_Average | :math:`m^{3}m^{-3}`    | | Average volumetric water content for entire soil  | SoilMoistureAv.tab       |
|                                   |                        | | profile                                           |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Water\_Content\_Up      | :math:`m^{3}m^{-3}`    | | Average volumetric water content for the two      | SoilMoistureUp.tab       |
|                                   |                        | | upper soil layers                                 |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Water\_Content\_L1      | :math:`m^{3}m^{-3}`    | | Volumetric water content for topmost              | SoilMoistureL1.tab       |
|                                   |                        | | soil layer                                        |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Water\_Content\_L2      | :math:`m^{3}m^{-3}`    | | Volumetric water content for second               | SoilMoistureL2.tab       |
|                                   |                        | | soil layer                                        |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Water\_Content\_L3      | :math:`m^{3}m^{-3}`    | | Volumetric water content for bottommost           | SoilMoistureL3.tab       |
|                                   |                        | | soil layer                                        |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_WaterTableDepth               | :math:`m`              | | Depth the equivalent water table using the        | WaterTableDepth.tab      |
|                                   |                        | | average soil moisture                             |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Sat\_Deficit            | :math:`m`              | | Water depth needed to saturate the cells          | SoilSatDef.tab           |
|                                   |                        | | identified in *TS\_mask*                          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Ground\_Water                 | :math:`m`              | | Meters of water above field capacity in the third | GroundWater.tab          |
|                                   |                        | | hydrologic layer                                  |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Net\_Rad                | :math:`Wm^{-2}`        | Soil net radiation integrated over the grid cell    | NetRadS.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_LE                      | :math:`Wm^{-2}`        | Latent heat for surface layer                       | LatHeat.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Sens\_Heat                    | :math:`Wm^{-2}`        | Sensible heat for surface layer                     | SensHeat.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Grnd\_Heat                    | :math:`Wm^{-2}`        | Ground heat                                         | GrndHeat.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Snow\_Heat                    | :math:`Wm^{-2}`        | Turbulent heat exchange with snowpack               | SnowHeat.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Soil\_Temperature             | :math:`^{\circ}C`      | | Soil temperature at the bottom of first thermal   | SoilTemp.tab             |
|                                   |                        | | layer                                             |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Skin\_Temperature             | :math:`^{\circ}C`      | Soil skin temperature                               | SkinTemp.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Total\_ET                     | :math:`m s^{-1}`       | Total evapotranspiration                            | Evap.tab                 |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Transpiration\_sum            | :math:`m s^{-1}`       | | Transpiration integrated over the grid cell using | EvapT.tab                |
|                                   |                        | | vegetation fractions                              |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Einterception\_sum            | :math:`m s^{-1}`       | | Evaporation of intercepted water integrated over  | EvapI.tab                |
|                                   |                        | | the grid cell using vegetation fractions          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Esoil\_sum                    | :math:`m s^{-1}`       | | Soil evaporation integrated over subcanopy and    | EvapS.tab                |
|                                   |                        | | bare soil fractions                               |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Net\_Rad\_sum                 | :math:`Wm^{-2}`        | | Top-of-canopy net radiation integrated over the   | NetRadtot.tab            |
|                                   |                        | | grid cell (including bare soil fraction)          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Veget\_frac                   | :math:`m^{2} m^{-2}`   | | Fraction of cell covered by canopy of vegetation  | p\_*n*.tab               |
|                                   |                        | | type *n*                                          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Stem\_Density                 | :math:`stems m^{-2}`   | Density of individuals of vegetation type *n*       | num\_of\_trees\_*n*.tab  |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Leaf\_Area\_Index             | :math:`m^{2} m^{-2}`   | Leaf area index of vegetation type *n*              | lai\_*n*.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_Conductance           | :math:`m s^{-1}`       | Canopy conductance for vegetation type *n*          | CanopyConduct\_*n*.tab   |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_GPP                           | :math:`gC m^{-2}`      | | Gross primary production for vegetation type *n*  | GPP\_*n*.tab             |
|                                   |                        | | during the time step                              |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_NPP                           | :math:`gC^{-1} m^{-2}` | | Net primary production for vegetation type *n*    | NPP\_*n*.tab             |
|                                   |                        | | during the time step                              |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Basal\_Area                   | :math:`m^{2}`          | Total basal area of vegetation type *n*             | BasalArea\_*n*.tab       |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Tree\_Height                  | :math:`m`              | Height of stand of vegetation type *n*              | TreeHeight\_*n*.tab      |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Root\_Mass                    | :math:`g m^{-3}`       | Root mass per volume of soil vegetation type *n*    | RootMass\_*n*.tab        |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_Temp                  | :math:`^{\circ}C`      | Canopy temperature of vegetation type *n*           | CanopyTemp\_*n*.tab      |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_NetR                  | :math:`W m^{-2}`       | Net radiation above the vegetation type *n*         | NetRadC\_*n*.tab         |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_LE\_E                 | :math:`W m^{-2}`       | | Latent heat for evaporation of canopy             | CanopyLatHeatEi\_*n*.tab |
|                                   |                        | | interception for vegetation type *n*              |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_LE\_T                 | :math:`W m^{-2}`       | Transpiration latent heat for vegetation type *n*   | CanopyLatHeatTr\_*n*.tab |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_Sens\_Heat            | :math:`W m^{-2}`       | Sensible heat, canopy layer of vegetation type *n*  | CanopySensHeat\_*n*.tab  |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Canopy\_Water\_Stor           | :math:`m`              | Intercepted water storage of vegetation type *n*    | CanopyWaterStor\_*n*.tab |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_species\_ET                   | :math:`m s^{-1}`       | Evapotranspiration within the vegetation type *n*   | ETc\_*n*.tab             |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Transpiration                 | :math:`m s^{-1}`       | Transpiration from vegetation type *n*              | EvapT\_*n*.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Einterception                 | :math:`m s^{-1}`       | | Evaporation of intercepted water for the          | EvapI\_*n*.tab           |
|                                   |                        | | vegetation type *n*                               |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Esoil                         | :math:`m s^{-1}`       | Soil evaporation under the vegetation type *n*      | EvapS\_*n*.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_GW\_to\_Channnel              | :math:`m`              |  Quantity of groundwater seeping in stream water    | GWtoChn.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Surface\_to\_Channel          | :math:`m`              | | Quantity of surface runoff contributing to        | SrftoChn.tab             |
|                                   |                        | | stream water                                      |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Infiltration                  | :math:`m`              | | Meters of water (re)infiltrated water in the      | Infilt.tab               |
|                                   |                        | | first hydrological layer                          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Return\_Flow\_Surface         | :math:`m`              | | Meters of water exfiltrated from the first        | ReturnSrf.tab            |
|                                   |                        | | hydrological layer                                |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Overland\_Inflow              | :math:`m`              | Surface run-on (excluding channel inflow)           | SrfLatI.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Stream\_Inflow                | :math:`m`              | Incoming stream water                               | ChnLatI.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Groundwater\_Inflow           | :math:`m`              | Lateral groundwater inflow                          | GWLatI.tab               |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Overland\_Outflow             | :math:`m`              | Surface run-off (excluding channel outflow)         | SrfLatO.tab              |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Groundwater\_Outflow          | :math:`m`              | Lateral groundwater outflow                         | GWLatO.tab               |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_GW\_to\_Channnel\_acc         | :math:`m`              | | Cumulated quantity of groundwater seeping in      | GWtoChnAcc.tab           |
|                                   |                        | | stream water                                      |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Surface\_to\_Channel\_acc     | :math:`m`              | | Cumulated quantity of surface runoff contributing | SrftoChnAcc.tab          |
|                                   |                        | |  to stream water                                  |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Infiltration\_acc             | :math:`m`              | | Cumulated meters of water (re)infiltrated water   | InfiltAcc.tab            |
|                                   |                        | | in the first hydrological layer                   |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Return\_Flow\_Surface\_acc    | :math:`m`              | | Cumulated meters of water exfiltrated from the    | ReturnSrfAcc.tab         |
|                                   |                        | | first hydrological layer                          |                          |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Overland\_Inflow\_acc         | :math:`m`              | Cumulated surface run-on (excluding channel inflow) | SrfLatIAcc.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Stream\_Inflow\_acc           | :math:`m`              | Cumulated lncoming stream water                     | ChnLatIAcc.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Groundwater\_Inflow\_acc      | :math:`m`              | Cumulated lateral groundwater inflow                | GWLatIAcc.tab            |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Overland\_Outflow\_acc        | :math:`m`              | Cumulated surface run-off (excludes discharge)      | SrfLatOAcc.tab           |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
| Ts\_Groundwater\_Outflow\_acc     | :math:`m`              | Cumulated lateral groundwater outflow               | GWLatOAcc.tab            |
+-----------------------------------+------------------------+-----------------------------------------------------+--------------------------+
