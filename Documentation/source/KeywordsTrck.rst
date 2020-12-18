Configuration (tracking) file keywords
======================================

**Boolean switches**

+-----------------+--------------+--------------------------------------------------------------+
| Keyword         | Description                                                                 |
+=================+==============+==============================================================+
| water\_dD       | Switch 1/0 to turn deuterium tracking on (1) or off (0)                     |
+-----------------+-----------------------------------------------------------------------------+
| water\_d18O     | Switch 1/0 to turn oxygen-18 tracking on (1) or off (0)                     |
+-----------------+-----------------------------------------------------------------------------+
| water\_Age      | Switch 1/0 to turn water age tracking on (1) or off (0)                     |
+-----------------+-----------------------------------------------------------------------------+
| water\_frac     | Switch 1/0 to turn evaporative fractionation in the soil on (1) or off (0)  |
+-----------------+-----------------------------------------------------------------------------+
| water\_lifo     | | Switch 1/0. If activated, same-timestep topsoil infiltration if the       |
|                 | | primary pool for soil evaporation, only the remainde (if any)             |
|                 | | subsequently mixes with older topsoil water (LastInFirstOut)              |
+-----------------+-----------------------------------------------------------------------------+

**Toggle switches** (used only if ``water_frac`` = 1)

+----------------------------------+--------------------------------------------------------------+
| Keyword                          | Description                                                  |
+==================================+==============================================================+
| Fractionation\_surface\_relhum   | | Relative humidity in air space between soil pores.         |
|                                  | | 0: *soilRH* =1; 1: *soilRH* follows Lee and Pielke (1992); |
|                                  | | 2: *soilRH* follows Soderberg et al. (2012)                |
+----------------------------------+--------------------------------------------------------------+
| Fractionation\_turbulent\_factor | | Choice of turbulent factor *n* in kinetic fractionation.   |
|                                  | | 0: *n* =1; 1: *n* depends on soil water content, following |
|                                  | | Mathieu and Bariac (1996)                                  |
+----------------------------------+--------------------------------------------------------------+

**Climate information for isotopes** 

Binary files must be placed in ``Clim_Maps_Folder``.
Only necessary if the corresponding switches (``water_dD``, ``water_d18O``) are set to 1. 

+--------------+--------------------------+------+--------------------+
| Keyword      | Type                     | Unit | Description        |
+==============+==========================+======+====================+
| dD\_precip   | Binary climate file name | ‰    | Deuterium input    |
+--------------+--------------------------+------+--------------------+
| d18O\_precip | Binary climate file name | ‰    | Oxygen 18 input    |
+--------------+--------------------------+------+--------------------+

**Initial conditions for water tracers** 

Files must be located in ``Maps_Folder``.
Only necessary if the corresponding switches (``water_dD``, ``water_d18O``, ``water_Age``) are set to 1. 

+-------------------------+---------------+--------------+----------------------------------------------------------+
| Keyword                 | Type          | Unit         | Description                                              |
+=========================+===============+==============+==========================================================+
| init\_dD\_snowpack      | Map file name | ‰            | Deuterium ratio in the snowpack                          |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_dD\_surface       | Map file name | ‰            | Deuterium ratio in the surface water (including channel) |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_dD\_soil1         | Map file name | ‰            | Deuterium ratio in the topsoil water                     |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_dD\_soil2         | Map file name | ‰            | Deuterium ratio in the water of layer 2                  |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_dD\_soil3         | Map file name | ‰            | | Deuterium ratio of the water beneath field capacity    |
|                         |               |              | | in bottommost hydrological layer                       |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_dD\_groundwater   | Map file name | ‰            | | Deuterium ratio in groundwater                         |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_snowpack    | Map file name | ‰            | Oxygen 18 ratio in the snowpack                          |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_surface     | Map file name | ‰            | Oxygen 18 ratio in the surface water (including channel) |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_soil1       | Map file name | ‰            | Oxygen 18 ratio in the topsoil water                     |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_soil2       | Map file name | ‰            | Oxygen 18 ratio in the water of layer 2                  |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_soil3       | Map file name | ‰            | | Oxygen 18 ratio of the water beneath field capacity    |
|                         |               |              | | in bottommost hydrological layer                       |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_d18O\_groundwater | Map file name | ‰            | | Oxygen 18 ratio in groundwater                         |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_snowpack     | Map file name | :math:`days` | Age of the snowpack                                      |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_surface      | Map file name | :math:`days` | Age of the surface water (including channel)             |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_soil1        | Map file name | :math:`days` | Age of the topsoil water                                 |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_soil2        | Map file name | :math:`days` | Age of the water in layer 2                              |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_soil3        | Map file name | :math:`days` | | Age of the water beneath field capacity in the         |
|                         |               |              | | in bottommost hydrological layer                       |
+-------------------------+---------------+--------------+----------------------------------------------------------+
| init\_Age\_groundwater  | Map file name | :math:`days` | | Groundwater age                                        |
+-------------------------+---------------+--------------+----------------------------------------------------------+


**Map report switches**

+----------------------+--------------+----------------------------------------------------------+--------------+
| Keyword              | Unit         | Description                                              | File root    |
+======================+==============+==========================================================+==============+
| Rep\_dDprecip        | ‰            | Deuterium ratio in precipitation (climate input)         | dDpcp        |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsnowpack      | ‰            | Deuterium ratio in the snowpack                          | dDsnw        |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsurface       | ‰            | Deuterium ratio in the surface water (including channel) | dDsrf        |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsoil1         | ‰            | Deuterium ratio in the topsoil water                     | dDsL1        |
+----------------------+--------------+----------------------------------------------------------+--------------+ 
| Rep\_dDsoil2         | ‰            | Deuterium ratio in the water of layer 2                  | dDsL2        |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsoilUp        | ‰            | | Deuterium ratio averaged over water content in the two | dDsUp        |
|                      |              | | upper hydrological layers                              |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsoil3         | ‰            | | Deuterium ratio of the water beneath field capacity    | dDsL3        |
|                      |              | | in bottommost hydrological layer                       |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDsoilAv        | ‰            | | Deuterium ratio averaged over water content in the     | dDsAv        |
|                      |              | | soil profile (groundwater included)                    |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDgroundwater   | ‰            | | Deuterium ratio in groundwater                         | dDgw         |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapS         | ‰            | | Deuterium ratio in water evaporated from topsoil under | dDeS\_ *n*   |
|                      |              | | vegetation type *n*                                    |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapS\_sum    | ‰            | | Deuterium ratio in water evaporated from topsoil,      | dDeS         |
|                      |              | | summed over the sub-canopy and bare soil fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapI         | ‰            | | Deuterium ratio in water evaporated from canopy        | dDeI\_ *n*   |
|                      |              | | interception for vegetation type *n*                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapI\_sum    | ‰            | | Deuterium ratio in water evaporated from canopy        | dDeI         |
|                      |              | | interception, summed over the vegetation fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapT         | ‰            | | Deuterium ratio in root water uptake for vegetation    | dDeT\_ *n*   |
|                      |              | | type *n*                                               |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_dDevapT\_sum    | ‰            | | Deuterium ratio in root water uptake, summed over      | dDeT         |
|                      |              | | vegetation fractions                                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Oprecip      | ‰            | Oxygen 18 ratio in precipitation (climate input)         | d18Opcp      |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Osnowpack    | ‰            | Oxygen 18 ratio in the snowpack                          | d18Osnw      |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Osurface     | ‰            | Oxygen 18 ratio in the surface water (including channel) | d18Osrf      |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Osoil1       | ‰            | Oxygen 18 ratio in the topsoil water                     | d18OsL1      |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Osoil2       | ‰            | Oxygen 18 ratio in the water of layer 2                  | d18OsL2      |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OsoilUp      | ‰            | | Oxygen 18 ratio averaged over water content in the two | d18OsUp      |
|                      |              | | upper hydrological layers                              |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Osoil3       | ‰            | | Oxygen 18 ratio of the water beneath field capacity    | d18OsL3      |
|                      |              | | in bottommost hydrological layer                       |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OsoilAv      | ‰            | | Oxygen 18 ratio averaged over water content in the     | d18OsAv      |
|                      |              | | soil profile (groundwater included)                    |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18Ogroundwater | ‰            | | Oxygen 18 ratio in groundwater                         | d18Ogw       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapS       | ‰            | | Oxygen 18 ratio in water evaporated from topsoil under | d18OeS\_ *n* |
|                      |              | | vegetation type *n*                                    |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapS\_sum  | ‰            | | Oxygen 18 ratio in water evaporated from topsoil,      | d18OeS       |
|                      |              | | summed over the sub-canopy and bare soil fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapI       | ‰            | | Oxygen 18 ratio in water evaporated from canopy        | d18OeI\_ *n* |
|                      |              | | interception for vegetation type *n*                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapI\_sum  | ‰            | | Oxygen 18 ratio in water evaporated from canopy        | d18OeI       |
|                      |              | | interception, summed over the vegetation fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapT       | ‰            | | Oxygen 18 ratio in root water uptake for vegetation    | d18OeT\_ *n* |
|                      |              | | type *n*                                               |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_d18OevapT\_sum  | ‰            | | Oxygen 18 ratio in root water uptake, summed over      | d18OeT       |
|                      |              | | vegetation fractions                                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agesnowpack     | :math:`days` | Age of the snowpack                                      | Agesnw       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agesurface      | :math:`days` | Age of the surface water (including channel)             | Agesrf       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agesoil1        | :math:`days` | Age of the topsoil water                                 | AgesL1       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agesoil2        | :math:`days` | Age of the water in layer 2                              | AgesL2       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgesoilUp       | :math:`days` | Average water age in the two upper hydrological layers   | AgesUp       |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agesoil3        | :math:`days` | | Age of the water beneath field capacity in the         | AgesL3       |
|                      |              | | in bottommost hydrological layer                       |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgesoilAv       | :math:`days` | | Average water age over the soil profile                | AgesAv       |
|                      |              | | (groundwater included)                                 |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_Agegroundwater  | :math:`days` | | Groundwater age                                        | Agegw        |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapS        | :math:`days` | | Age of water evaporated from topsoil under             | AgeeS\_ *n*  |
|                      |              | | vegetation type *n*                                    |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapS\_sum   | :math:`days` | | Age of in water evaporated from topsoil,               | AgeeS        |
|                      |              | | summed over the sub-canopy and bare soil fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapI        | :math:`days` | | Age of in water evaporated from canopy                 | AgeeI\_ *n*  |
|                      |              | | interception for vegetation type *n*                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapI\_sum   | :math:`days` | | Age of in water evaporated from canopy                 | AgeeI        |
|                      |              | | interception, summed over the vegetation fractions     |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapT        | :math:`days` | | Age of in root water uptake for vegetation             | AgeeT\_ *n*  |
|                      |              | | type *n*                                               |              |
+----------------------+--------------+----------------------------------------------------------+--------------+
| Rep\_AgeevapT\_sum   | :math:`days` | | Age of in root water uptake, summed over               | AgeeT        |
|                      |              | | vegetation fractions                                   |              |
+----------------------+--------------+----------------------------------------------------------+--------------+


**Time series report switches**

Written outputs file are time series tables at cells identified in ``TS_mask`` (see main configuration file).

+---------------------+--------------+----------------------------------------------------------+----------------------+
| Keyword             | Unit         | Description                                              | File root            |
+=====================+==============+==========================================================+======================+
| Ts\_dDprecip        | ‰            | Deuterium ratio in precipitation (climate input)         | dD_precip.tab        |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsnowpack      | ‰            | Deuterium ratio in the snowpack                          | dD_snowpack.tab      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsurface       | ‰            | Deuterium ratio in the surface water (including channel) | dD_surface.tab       |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsoil1         | ‰            | Deuterium ratio in the topsoil water                     | dD_soilL1.tab        |
+---------------------+--------------+----------------------------------------------------------+----------------------+ 
| Ts\_dDsoil2         | ‰            | Deuterium ratio in the water of layer 2                  | dD_soilL2.tab        |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsoilUp        | ‰            | | Deuterium ratio averaged over water content in the two | dD_soilUp.tab        |
|                     |              | | upper hydrological layers                              |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsoil3         | ‰            | | Deuterium ratio of the water beneath field capacity    | dD_soilL3.tab        |
|                     |              | | in bottommost hydrological layer                       |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDsoilAv        | ‰            | | Deuterium ratio averaged over water content in the     | dD_soilAv.tab        |
|                     |              | | soil profile (groundwater included)                    |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDgroundwater   | ‰            | Deuterium ratio in groundwater                           | dD_groundwater.tab   |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapS         | ‰            | | Deuterium ratio in water evaporated from topsoil under | dDevapS\_ *n*.tab    |
|                     |              | | vegetation type *n*                                    |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapS\_sum    | ‰            | | Deuterium ratio in water evaporated from topsoil,      | dD_evapS.tab         |
|                     |              | | summed over the sub-canopy and bare soil fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapI         | ‰            | | Deuterium ratio in water evaporated from canopy        | dDevapI\_ *n*.tab    |
|                     |              | | interception for vegetation type *n*                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapI\_sum    | ‰            | | Deuterium ratio in water evaporated from canopy        | dD_evapI.tab         |
|                     |              | | interception, summed over the vegetation fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapT         | ‰            | | Deuterium ratio in root water uptake for vegetation    | dDevapT\_ *n*.tab    |
|                     |              | | type *n*                                               |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_dDevapT\_sum    | ‰            | | Deuterium ratio in root water uptake, summed over      | dD_evapT.tab         |
|                     |              | | vegetation fractions                                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Oprecip      | ‰            | Oxygen 18 ratio in precipitation (climate input)         | d18O_precip.tab      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Osnowpack    | ‰            | Oxygen 18 ratio in the snowpack                          | d18O_snowpack.tab    |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Osurface     | ‰            | Oxygen 18 ratio in the surface water (including channel) | d18O_surface.tab     |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Osoil1       | ‰            | Oxygen 18 ratio in the topsoil water                     | d18O_soilL1.tab      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Osoil2       | ‰            | Oxygen 18 ratio in the water of layer 2                  | d18O_soilL2.tab      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OsoilUp      | ‰            | | Oxygen 18 ratio averaged over water content in the two | d18O_soilUp.tab      |
|                     |              | | upper hydrological layers                              |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Osoil3       | ‰            | | Oxygen 18 ratio of the water beneath field capacity    | d18O_soilL3.tab      |
|                     |              | | in bottommost hydrological layer                       |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OsoilAv      | ‰            | | Oxygen 18 ratio averaged over water content in the     | d18O_soilAv.tab      |
|                     |              | | soil profile (groundwater included)                    |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18Ogroundwater | ‰            | Oxygen 18 ratio in groundwater                           | d18O_groundwater.tab |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapS       | ‰            | | Oxygen 18 ratio in water evaporated from topsoil under | d18OevapS\_ *n*.tab  |
|                     |              | | vegetation type *n*                                    |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapS\_sum  | ‰            | | Oxygen 18 ratio in water evaporated from topsoil,      | d18O_evapS.tab       |
|                     |              | | summed over the sub-canopy and bare soil fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapI       | ‰            | | Oxygen 18 ratio in water evaporated from canopy        | d18OevapI\_ *n*.tab  |
|                     |              | | interception for vegetation type *n*                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapI\_sum  | ‰            | | Oxygen 18 ratio in water evaporated from canopy        | d18O_evapI.tab       |
|                     |              | | interception, summed over the vegetation fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapT       | ‰            | | Oxygen 18 ratio in root water uptake for vegetation    | d18OevapT\_ *n*.tab  |
|                     |              | | type *n*                                               |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_d18OevapT\_sum  | ‰            | | Oxygen 18 ratio in root water uptake, summed over      | d18O_evapT.tab       |
|                     |              | | vegetation fractions                                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agesnowpack     | :math:`days` | Age of the snowpack                                      | Age_snowpack.tab     |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agesurface      | :math:`days` | Age of the surface water (including channel)             | Age_surface.tab      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agesoil1        | :math:`days` | Age of the topsoil water                                 | Age_soilL1.tab       |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agesoil2        | :math:`days` | Age of the water in layer 2                              | Age_soilL2.tab       |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgesoilUp       | :math:`days` | Average water age in the two upper hydrological layers   | Age_soilUp.tab       |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agesoil3        | :math:`days` | | Age of the water beneath field capacity in the         | Age_soilL3.tab       |
|                     |              | | in bottommost hydrological layer                       |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgesoilAv       | :math:`days` | | Average water age over the soil profile                | Age_soilAv.tab       |
|                     |              | | (groundwater included)                                 |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_Agegroundwater  | :math:`days` | Groundwater age                                          | Age_groundwater.tab  |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapS        | :math:`days` | | Age of water evaporated from topsoil under             | AgeevapS\_ *n*.tab   |
|                     |              | | vegetation type *n*                                    |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapS\_sum   | :math:`days` | | Age of water evaporated from topsoil,                  | Age_evapS.tab        |
|                     |              | | summed over the sub-canopy and bare soil fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapI        | :math:`days` | | Age of water evaporated from canopy                    | AgeevapI\_ *n*.tab   |
|                     |              | | interception for vegetation type *n*                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapI\_sum   | :math:`days` | | Age of water evaporated from canopy                    | Age_evapI.tab        |
|                     |              | | interception, summed over the vegetation fractions     |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapT        | :math:`days` | | Age of root water uptake for vegetation                | AgeevapT\_ *n*.tab   |
|                     |              | | type *n*                                               |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+
| Ts\_AgeevapT\_sum   | :math:`days` | | Age of root water uptake, summed over                  | Age_evapT.tab        |
|                     |              | | vegetation fractions                                   |                      |
+---------------------+--------------+----------------------------------------------------------+----------------------+


