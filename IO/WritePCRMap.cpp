/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta
 *******************************************************************************/
/*
 * WritePCRMap.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#include <fstream>
#include "Grid.h"

int grid::grid2PCRMap(std::string fname,
                        CSF_CR cellRepr,
                        CSF_VS dataType)
{
  MAP *out; REAL4 data = 0; UINT4 nrCells = 0; CSF_PT proj;

  if (north > south)
    proj = PT_YDECT2B;
  else
    proj = PT_YINCT2B;

  out = Rcreate(fname.c_str(), r, c, cellRepr, dataType, proj, west, north,0, dx);

  if (Merrno != 0)
    {cout << MstrError(); return -1;}
  RuseAs(out, CR_REAL4);

  for(unsigned i = 0; i < r; i++) {
    for(unsigned j = 0; j < c; j++) {
      data = (REAL4)matrix[i][j];
      if (data == nodata)
	SET_MV_REAL4(&data);
      nrCells += RputCell(out, i, j, &data);
    }
  }
  RuseAs(out, cellRepr);
  Mclose(out);
  return nrCells;
}
