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
 * CheckMaps.cpp
 *
 *  Created on: Feb 10, 2011
 *      Author: Marco.Maneta
 */

#include "Tracking.h"

void Tracking::CheckMapsTrck(Control &ctrl, Basin &bsn) {

  UINT4 r, c;
  UINT4 j = 0;
  bool excep_thrown = false; //  poor man way to rethrow exception outside omp pragma
  UINT4 length = bsn.getSortedGrid().cells.size();
#pragma omp parallel for			\
  default(shared) private(r,c,j)
  for (j = 0; j < length; j++) {
    r = bsn.getSortedGrid().cells[j].row;
    c = bsn.getSortedGrid().cells[j].col;
    try {

      if(ctrl.sw_trck && ctrl.sw_2H){
	if (getd2Hsnowpack()->matrix[r][c] == getd2Hsnowpack()->nodata) {
	  string e("Initial d2H map in snowpack contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd2Hsurface()->matrix[r][c] == getd2Hsurface()->nodata) {
	  string e("Initial d2H map in channel contains no data values inside the valid channel domain...\n");
	  throw e;}
	if (getd2Hsoil1()->matrix[r][c] ==getd2Hsoil1()->nodata) {
	  string e("Initial d2H map in layer 1 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd2Hsoil2()->matrix[r][c] ==getd2Hsoil2()->nodata) {
	  string e("Initial d2H map in layer 2 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd2Hsoil3()->matrix[r][c] == getd2Hsoil3()->nodata) {
	  string e("Initial d2H map in layer 3 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd2Hgroundwater()->matrix[r][c] == getd2Hgroundwater()->nodata) {
	  string e("Initial d2H map in groundwater contains no data values inside the valid domain...\n");
	  throw e;}
      }

      if(ctrl.sw_trck && ctrl.sw_18O){
	if (getd18Osnowpack()->matrix[r][c] == getd18Osnowpack()->nodata) {
	  string e("Initial d18O map in snowpack contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd18Osurface()->matrix[r][c] == getd18Osurface()->nodata) {
	  string e("Initial d18O map in channel contains no data values inside the valid channel domain...\n");
	  throw e;}
	if (getd18Osoil1()->matrix[r][c] == getd18Osoil1()->nodata) {
	  string e("Initial d18O map in layer 1 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd18Osoil2()->matrix[r][c] == getd18Osoil2()->nodata) {
	  string e("Initial d18O map in layer 2 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd18Osoil3()->matrix[r][c] == getd18Osoil3()->nodata) {
	  string e("Initial d18O map in layer 3 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getd18Ogroundwater()->matrix[r][c] == getd18Ogroundwater()->nodata) {
	  string e("Initial d18O map in groundwater contains no data values inside the valid domain...\n");
	  throw e;}
      }

      if(ctrl.sw_trck && ctrl.sw_Age){
	if (getAgesnowpack()->matrix[r][c] == getAgesnowpack()->nodata) {
	  string e("Initial Age map in snowpack contains no data values inside the valid domain...\n");
	  throw e;}
	if (getAgesurface()->matrix[r][c] == getAgesurface()->nodata) {
	  string e("Initial Age map in channel contains no data values inside the valid channel domain...\n");
	  throw e;}
	if (getAgesoil1()->matrix[r][c] == getAgesoil1()->nodata) {
	  string e("Initial Age map in layer 1 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getAgesoil2()->matrix[r][c] == getAgesoil2()->nodata) {
	  string e("Initial Age map in layer 2 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getAgesoil3()->matrix[r][c] == getAgesoil3()->nodata) {
	  string e("Initial Age map in layer 3 contains no data values inside the valid domain...\n");
	  throw e;}
	if (getAgegroundwater()->matrix[r][c] == getAgegroundwater()->nodata) {
	  string e("Initial Age map in groundwater contains no data values inside the valid domain...\n");
	  throw e;}
      }

    } catch (string &e) {
      cout << e;
      cout << "In row " << r << " col " << c << endl;
    }

  }

  if (excep_thrown)
    throw;
}
