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
 * ReportConstruct.cpp
 *
 *  Created on: Sep 21, 2010
 *      Author: Marco.Maneta
 */

#include"Grid.h"
#include"Report.h"


Report::Report(Control &ctrl){

	grid *temp;

	temp = new grid(ctrl.path_BasinFolder + ctrl.fn_rep_mask, ctrl.MapType);

	 UINT4 r, c, data = 0;
	 r = temp->r;
	 c = temp->c;

	 for (UINT4 i=0; i<r; i++)
	    {
	        for (UINT4 j=0; j<c; j++)
	        {
	           data = (UINT4)temp->matrix[i][j];
	        	if( data <= 0 || data >= 32) continue;
	           mask.cells.push_back(cell(i,j,(int)data));
	        }
	    }


	if(temp)
		delete temp;

}
