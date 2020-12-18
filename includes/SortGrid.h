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
 * SortedGrid.h
 *
 *  Created on: May 23, 2009
 *      Author: Marco Maneta
 */

#ifndef SORTEDGRID_H_
#define SORTEDGRID_H_

#include "Grid.h"
#include <vector>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
//#define nodata (REAL8)-999.000000000000000000000000000000000

struct cell {
	UINT4 row, col;
	int dir;
	REAL8 val;

	//ctor
	cell() {};
	cell(UINT4 Row, UINT4 Col, int Dir) :
			row(Row), col(Col), dir(Dir) {};
	cell(UINT4 Row, UINT4 Col, REAL8 Val) :
			row(Row), col(Col), val(Val) {};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & row;
		ar & col;
		ar & dir;
		ar & val;
	}
};

struct vectCells {
	UINT4 zone;
	UINT4 numCells;
	vector<cell> cells;
	//ctor
	vectCells() {};
	vectCells(UINT4 Zone) :
			zone(Zone) {};

	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & zone;
		ar & numCells;
		ar & cells;
	}




};

vectCells SortLdd(grid *ldd, vectCells &nonsorted, REAL8 nd);

void loadSortedGrid(vectCells &v, const char *filename);
void saveSortedGrid(vectCells &v, const char *filename);


#endif /* SORTEDGRID_H_ */
