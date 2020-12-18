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
 * FileRoot.cpp
 *
 *  Created on: Oct 12, 2009
 *      Author: Marco Maneta
 */

#include <string.h>
#include "ParsingFunctions.h"

string FileRoot(string S)
{
    char p[512];
    int k;
    S.erase(S.find('.'), 4);
    strcpy (p, S.c_str());
    k = S.length();
    while(p[k-1]>='0' && p[k-1]<='9') {
      S.erase(k,1);
      k--;
    }
    return S;
}
