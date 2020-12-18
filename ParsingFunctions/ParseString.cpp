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
 * ParseString.cpp
 *
 *  Created on: Oct 12, 2009
 *      Author: Marco Maneta
 */

#include "ParsingFunctions.h"

string ParseString(string FileRoot, int TimeStep)
{
    string S = "0000000000";
    std::stringstream ss;
    ss << TimeStep;
    string time = ss.str();
    S.insert(0, FileRoot);
    S.erase(FileRoot.length()+1, FileRoot.length());
    int length = S.length();
    S.insert((length-time.length()+1), time);
    length = S.length();
    S.erase(S.length()-time.length()+1, time.length());
    S.insert(S.length()-3, ".");

    return S;
}
