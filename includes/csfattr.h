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
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef CSF__ATTR_H
#define CSF__ATTR_H

#ifndef lint
#define RCS_ID_CSFLEGEN_H "$Header: /Repositories/EcoHydro/includes/csfattr.h,v 1.1.1.1 2009/11/18 19:30:50 mmaneta Exp $"
#endif

#ifdef __cplusplus
 extern "C" {
#endif

typedef enum CSF_ATTR_ID {
        ATTR_ID_LEGEND_V1=1,   /* version 1 legend */
        ATTR_ID_HISTORY=2,     /* history fields */
        ATTR_ID_COLOUR_PAL=3,  /* colour palette */
        ATTR_ID_GREY_PAL=4,    /* grey palette */
        ATTR_ID_DESCRIPTION=5, /* description */
        ATTR_ID_LEGEND_V2=6    /* version 2 legend */
} CSF_ATTR_ID;

#define CSF_LEGEND_ENTRY_SIZE  64
#define CSF_LEGEND_DESCR_SIZE  60

typedef struct CSF_LEGEND {
	INT4    nr;
	char    descr[60];
} CSF_LEGEND;

#ifdef __cplusplus
 }
#endif

#endif /*  CSF__ATTR_H */
