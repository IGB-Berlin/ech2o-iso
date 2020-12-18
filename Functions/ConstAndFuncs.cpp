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
 * ConstAndFuncs.cpp
 *
 *  Created on: Feb 23, 2010
 *      Author: Marco Maneta
 */

#include "ConstAndFuncs.h"

//Calculates density of air
//double AirDensity(const double &T /*, const double &P*/){ //TODO: include readings of air pressure to improve calculation of air density
//
//	double P = 101325; // air pressure in Pa
//	return P/(Ra * (T + 273.2)); // air density in Kgm-3
//}

/*Calculates air emissivity using Swinbank(1963) empirical formula with Air Temperature in Celsius*/
double AirEmissivity(const double &AirTemperature){
	return 0.92e-5 * powl(273.2 + AirTemperature, 2);
}

/*calculates the vapor pressure for a given temperature.
T in C and returns vapor pressure in Pa*/
double SatVaporPressure(const double &T){ // T in C and e in Pa

	return 611 * expl((17.3 * T)/(T + 237.3));
 }

//double PsychrometricConst(const double &P, const double &z){ //psychrometric constant air pressure P in Pa
//	//adjust P for elevation as per Allen FAO
//	double Pz = P * powl( ( 293-0.0065*z )/293, 5.26 );
//	return spec_heat_air * Pz / (lat_heat_vap * 0.622); // P in Pa and psychrometric constant in Pa C-1
//}

// Convert isotopic delta to particle ratio, using VSMOW
double Delta2Ratio(const double &di, int iso)
{
  double VSMOW_D = 0.00015576;
  double VSMOW_18O = 0.0020052;
  double Ro = 0;

  // Deuterium
  if (iso == 0)
    Ro = VSMOW_D* (di/1000 + 1);
  // Oxygen 18
  if (iso == 1)
    Ro = VSMOW_18O* (di/1000 + 1);
  // Else(e.g. age)
  if(iso > 1)
    Ro = di;

  return Ro;
}

// Convert isotopic delta to particle ratio, using VSMOW
double Ratio2Delta(const double &Ri, int iso)
{
  double VSMOW_D = 0.00015576;
  double VSMOW_18O = 0.0020052;
  double d = 0;

  // Deuterium
  if (iso == 0)
    d = (Ri/VSMOW_D - 1)* 1000;
  // Oxygen 18
  if (iso == 1)
    d = (Ri/VSMOW_18O - 1)* 1000;
  // Else(e.g. age)
  if(iso > 1)
    d = Ri;

  return d;
}


extern double SoilHeatCapacity(const double &DrySoilHeatCap, const double &Porosity, const double &Theta, const double &SoilTemp){// Dry soil heat capacity in Joules m-3 C-1
	return (1 - Porosity) * DrySoilHeatCap
			+ (Theta) * spec_heat_water * rho_w
			+ (Porosity - Theta) * spec_heat_air * AirDensity(SoilTemp);
}

extern double SoilHeatConductivity(const double &DrySoilHeatCond,
								const double &Porosity,
								const double &Theta){
		return (1 - Porosity) * DrySoilHeatCond
				+ (Theta) * thermal_conduct_water
				+ (Porosity - Theta) * thermal_conduct_air;
}

extern double Calculate_fa(const double &MaxAge, const double &Age){

	double Far;
	double fa;
	Far = std::max<double>(0.0, (Age - 0.2*MaxAge)/MaxAge);
	fa = (Age < 0.2*MaxAge) ?
			0.7 + 0.3 * (Age / (0.2 * MaxAge) ) :
			1.0 / (1.0 + pow((Far/0.95),3));
	return fa;

}

extern double Calculate_ft(const double &Ta, const double &Tmax, const double &Tmin, const double &Topt ){

	double ft;
	if((Ta < Tmin) || (Ta > Tmax))
			return 0;

	ft = std::min<double>(1.0,
			std::max<double>(0.0,
					powl(((Ta - Tmin)/(Topt - Tmin))*((Tmax - Ta)/(Tmax - Topt)),((Tmax - Topt)/(Topt - Tmin)))
			)
		);
	return ft;
}

extern double Calculate_fw(const double &fgc, const double &gsmax, const double &Wr, const double &Wc, const double &Wp){
	double fw;
	fw = std::min<double>(fgc/gsmax, 1/(1 + powl((1-Wr)/Wc, Wp)));
	return fw;

}

extern double Calculate_fd(const double &PropFrostDays){
	return (1 - PropFrostDays);
}

extern double Calculate_gs_light(const double &SolarRadiation, const double &Coeff){

	double f_light;
	f_light = SolarRadiation / ( Coeff + SolarRadiation );
		return(f_light);
}

extern double Calculate_gs_vpd(const double &vpd, const double &Coeff){

		return (expl(-Coeff*vpd));
}

extern double Calculate_gs_lwp(const double &lwp, const double &lwp_max, const double &lwp_min){


	return std::max<double>(0, std::min<double>(1, (lwp_min - lwp)/(lwp_min - lwp_max)));
	//return 1 / (1 + powl(lwp/lwp_min, lwp_max));

}

/*extern double Calculate_gs_theta(const double &theta, const double &fieldcap, const double &wiltpnt, const double &Coeff){

	double Beta;

	Beta = std::max<double>(0, std::min<double>(1, (theta - wiltpnt)/(fieldcap - wiltpnt)));

	return ( 1-expl(-Coeff*Beta) ) / ( 1-expl(-Coeff) ); //from Cox et al 1998 J hydrol "A canopy conductance and photosynt model for use in GCM land surf scheme"

	//return ( Coeff*Beta - powl(Beta,2) );
agui
}*/

extern double rlog(double u_za, double z_a, double z_d, double z_0, int option) {
	double k2 = vonkarman * vonkarman;
		if (option == 0) {
			if (u_za < 0.01)
				u_za = 0.01; //minimum wind velocity to avoid division by zero

			return (powl(log((z_a - z_d) / z_0), 2) / (k2 * u_za));

		} else if (option == 1) {

			return (4.72 * powl(log((z_a - z_d) / z_0), 2) * (1 / (1 + 0.54	* u_za))); //accounts for free convection when u_za->0. From Thom and Oliver 1977 and Jackson et al 1988

		} else {
			std::cout << "Wrong option in the Aerodyn_resist_opt toggle switch. Please verify the configuration file" << std::endl;
			exit(EXIT_FAILURE);

		}
}

extern double Calculate_beta(const double &theta, const double &fieldcap){

	//double Beta;

	if (theta >= fieldcap)
		return 1.0;
	else
		return 0.25*powl( 1 - cos( theta * PI / fieldcap ) , 2 ) ;

	//return ( 1-expl(-Coeff*Beta) ) / ( 1-expl(-Coeff) ); //from Lee and Pielke 1992 or Kondo et al 1990

	//return ( Coeff*Beta - powl(Beta,2) );

}

extern double Calculate_leaf_turnover_TempStress(const double &maxstress, const double &shape, const double &Tcold, const double &T){

	double beta = std::max<double>(0.0, std::min<double>(1, (T - (Tcold - 5))/5 ) );

	return maxstress * powl((1 - beta), shape );

}

extern double Calculate_leaf_turnover_WaterStress(const double &maxstress, const double &shape, const double &UsableWater){

	return maxstress * powl((1 - UsableWater), shape );
}

extern void printProgBar( int percent ){
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}

