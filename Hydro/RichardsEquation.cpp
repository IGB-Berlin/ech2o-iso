/*
 * RichardsEquation.cpp
 *
 *  Created on: Nov 3, 2020
 *      Author: marco
 */

//#define ARMA_NO_DEBUG //disables armadillo bound checks for speed optimization
#include <armadillo>
#include "Basin.h"
using namespace arma;

colvec phi_from_head(const colvec &h, const colvec &bc_lambda){
        // Kirchoff transform of permeability function phi = integral_-inf^h k(s)ds

  colvec phi(h);
  phi.elem(find(h < 0)) = arma::exp( bc_lambda.elem(find(h < 0)) % h.elem(find(h < 0)) ) / bc_lambda.elem(find(h < 0));
  phi.elem( arma::find(h >= 0) ) = (1 / bc_lambda.elem(find(h >= 0)) + h.elem(find(h >= 0)));
  return phi;
}

colvec head_from_phi(const colvec &phi, const colvec &bc_lambda){

  uvec pos_a = arma::find(phi < (1 / bc_lambda));
  uvec pos_b = arma::find(phi >= (1 / bc_lambda));

  colvec a = arma::log(bc_lambda % phi) / bc_lambda;
  colvec h = phi - (1 / bc_lambda);

  h.elem(pos_a) = a.elem(pos_a);

  return h;
}

colvec theta_from_phi(const colvec &phi,
		const colvec &bc_lambda,
		const colvec &theta_r,
		const colvec &theta_s){
  colvec h = head_from_phi(clamp(phi, 1e-14, phi.max()), bc_lambda);
  colvec theta = theta_s;
  colvec a = theta_r + (theta_s - theta_r) % exp(bc_lambda % h);
  theta.elem(find(h < 0)) = a.elem(find(h < 0));
  return theta;
}

colvec theta_from_head(const colvec &h,
		const colvec &theta_r,
		const colvec &theta_s,
		const colvec &bc_lambda){
  colvec theta = theta_r + (theta_s - theta_r) % exp(bc_lambda % h);
  return theta;
}

colvec head_from_theta(const colvec &theta,
		const colvec &theta_r,
		const colvec &theta_s,
		const colvec &bc_lambda){
  colvec h =  log( (theta - theta_r) /  (theta_s - theta_r) ) / bc_lambda;
  return h;
}


colvec C(const colvec &phi,
		const colvec &theta_r,
		const colvec &theta_s,
		const colvec &bc_lambda){
  //"""returns transformed capacitance"""
  colvec a = bc_lambda % (theta_s - theta_r);
  colvec C = a;
  C.elem(find( phi >= (1 / bc_lambda))).zeros();
  return C;
}

double C_surf(const double &h, const double &bc_lambda){
  if (h >= 0)
    return 1;
  else
    return 1 / expl(bc_lambda * h);
}

colvec B(const colvec &internode_phi,
		const colvec &internode_Ksat,
		const colvec &internode_lambda){
  //"""Operator B to switch from gravitiy driving flow in unsaturated and saturated conditions"""
  colvec tmp = internode_Ksat % internode_lambda % internode_phi;
  colvec B = internode_Ksat;
  B.elem(find(internode_phi < (1/internode_lambda))) = tmp.elem(find(internode_phi < (1/internode_lambda)));
  return B;
}

colvec q_ds(const colvec &phi,
		const colvec &Ksat,
		const colvec &bc_lambda,
		const double &slope){
  //"""calculates kinematic lateral flow out of soil column"""

  colvec tmp = Ksat % bc_lambda % phi;
  colvec qds = Ksat;
  qds.elem(find( phi < (1/bc_lambda) )) = tmp.elem(find( phi < (1/bc_lambda) ));

  qds = qds * sin(atan(slope));
 
  return qds;
}

colvec Ks_H(const colvec &Ksat){
  //"""returns interblock saturated hydraulic conductivity using harmonic mean"""
  colvec KH = 2 * Ksat % shift(Ksat, -1) / (Ksat + shift(Ksat, -1));

  //# drop mean KH from the bottom boundary
  return KH.head( KH.n_elem - 1 );
}

colvec internode_ave(const colvec &arr){

  colvec temp = (arr + shift(arr, -1)) * 0.5;
  return temp.head(temp.n_elem-1);
}

mat create_Gplus(unsigned int num_nodes, const colvec &internode_distance){
  mat Gplus = diagmat(-ones(num_nodes)) + diagmat(ones((num_nodes - 1)), 1);
  Gplus = Gplus % repelem( (1 / internode_distance), 1, Gplus.n_cols);

  return Gplus;
}

mat create_Gminus(unsigned int num_nodes, const colvec &internode_distance){
  mat Gminus = diagmat(-ones(num_nodes - 1), -1) + diagmat(ones((num_nodes)));

  Gminus = Gminus % repelem((1 / shift(internode_distance,1)), 1, Gminus.n_cols);
  Gminus.row(0) *= 0;

  return Gminus;
}

mat create_diff_op(unsigned int num_nodes){
  mat D = diagmat(-ones(num_nodes)) + diagmat(ones((num_nodes - 1)), -1);
//    D(0,0)=0;
//    D = join_horiz(D, zeros(num_nodes));
//    D(D.n_rows - 1, D.n_cols - 1) = 1;

  return D;
}

colvec pad(const colvec &vec, unsigned int n=1, string type="constant", string side="both", double val = 0){

  colvec temp = ones(vec.size()+n*2) * val;
  temp.subvec(n, n+vec.size()-1) = vec;
  if (type.compare(string("edge")) == 0){
    temp.head(n) = repelem(vec.head(1), n, 1);
    temp.tail(n) = repelem(vec.tail(1), n, 1);
  }

  if (side.compare(string("right")) == 0){
    temp = ones(vec.size()+n) * val;
    temp.head(vec.size()) = vec;
    if (type.compare(string("edge")) == 0)
      temp.tail(n) = repelem(vec.tail(1), n, 1);
  }
  else if (side.compare(string("left")) == 0){
    temp = ones(vec.size()+n) * val;
    temp.tail(vec.size()) = vec;
    if (type.compare(string("edge")) == 0)
      temp.head(n) = repelem(vec.head(1), n, 1);
  }

  return temp;

}

REAL8 q_ovlnd_out(const REAL8 &h, const REAL8 &slope, const REAL8 &n, const REAL8 &dt, const REAL8 &dx){

  //return h/dt;
  if(h<=0)
    return 0;
  else
		//return (sqrt(slope) * powl(h, 2./3))/n * h;
    return std::min<REAL8>(sqrt(slope)/n * powl(h, 5./3), h/dt);

}

void Basin::RichardsEquation(Control &ctrl, double &f,
		double &F,
		double &theta,
		double &theta2,
		double &theta3,
		double &pond,
		double dt,
		int r,
		int c,
		int d) //time step
{

  UINT4 num_nodes = 4; // Number of layers
  UINT4 rind = r;
  UINT4 cind = c;
  REAL8 deltax = _dx;
  REAL8 slope = _slope->matrix[r][c];
  REAL8 leakance = _bedrock_leak->matrix[r][c];
  REAL8 topsoil_ksat = _KsatTopSoil->matrix[r][c];
  //	REAL8 rand_rough = _random_roughness->matrix[r][c] * 100; //random rough is required in cm
  REAL8 mannings_n = _Manningn->matrix[r][c];

  colvec thetan = { theta, theta2, theta3 };
  colvec thetan1 = thetan;
  cout << "theta " << thetan << endl;
  // calculate vertical distance between nodes and add ghost nodes for BC
  colvec deltaz = {1, _depth_layer1->matrix[r][c],
                      _depth_layer2->matrix[r][c],
		      (_soildepth->matrix[r][c] -
		       _depth_layer1->matrix[r][c] -
		       _depth_layer2->matrix[r][c])};
  
  colvec dz_ghost_nodes = pad(deltaz.tail(num_nodes - 1), 1, string("edge"), string("right"));
  colvec internode_distance = (dz_ghost_nodes + shift(dz_ghost_nodes, 1)) * 0.5;
  internode_distance(0) = dz_ghost_nodes(0) * 0.5;
  internode_distance.tail(1) = dz_ghost_nodes.tail(1);

  colvec thetar = ones<colvec>(num_nodes - 1) * _theta_rL1->matrix[r][c];
  colvec thetas = ones<colvec>(num_nodes - 1) * _porosity0->matrix[r][c];
  colvec bclmbda = ones<colvec>(num_nodes) * _BClambda->matrix[r][c];
  //colvec ghost_nodes = pad(bclmbda, 1, "edge");
  colvec internode_lambda = pad(internode_ave(bclmbda), 1, string("edge"), string("right"));

  colvec Ksat = ones<colvec>(thetan.size()) * _Ksat0->matrix[r][c] * _KvKs->matrix[r][c];
  //padd Ksat to include bottom ghost node
  Ksat = pad(Ksat, 1, "edge", "right");
  colvec internode_ksat = pad(internode_ave(Ksat), 1, "edge", "right");

  internode_ksat.head(1) = topsoil_ksat;
  internode_ksat.tail(1) *= leakance;

  // Calculate interblock hydraulic conductivity
  mat Kbarplus = diagmat(internode_ksat);
  mat Kbarminus = diagmat(pad(internode_ksat, 1, "constant", "left", 0).head(num_nodes));

  // Calculate vertical distance between nodes
  colvec dzplus = internode_distance;
  colvec dzminus = shift(internode_distance, 1);
  dzminus(0)=1;
  
  // # Generate gradient operators
  mat Gplus = create_Gplus(num_nodes, internode_distance);
  mat Gminus = create_Gminus(num_nodes, internode_distance);

  // # Generate differential operator for Kbar
  mat D = create_diff_op(num_nodes);
  //calculate heads excluding the overland node
  colvec h = head_from_theta(thetan, thetar, thetas, bclmbda.tail(num_nodes-1));
  
  //add overland flow head
  h = pad(h, 1, "edge", "left");
  h(0) = pond;
  colvec phi = phi_from_head(h, bclmbda);
  
  colvec aug_thetan = zeros(num_nodes);
  aug_thetan.tail(num_nodes-1) = thetan;
  aug_thetan(0) = h(0);

  colvec aug_thetan1 = aug_thetan;

  REAL8 hn = as_scalar(h(0));
  REAL8 hn1 = hn;

  // generate a zero-padded augmented version of phi with one element at the end (right)
  colvec aug_phi = pad(phi, 1, string("constant"), "right");

  // zero vector to hold vertical water fluxes
  colvec qz(num_nodes, fill::zeros);

  colvec Cs(num_nodes, fill::zeros);

  colvec q_ups = _matGWupstreamBC.tube(r,c);
  q_ups /= deltaz;

  colvec q_gw(num_nodes - 1, fill::zeros);
  q_gw = q_ds(phi.tail(num_nodes-1),
			Ksat.tail(num_nodes-1),
			bclmbda.tail(num_nodes-1),
			slope);
  q_gw.insert_rows(0,1);

  if(ctrl.toggle_hydrologic_engine == 2)
    q_gw(0) = q_ovlnd_out(hn1, slope, mannings_n, dt, deltax);

  colvec BCplus(num_nodes, fill::zeros);
  colvec BCminus(num_nodes, fill::zeros);
  
  colvec delta_phi(size(phi), fill::zeros);

  BCminus.head(1) = _incident_water_depth->matrix[r][c]/dt;
  // start Picard iteration

  unsigned int k = 0;
  cout << "phi: "   << phi << endl;
  cout << "h  : "   << h << endl;
  do{
    // Calculate capacitance (dpsi/dthedta)
    Cs.tail(num_nodes-1) = C(phi.tail(num_nodes-1), thetar, thetas, bclmbda.tail(num_nodes-1));
    Cs(0) = C_surf(hn1, as_scalar(bclmbda(0)));
    aug_phi.subvec(0, phi.size()-1) = phi;
    aug_phi.tail(1) = phi.tail(1); //bcplus;

    BCplus.head(1) = 0;
    BCplus.tail(1) = phi.tail(1);

    internode_ksat.head(1) = topsoil_ksat;

    Kbarplus = diagmat(internode_ksat);
    Kbarminus = diagmat(pad(internode_ksat, 1, "constant", "left", 0).head(num_nodes));

    mat lhs = ((1/dt) * diagmat(Cs) -
				repelem((1/deltaz), 1, num_nodes) %
		        (Kbarplus * Gplus -
		        Kbarminus * Gminus));

    colvec rhs = ((1/deltaz) %
		  (
		   (Kbarplus * ((Gplus * phi) + BCplus/dzplus)) -
		   (Kbarminus * (Gminus * phi) - BCminus) +
		   (D * B(internode_ave(aug_phi), internode_ksat, internode_lambda))
		   ) - (q_gw - q_ups)/deltax - (aug_thetan1 - aug_thetan)/dt);

    delta_phi = solve(lhs, rhs);

    phi+=delta_phi;
    //calculate overland flow pressure
    hn1 = as_scalar(head_from_phi(phi, bclmbda)(0));
    cout << "hn1("<<k<<"): " << hn1 << endl;
    // calculate soil moisture at current tension states
    thetan1 = theta_from_phi(phi.tail(num_nodes-1), bclmbda.tail(num_nodes-1), thetar, thetas);
    aug_thetan1.tail(num_nodes - 1) = thetan1;
    aug_thetan1(0) = hn1;

    if(ctrl.toggle_hydrologic_engine == 2)
      q_gw(0) = q_ovlnd_out(hn1, slope, mannings_n, dt, deltax);

    q_gw.tail(num_nodes-1) = q_ds(phi.tail(num_nodes-1),
				Ksat.tail(num_nodes-1),
				bclmbda.tail(num_nodes-1),
				slope);
    qz = -(Kbarplus * ((Gplus * phi) + BCplus/dzplus) - B(internode_ave(aug_phi), internode_ksat, internode_lambda));
    k++;

  }while(norm(delta_phi, "inf") > ctrl.closure_tolerance && k < MAX_ITER);
  if(k>=MAX_ITER)
    cout << "chit" << endl;

  thetan1 = theta_from_phi(phi.tail(num_nodes-1), bclmbda.tail(num_nodes-1), thetar, thetas);
  theta = thetan1(0);
  theta2 = thetan1(1);
  theta3 = thetan1(2);
  _BedrockLeakageFlux->matrix[r][c] = as_scalar(qz.tail(1)); // [m/s]
  if(ctrl.toggle_hydrologic_engine == 2)
    pond = hn1;
  else{
    cout << "hn1 " << hn1 << endl;
    REAL8 pnd = hn1 > 0 ? hn1 : 0;
    q_gw(0) = pnd * deltax / dt;
    pond = hn1 > 0 ? 0 : hn1;
  }

  cout << "Pond "  << pond << endl;
  
  rind = r + rr(d);
  cind = c + cc(d);

  _Disch_old->matrix[r][c] = q_gw(0) * _dx;
  if(d==5){
    _dailyGwtrOutput.cells.push_back(cell(r, c, sum(q_gw.tail(num_nodes-1) % deltaz.tail(num_nodes-1)) * _dx ));
    _dailyOvlndOutput.cells.push_back(cell(r, c, q_gw(0) * _dx));
  }
  else{
    _matGWupstreamBC.tube(rind, cind) += q_gw % deltaz;
  }

  _ponding->matrix[r][c] = pond;
  _incident_water_depth->matrix[r][c] = 0;
} //end richards


