//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//
// Edited by Katie Schram 2018
//========================================================================================
//! \file taylor_green.cpp
//  \brief Problem generator for Taylor-Green vortex problem.
//
// REFERENCE:
//========================================================================================

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Orszag-Tang test.  The initial conditions are
//  constructed assuming the domain extends over [-0.5x0.5, -0.5x0.5], so that exact
//  symmetry can be enforced across x=0 and y=0.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real v0 = pin->GetReal("problem","v0");
  Real d0 = pin->GetReal("problem","d0");
  Real p0 = pin->GetReal("problem","p0");
  Real L = pin->GetReal("problem","L");

  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx3,nx2,nx1);

  Real R = 8.314;
  Real T0 = p0/(d0*R);
    
  // Initialize density and momentum
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real p = p0 + (d0*v0*v0/16)*(cos(2*pcoord->x1v(i)/L)+cos(2*pcoord->x2v(j)/L))*(cos(2*pcoord->x3v(k)/L)+2);
    phydro->u(IDN,k,j,i) = p/(R*T0); //density
    phydro->u(IM1,k,j,i) = -d0*v0*cos(pcoord->x1v(i)/L)*sin(pcoord->x2v(j)/L)*cos(pcoord->x3v(k)/L); //y momentum, L = 2pi
    phydro->u(IM2,k,j,i) = d0*v0*sin(pcoord->x1v(i)/(2.0*PI))*cos(pcoord->x2v(j)/(2.0*PI))*cos(pcoord->x3v(k)/L); //x momentum, L = 2pi
    phydro->u(IM3,k,j,i) = 0.0; //z momentum
  }}}

  az.DeleteAthenaArray();
  return;
}
