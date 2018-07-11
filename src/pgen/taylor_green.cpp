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
  Real M0 = pin->GetReal("problem","M0");
  Real rho0 = pin->GetReal("problem","rho0");
  Real p0 = pin->GetReal("problem","p0");
  Real L = pin->GetReal("mesh","x1max")/PI;
  Real gamma = pin->GetReal("hydro","gamma");
    
  Real c = sqrt(gamma*p0/rho0);
  Real v0 = c*M0;
    
  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx3,nx2,nx1);

  // Initialize density and momentum
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    Real p = p0 + (rho0*v0*v0/16)*(cos(2*pcoord->x1v(i)/L)+cos(2*pcoord->x2v(j)/L))*(cos(2*pcoord->x3v(k)/L)+2);
    Real rho = p*rho0/p0;
    
    phydro->u(IDN,k,j,i) = rho; //density
    phydro->u(IM1,k,j,i) = rho*v0*sin(pcoord->x1v(i)/L)*cos(pcoord->x2v(j)/L)*cos(pcoord->x3v(k)/L); //x momentum
    phydro->u(IM2,k,j,i) = -1*rho*v0*cos(pcoord->x1v(i)/L)*sin(pcoord->x2v(j)/L)*cos(pcoord->x3v(k)/L); //y momentum
    phydro->u(IM3,k,j,i) = 0.0; //z momentum
      
    if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = 0.5*rho*
        ((phydro->u(IM1,k,j,i)/rho)*(phydro->u(IM1,k,j,i)/rho) +
         (phydro->u(IM2,k,j,i)/rho)*(phydro->u(IM2,k,j,i)/rho) +
         (phydro->u(IM3,k,j,i)/rho)*(phydro->u(IM3,k,j,i)/rho)) +
        rho*p/(rho*(gamma-1.0));
        }
    }}}
    
  az.DeleteAthenaArray();
  return;
}
