///////////////////////////////////////////////////////////////////////////////
/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact for Palabos:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 * 
 * Contact for npFEM:
 * Christos Kotsalos
 * kotsaloscv@gmail.com
 * Computer Science Department
 * University of Geneva
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
///////////////////////////////////////////////////////////////////////////////
#ifndef CONSTRAINT_FLATTENING_H
#define CONSTRAINT_FLATTENING_H
///////////////////////////////////////////////////////////////////////////////
#include "common.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
struct Constraint_flat 
{
  int ConstraintType_ = -1;
  int idO_ = -1;
  ShapeOpScalar rangeMin_ = 1;
  ShapeOpScalar rangeMax_ = 1;
  ShapeOpScalar Scalar1_ = 1;
  ShapeOpScalar weight_ = 100;
  ShapeOpScalar E_nonePD_ = 0;
  int    *idI_;
  ShapeOpScalar *vectorx_;
  // Matrices are stored in ColMajor order which is also the default of Eigen
  ShapeOpScalar *matrix22_;
  ShapeOpScalar *matrix33_;


  /*
  int    idI_[4];
  ShapeOpScalar vectorx_[4] = {0,0,0,0};
  // Matrices are stored in ColMajor order which is also the default of Eigen
  ShapeOpScalar matrix22_[4] = { 0,0,0,0 };
  ShapeOpScalar matrix33_[9] = { 0,0,0, 0,0,0 ,0,0,0 };
  */
  
  Constraint_flat() 
  {
    ConstraintType_ = -1;
    idO_            = -1;
    rangeMin_       =  1.;
    rangeMax_       =  1.;
    Scalar1_        =  1.;
    weight_         =  100.;
    E_nonePD_       =  0.;
    idI_            = new int    [4];
    vectorx_        = new ShapeOpScalar [4];
    matrix22_       = new ShapeOpScalar [4];
    matrix33_       = new ShapeOpScalar [9];
    
    idI_[0] = -1     ; idI_[1] = -1     ; idI_[2] = -1     ; idI_[3] = -1     ;
    vectorx_[0]  = 0.; vectorx_[1]  = 0.; vectorx_[2]  = 0.; vectorx_[3]  = 0.;
    matrix22_[0] = 0.; matrix22_[1] = 0.; matrix22_[2] = 0.; matrix22_[3] = 0.;
    matrix33_[0] = 0.; matrix33_[1] = 0.; matrix33_[2] = 0.; matrix33_[3] = 0.; 
    matrix33_[4] = 0.; matrix33_[5] = 0.; matrix33_[6] = 0.; matrix33_[7] = 0.; 
    matrix33_[8] = 0.;
  }

  ~Constraint_flat() 
  {
    delete idI_;
    delete vectorx_;
    delete matrix22_;
    delete matrix33_;
  }
  
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // CONSTRAINT_FLATTENING_H
///////////////////////////////////////////////////////////////////////////////