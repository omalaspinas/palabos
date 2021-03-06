/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
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

/** \file
 * Base class for scalar, vector and tensor fields for 3D data analysis -- header file.
 */

#ifndef DATA_FIELD_BASE_3D_H
#define DATA_FIELD_BASE_3D_H

#include "core/array.h"
#include "core/block3D.h"
#include "core/dataFieldBase2D.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

namespace plb {

/// Interface for the variants of 3D scalar fields.
template <typename T>
class ScalarFieldBase3D {
public:
    virtual ~ScalarFieldBase3D() { }

public:
    virtual void reset() = 0;
    virtual T &get(plint iX, plint iY, plint iZ) = 0;
    virtual T const &get(plint iX, plint iY, plint iZ) const = 0;
};

/// Interface for the variants of 3D vector and tensor fields.
template <typename T, int nDim>
class TensorFieldBase3D {
public:
    virtual ~TensorFieldBase3D() { }

public:
    virtual void reset() = 0;
    virtual Array<T, nDim> &get(plint iX, plint iY, plint iZ) = 0;
    virtual Array<T, nDim> const &get(plint iX, plint iY, plint iZ) const = 0;
};

/// Interface for the variants of generic-sized 2D vector and tensor fields.
/** The main purpose for these classes is use in dynamically typed languages
 *  like Python. In C++ it's most often better to use the static-sized
 *  TensorField to guarantee type safety.
 */
template <typename T>
class NTensorFieldBase3D {
public:
    NTensorFieldBase3D(int ndim_) : ndim(ndim_) { }
    NTensorFieldBase3D(NTensorFieldBase3D<T> const &rhs) : ndim(rhs.ndim) { }
    void swap(NTensorFieldBase3D &rhs)
    {
        std::swap(ndim, rhs.ndim);
    }
    virtual ~NTensorFieldBase3D() { }

public:
    virtual void reset() = 0;
    virtual T *get(plint iX, plint iY, plint iZ) = 0;
    virtual T const *get(plint iX, plint iY, plint iZ) const = 0;
    plint getNdim() const
    {
        return ndim;
    }

private:
    plint ndim;
};

}  // namespace plb

#endif  // DATA_FIELD_BASE_3D_H
