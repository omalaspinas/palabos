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

#ifndef REAL_FUNCTION_3D_H
#define REAL_FUNCTION_3D_H

#include <cmath>

#include "core/array.h"
#include "core/functions.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockLattice3D.hh"

namespace plb {

template <typename T>
struct RealFunction3D {
    virtual ~RealFunction3D() { }
    virtual T operator()(Array<T, 3> const &position) const = 0;
    virtual RealFunction3D<T> *clone() const = 0;
};

template <typename T>
class ConstantRealFunction3D : public RealFunction3D<T> {
public:
    ConstantRealFunction3D(T constantValue_) : constantValue(constantValue_) { }

    virtual T operator()(Array<T, 3> const &) const
    {
        return constantValue;
    }

    virtual ConstantRealFunction3D<T> *clone() const
    {
        return new ConstantRealFunction3D<T>(*this);
    }

private:
    T constantValue;
};

template <typename T, template <typename U> class Descriptor>
class IncreasingRealFunction3D : public RealFunction3D<T> {
public:
    IncreasingRealFunction3D(
        T maxValue_, MultiBlockLattice3D<T, Descriptor> const &lattice_, T tOffset_, T maxT_) :
        maxValue(maxValue_), lattice(lattice_), tOffset(tOffset_), maxT(maxT_)
    { }

    virtual T operator()(Array<T, 3> const &) const
    {
        T t = (T)lattice.getTimeCounter().getTime() + tOffset;
        T value = util::sinIncreasingFunction<T>(t, maxT) * maxValue;
        return value;
    }

    virtual IncreasingRealFunction3D<T, Descriptor> *clone() const
    {
        return new IncreasingRealFunction3D<T, Descriptor>(*this);
    }

private:
    T maxValue;
    MultiBlockLattice3D<T, Descriptor> const &lattice;
    T tOffset, maxT;
};

template <typename T, template <typename U> class Descriptor>
class HarmonicRealFunction3D : public RealFunction3D<T> {
public:
    HarmonicRealFunction3D(
        T amplitude_, T angularFrequency_, T phase_,
        MultiBlockLattice3D<T, Descriptor> const &lattice_, T tOffset_) :
        amplitude(amplitude_),
        angularFrequency(angularFrequency_),
        phase(phase_),
        lattice(lattice_),
        tOffset(tOffset_)
    { }

    virtual T operator()(Array<T, 3> const &) const
    {
        T t = (T)lattice.getTimeCounter().getTime() + tOffset;
        return amplitude * std::cos(angularFrequency * t + phase);
    }

    virtual HarmonicRealFunction3D<T, Descriptor> *clone() const
    {
        return new HarmonicRealFunction3D<T, Descriptor>(*this);
    }

private:
    T amplitude, angularFrequency, phase;
    MultiBlockLattice3D<T, Descriptor> const &lattice;
    T tOffset;
};

}  // namespace plb

#endif  // REAL_FUNCTION_3D_H
