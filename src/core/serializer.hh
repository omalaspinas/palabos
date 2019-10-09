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


#ifndef SERIALIZER_HH
#define SERIALIZER_HH

#include "parallelism/mpiManager.h"
#include "core/serializer.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {


/* *************** Class WriteToSerialArray ********************************* */

template<typename T>
WriteToSerialArray<T>::WriteToSerialArray(T* array_, plint typedSize_)
    : array((char*)(array_)),
      typedSize(typedSize_),
      pos(0)
{ }

template<typename T>
WriteToSerialArray<T>* WriteToSerialArray<T>::clone() const {
    return new WriteToSerialArray<T>(*this);
}

template<typename T>
void WriteToSerialArray<T>::writeHeader(pluint dataSize)
{ }

template<typename T>
void WriteToSerialArray<T>::writeData (
        char const* dataBuffer, pluint bufferSize )
{
    plint charSize = (plint)bufferSize;
    std::copy(dataBuffer, dataBuffer+charSize, array+pos);
    pos += charSize;
}


/* *************** Class ReadFromSerialArray ********************************* */

template<typename T>
ReadFromSerialArray<T>::ReadFromSerialArray(T const* array_, plint typedSize_)
    : array((char const*)(array_)),
      typedSize(typedSize_),
      pos(0)
{ }

template<typename T>
ReadFromSerialArray<T>* ReadFromSerialArray<T>::clone() const {
    return new ReadFromSerialArray<T>(*this);
}

template<typename T>
void ReadFromSerialArray<T>::readHeader(pluint dataSize) const
{ }

template<typename T>
void ReadFromSerialArray<T>::readData(char* dataBuffer, pluint bufferSize) const
{
    plint charSize = (plint)bufferSize;
    std::copy(array+pos, array+pos+charSize, dataBuffer);
    pos += charSize;
}

}  // namespace plb

#endif  // SERIALIZER_HH
