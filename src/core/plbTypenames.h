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

#ifndef PLB_TYPENAMES_H
#define PLB_TYPENAMES_H

#include <string>

namespace plb {

struct DynamicNativeType {
    virtual ~DynamicNativeType() { }
    virtual DynamicNativeType* clone() const =0;
    virtual int getTypeSize() const =0;
};

template <typename T>
struct NativeType : public DynamicNativeType {
    static char const* getName();
    virtual NativeType<T>* clone() const { return new NativeType<T>(*this); }
    virtual int getTypeSize() const { return sizeof(T); }
};

class NativeTypeConstructor {
public:
    NativeTypeConstructor(std::string typeName);
    ~NativeTypeConstructor();
    NativeTypeConstructor(NativeTypeConstructor const& rhs);
    NativeTypeConstructor& operator=(NativeTypeConstructor const& rhs);
    void swap(NativeTypeConstructor& rhs);
    int getTypeSize() const { return nativeType->getTypeSize(); }
private:
    DynamicNativeType* nativeType;
};

}  // namespace plb

#endif  // PLB_TYPENAMES_H
