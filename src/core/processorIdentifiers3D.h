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

#ifndef PROCESSOR_IDENTIFIERS_3D_H
#define PROCESSOR_IDENTIFIERS_3D_H

#include <map>
#include <string>
#include <vector>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "core/plbTypenames.h"
#include "core/util.h"

namespace plb {

namespace meta {

struct ProcessorFactory3D {
    virtual ~ProcessorFactory3D() { }
    virtual BoxProcessingFunctional3D *create(std::string data) const = 0;
};

class ProcessorRegistration3D {
public:
    struct Entry {
        Entry(std::string name_, ProcessorFactory3D *factory_) : name(name_), factory(factory_) { }
        std::string name;
        ProcessorFactory3D *factory;
    };
    struct EntryLessThan {
        bool operator()(Entry const &entry1, Entry const &entry2) const
        {
            return entry1.name < entry2.name;
        }
    };
    typedef std::map<Entry, int, EntryLessThan> EntryMap;

public:
    ProcessorRegistration3D() { }
    ~ProcessorRegistration3D();
    ProcessorRegistration3D(ProcessorRegistration3D const &) = delete;
    ProcessorRegistration3D &operator=(ProcessorRegistration3D const &) = delete;
    int announce(std::string nameOfProcessor, ProcessorFactory3D *factory_ = 0);
    int getId(std::string name) const;
    int getNumId() const;
    std::string getName(int id) const;
    BoxProcessingFunctional3D *create(std::string procName, std::string data);
    EntryMap::const_iterator begin() const;
    EntryMap::const_iterator end() const;

private:
    EntryMap processorByName;
    std::vector<Entry> processorByNumber;
};

ProcessorRegistration3D &processorRegistration3D();

template <class FunctionalT>
class FunctionalFactory3D : public ProcessorFactory3D {
    virtual BoxProcessingFunctional3D *create(std::string data) const
    {
        BoxProcessingFunctional3D *functional = new FunctionalT;
        functional->unserialize(data);
        return functional;
    }
};

template <class FunctionalT>
int registerProcessor3D(std::string name)
{
    return processorRegistration3D().announce(name, new FunctionalFactory3D<FunctionalT>);
}

template <class FunctionalT, typename T>
int registerProcessor3D(std::string name)
{
    std::string fullName = name + "_" + std::string(NativeType<T>::getName());
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

template <class FunctionalT, typename T1, typename T2>
int registerProcessor3D(std::string name)
{
    std::string fullName = name + "_" + std::string(NativeType<T1>::getName()) + "_"
                           + std::string(NativeType<T2>::getName());
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

template <class FunctionalT, typename T, int val>
int registerProcessor3D(std::string name)
{
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" + util::val2str(val);
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

template <class FunctionalT, typename T, template <typename U> class Descriptor>
int registerProcessor3D(std::string name)
{
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" + std::string(Descriptor<T>::name);
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

template <class FunctionalT, typename T, template <typename U> class Descriptor, int val1, int val2>
int registerProcessor3D(std::string name)
{
    std::string fullName = name + "_" + std::string(NativeType<T>::getName()) + "_"
                           + std::string(Descriptor<T>::name) + "_" + util::val2str(val1) + "_"
                           + util::val2str(val2);
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

template <
    class FunctionalT, typename T, template <typename U> class Descriptor, int val1, int val2,
    int val3>
int registerProcessor3D(std::string name)
{
    std::string fullName = name + "_" + std::string(NativeType<T>::getName()) + "_"
                           + std::string(Descriptor<T>::name) + "_" + util::val2str(val1) + "_"
                           + util::val2str(val2) + "_" + util::val2str(val3);
    return processorRegistration3D().announce(fullName, new FunctionalFactory3D<FunctionalT>);
}

}  // namespace meta

}  // namespace plb

#endif  // PROCESSOR_IDENTIFIERS_3D_H
