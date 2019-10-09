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


#include "core/processorIdentifiers2D.h"
#include "core/runTimeDiagnostics.h"
#include <sstream>

namespace plb {

namespace meta {

ProcessorRegistration2D::~ProcessorRegistration2D()
{
    for (pluint iEntry=0; iEntry<processorByNumber.size(); ++iEntry) {
        delete processorByNumber[iEntry].factory;
    }
}

int ProcessorRegistration2D::announce (
        std::string nameOfProcessor,
        ProcessorFactory2D* factory_ )
{
    Entry entry(nameOfProcessor, factory_);
    EntryMap::iterator it = processorByName.find(entry);
    if (it != processorByName.end()) {
        plbLogicError( std::string("The processor ") + nameOfProcessor +
                       std::string(" was registered twice") );
    }
    processorByNumber.push_back(entry);
    int nextId = processorByNumber.size();
    processorByName[entry] = nextId;
    return nextId;
}

int ProcessorRegistration2D::getId(std::string name) const
{
    Entry entry(name, 0);
    EntryMap::const_iterator it = processorByName.find(entry);
    if (it == processorByName.end()) {
        return 0;
    }
    else {
        return it->second;
    }
}

int ProcessorRegistration2D::getNumId() const
{
    return (int)(processorByNumber.size());
}

std::string ProcessorRegistration2D::getName(int id) const
{
    if (id==0) {
        return std::string("Undefined");
    }
    if (id < 0 || id > (int)processorByNumber.size()) {
        std::stringstream message;
        message << "A processor with ID " << id << " doesn't exist.";
        plbLogicError(message.str());
    }
    return processorByNumber[id-1].name;
}

BoxProcessingFunctional2D* ProcessorRegistration2D::create(std::string procName, std::string data)
{
    int id = getId(procName);
    if (id==0) {
        plbLogicError(std::string("A processor with name ")+procName+" does not exits.");
    }
    return processorByNumber[id-1].factory->create(data);
}

ProcessorRegistration2D::EntryMap::const_iterator
    ProcessorRegistration2D::begin() const
{
    return processorByName.begin();
}

ProcessorRegistration2D::EntryMap::const_iterator
    ProcessorRegistration2D::end() const
{
    return processorByName.end();
}


ProcessorRegistration2D& processorRegistration2D() {
    static ProcessorRegistration2D instance;
    return instance;
}

}  // namespace meta

}  // namespace plb
