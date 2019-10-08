/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2019 FlowKit-Numeca Group Sarl
 * Copyright (C) 2011-2019 University of Geneva
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
 * I/O utilities -- implementation.
 */

#include "core/globalDefs.h"
#include "core/util.h"
#include "io/multiBlockWriter2D.h"
#include "io/multiBlockReader2D.h"
#include "io/imageWriter.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "io/plbFiles.h"
#include "parallelism/mpiManager.h"
#include "io/utilIO_2D.h"

#include <vector>

namespace plb {

void saveState(std::vector<MultiBlock2D*> blocks, plint iteration, bool saveDynamicContent,
        FileName xmlFileName, FileName baseFileName, plint fileNamePadding)
{
    std::string fname_base = createFileName(baseFileName.get(), iteration, fileNamePadding);
    for (pluint i = 0; i < blocks.size(); i++) {
        std::string fname(fname_base+"_"+util::val2str(i));
        parallelIO::save(*blocks[i], fname, saveDynamicContent);
    }
    XMLwriter restart;
    XMLwriter& entry = restart["continue"];
    entry["name"].setString(FileName(fname_base).defaultPath(global::directories().getOutputDir()));
    entry["num_blocks"].set(blocks.size());
    entry["iteration"].set(iteration);
    restart.print(xmlFileName);
}

void loadState(std::vector<MultiBlock2D*> blocks, plint& iteration, bool saveDynamicContent,
        FileName xmlFileName)
{
    XMLreader restart(xmlFileName.get());
    std::string fname_base;
    restart["continue"]["name"].read(fname_base);
    pluint numBlocks;
    restart["continue"]["num_blocks"].read(numBlocks);
    restart["continue"]["iteration"].read(iteration);
    PLB_ASSERT(blocks.size() == numBlocks);
    for (pluint i = 0; i < blocks.size(); i++) {
        std::string fname(fname_base+"_"+util::val2str(i));
        parallelIO::load(fname, *blocks[i], saveDynamicContent);
    }
}

bool abortExecution(FileName abortFileName, std::vector<MultiBlock2D*> blocks, plint iteration,
        bool saveDynamicContent, FileName xmlFileName, FileName baseFileName,
        plint fileNamePadding)
{
    int stop = 0;
    bool stopExecution = false;
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(abortFileName.get().c_str(), "rb");
        if (fp != NULL) {
            stop = 1;
            fclose(fp);
            remove(abortFileName.get().c_str());
        }
    }
    global::mpi().bCast(&stop, 1);
    if (stop) {
        stopExecution = true;
        saveState(blocks, iteration, saveDynamicContent, xmlFileName, baseFileName, fileNamePadding);
    }

    return stopExecution;
}

}  // namespace plb
