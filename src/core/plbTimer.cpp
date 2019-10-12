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

#include "parallelism/mpiManager.h"
#include "core/plbTimer.h"
#include <iostream>
#include <map>

#include <ctime>

#ifdef PLB_USE_POSIX
#include <unistd.h>
#endif

namespace plb {

namespace global {

/* ************** Timer ***************************************** */

PlbTimer::PlbTimer()
    : cumulativeTime(0.),
      isOn(false)
{ }

void PlbTimer::start() {
#ifdef PLB_MPI_PARALLEL
    startTime = mpi().getTime();
#else
#if defined PLB_USE_POSIX && defined _POSIX_TIMERS && (_POSIX_TIMERS > 0) && !defined(PLB_NGETTIME)
    clock_gettime(CLOCK_REALTIME, &startTime);
#else
    startClock = clock();
#endif
#endif
    isOn = true;
}

void PlbTimer::restart() {
    reset();
    start();
}

double PlbTimer::stop() {
    cumulativeTime = getTime();
    isOn = false;
    return cumulativeTime;
}

void PlbTimer::reset() {
    cumulativeTime = 0.;
}

double PlbTimer::getTime() const {
    if (isOn) {
#ifdef PLB_MPI_PARALLEL
        return cumulativeTime + mpi().getTime()-startTime;
#else
#if defined PLB_USE_POSIX && defined _POSIX_TIMERS && (_POSIX_TIMERS > 0) && !defined(PLB_NGETTIME)
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        long seconds = ts.tv_sec - startTime.tv_sec; 
        long ns = ts.tv_nsec - startTime.tv_nsec; 

        if (startTime.tv_nsec > ts.tv_nsec) { // clock underflow 
            --seconds; 
            ns += 1e9; 
        } 

        double deltaTime = (double) seconds + (double)ns * (double)1e-9;

        return cumulativeTime + deltaTime;
#else
        return cumulativeTime + (double)(clock()-startClock)
                              / (double)CLOCKS_PER_SEC;
#endif
#endif
    }
    else {
        return cumulativeTime;
    }
}

PlbTimer& timer(std::string nameOfTimer) {
    static std::map<std::string, PlbTimer> timerCollection;
    return timerCollection[nameOfTimer];
}

PlbTimer& plbTimer(std::string nameOfTimer) {
    static std::map<std::string, PlbTimer> timerCollection;
    PlbTimer& answer=timerCollection[nameOfTimer];
    return answer;
}

/* ************** Counter ***************************************** */

PlbCounter::PlbCounter()
    : count(0)
{ }

plint PlbCounter::increment(plint value) {
    count += value;
    return count;
}

void PlbCounter::reset() {
    count = 0;
}

plint PlbCounter::getCount() const {
    return count;
}

PlbCounter& counter(std::string nameOfCounter) {
    static std::map<std::string, PlbCounter> counterCollection;
    return counterCollection[nameOfCounter];
}

PlbCounter& plbCounter(std::string nameOfCounter) {
    static std::map<std::string, PlbCounter> counterCollection;
    return counterCollection[nameOfCounter];
}

}  // namespace global

}  // namespace plb
