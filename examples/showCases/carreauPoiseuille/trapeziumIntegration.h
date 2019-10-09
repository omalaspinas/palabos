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

#include <cmath>

#include "functions.h"

#ifndef TRAPEZIUM_INTEGRATION_H
#define TRAPEZIUM_INTEGRATION_H


/** This class is a Newton-Raphson root finding solver.
 * It find y such that F(y)=0. It takes the function type as template
 * parameter and F(y), F'(y), the tolerance and the maximum iterations
 * the user allows in order for the NR solver to find the root.
 */

template <typename T>
class TrapeziumIntegration
{
public :
	TrapeziumIntegration(Function<T> *function_, T y0_, int numberSteps_)
		: function(function_), y0(y0_), numberSteps(numberSteps_)
	{   }
	
	T operator()(T y) const 
	{	
		T dy = (y-y0)/(T)numberSteps;
		
		if (dy <= T() || y==y0)
		{
			return T();
		}
		T integral = ((*(function))(y0)+(*(function))(y))/(T)2;
// 		std::cout << y << ", " << integral << std::endl;
		for (int iS = 1; iS < numberSteps; ++iS)
		{
// 			std::cout << integral << std::endl;
			T ty = y0+(T)iS*dy;
			integral += (*(function))(ty);
		}
		integral *= dy;
        
// 		std::cout << y << ", " << integral << std::endl;
		return integral;
    }
	
private :
	Function<T> *function;
	T y0;
	int numberSteps;
};

#endif
