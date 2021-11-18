#ifndef CLOT_H
#define CLOT_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "simParams.h"
#include <vector>
#include <cstdlib>      /* srand, rand */
#include <ctime>        /* time */
#include <random>

#define DESCRIPTOR descriptors::RhoBarJD3Q19Descriptor
extern SimulationParameters simParam;

using namespace plb;
using namespace std;

// Data processor to apply the solid fraction to the cell density (that will be used for partial collision)
template< typename T,template<typename U> class Descriptor>
class ApplySolidFractionToRhoBar : public BoxProcessingFunctional3D 
{
public:
  virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks)
  {
  	PLB_ASSERT(atomicBlocks.size() == 2);
    BlockLattice3D<T,Descriptor>& lattice = *dynamic_cast<BlockLattice3D<T,Descriptor> *>(atomicBlocks[0]);
    ScalarField3D<T>& solidFraction = *dynamic_cast<ScalarField3D<T>*>(atomicBlocks[1]);
    Dot3D offsetSF = computeRelativeDisplacement(lattice, solidFraction);
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint posX = iX + offsetSF.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint posY = iY + offsetSF.y;
        	for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            	plint posZ = iZ + offsetSF.z;

                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt) 
                = solidFraction.get(posX,posY,posZ); // apply here the solid fraction value to the lattice density
            }
        }
    }
  }

  virtual ApplySolidFractionToRhoBar<T,Descriptor>* clone() const
  {
  	return new ApplySolidFractionToRhoBar<T,Descriptor>(*this);
  }
  virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
  	modified[0] = modif::dataStructure;
      modified[1] = modif::nothing;
  }
};

  // data processor that generates the shape of the clot for analytical cylinder geometry
template<typename T, template<typename U> class Descriptor>
class readClotAnaFromFile : public BoxProcessingFunctional3D
{
public:
  readClotAnaFromFile(double frac_) : frac(frac_) {}
  virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields)
  {
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<plint>& clotFlags = *dynamic_cast<ScalarField3D<plint>*>(fields[0]);
    ScalarField3D<T>& clotSolidFraction = *dynamic_cast<ScalarField3D<T>*>(fields[1]);
    ScalarField3D<T>& clotSolidFractionPhys = *dynamic_cast<ScalarField3D<T>*>(fields[2]);

    random_device rd;       
    // random-number engine used (Mersenne-Twister in this case)
    mt19937 rng(rd());		// comment '(rd())' to use always the same seed (for graphs) 
    Dot3D offset = clotFlags.getLocation();
    Dot3D offsetSolidFraction = computeRelativeDisplacement(clotFlags,clotSolidFraction);
    Dot3D offsetSolidFractionPhys = computeRelativeDisplacement(clotFlags,clotSolidFractionPhys);
    
    plint totalFnQty = 0;
    plint totalNbVoxels = 0;
  	std::ifstream clotGeoFile;
  	clotGeoFile.open(simParam.clotGeometryFile.c_str(),std::fstream::in);
  	std::string line;
    std::string valsF;  

    for(plint iZ=domain.z0; iZ<=domain.z1;++iZ)
    {
    	// position the reader at correct Z slice
    	// iZ + offset.z : absolute z position ; read until position of z slice starting for this core (disregard previous z slices)
    	for (size_t lineNb = 0 ; lineNb < int(iZ+offset.z-simParam.clotBeginZ*simParam.nz)*(countNXY(simParam.clotGeometryFile)) 
    		&& line == ""; ++lineNb)
    		getline(clotGeoFile,line);
    	if (line == "")
    		getline(clotGeoFile,line);
    	// position reader at coordinates relevant for this atomic block (in terms of X and Y domains)	
    	for (size_t lineNb = 0; lineNb <= domain.y0+offset.y; ++lineNb)
    		getline(clotGeoFile,line);
    	

    	for(plint iY=domain.y0; iY <= domain.y1; ++iY)
    	{
		    if (getline(clotGeoFile,line) && line!="")
		    {
			    std::stringstream ss(line);

  				for (size_t tokenNb = 0;tokenNb <= domain.x0+offset.x; ++tokenNb)
  			    	ss >> valsF;
	    		for(plint iX=domain.x0; iX <= domain.x1; ++iX)
				  {
            ss >> valsF;
				    double ns;
				    // apply here transformation ns(ns*). Now using Davies equation
				    if (std::stod(valsF) == 0.)
				    	ns = 0;
				    else
            {
				    	if (simParam.permeModel == "Davies")
  							ns = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)/(16*pow(std::stod(valsF),1.5)*(1+56*pow(std::stod(valsF),3.)))/(simParam.nu_LB *pow(simParam.dx,2.)));
  						else if (simParam.permeModel == "Clague")
  							ns = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)*0.50941*pow(pow(simParam.pi/T(4.*std::stod(valsF)), 0.5)-1., 2.)*exp(-1.8042*std::stod(valsF)) /(simParam.nu_LB *pow(simParam.dx,2.)));
  						else if (simParam.permeModel == "JJ")
  							ns = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)*(3./T(20*std::stod(valsF)))*(-log(std::stod(valsF))-0.931)/(simParam.nu_LB *pow(simParam.dx,2.)));
  						else	// Walsh by default, no rescaling
  							ns = std::stod(valsF);
				    }
				    clotSolidFraction.get(iX+offsetSolidFraction.x,iY+offsetSolidFraction.y,iZ+offsetSolidFraction.z) = ns*frac;
				    clotSolidFractionPhys.get(iX+offsetSolidFractionPhys.x,iY+offsetSolidFractionPhys.y,iZ+offsetSolidFractionPhys.z) = std::stod(valsF)*frac;
				    
  					clotFlags.get(iX,iY,iZ) = plint(simParam.p0*std::stod(valsF)*frac*simParam.dx*simParam.dx*simParam.dx*1e27*(2./double(45.)));
  					totalFnQty += clotFlags.get(iX,iY,iZ);
  					totalNbVoxels += 1;				
					
  				}
  			}
    	}
    	while (line !="")
    	// flush out the empty line between Z slices
    		getline(clotGeoFile,line);
    }
  }
  void getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
  	modified[0] = modif::staticVariables;
  	modified[1] = modif::staticVariables;
  	modified[2] = modif::staticVariables;
  }
  virtual readClotAnaFromFile<T,Descriptor>* clone() const  {
    return new readClotAnaFromFile<T,Descriptor>(*this);
  }
private:
	double frac;
};
////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////// 
//Data processor that computes reaction between anti-fibrin and fibrin
template<typename T, template<typename U> class Descriptor>
class clotParticleInteraction : public BoxProcessingFunctional3D
{
public:
	clotParticleInteraction(plint currentIter_): currentIter(currentIter_) {}
	virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields)
	{
    PLB_PRECONDITION( fields.size()==4 );
		ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(fields[0]);
    ScalarField3D<plint>& clotFlags = *dynamic_cast<ScalarField3D<plint>*>(fields[1]);
  	ScalarField3D<T>& clotSolidFraction = *dynamic_cast<ScalarField3D<T>*>(fields[2]);
  	ScalarField3D<T>& clotSolidFractionPhys = *dynamic_cast<ScalarField3D<T>*>(fields[3]);
    Dot3D offsetFlags = computeRelativeDisplacement(particleField, clotFlags);
    Dot3D offsetSF = computeRelativeDisplacement(particleField, clotSolidFraction);
    Dot3D offsetSFPhys = computeRelativeDisplacement(particleField, clotSolidFractionPhys);

    for(plint iX=domain.x0; iX <= domain.x1; ++iX)
    {
    	for(plint iY=domain.y0; iY <= domain.y1; ++iY)
    	{
  			for(plint iZ = domain.z0; iZ<=domain.z1;++iZ)
  			{
  				if(clotFlags.get(iX+offsetFlags.x, iY+offsetFlags.y, iZ+offsetFlags.z) > 0)	// if current position is a clot (flag > 0) 
  				{
  					// ------------- "t" ----------------
    				std::vector<Particle3D<T,Descriptor>*> foundParticles;		// vector that will contain the particles found in the cellDomain
          	Box3D cellDomain;
          	if (simParam.interactNeighbor)
          		cellDomain = Box3D(iX-1,iX+1, iY-1,iY+1, iZ-1,iZ+1);						// cellDomain to check = current position
          	else
          		cellDomain = Box3D(iX,iX, iY,iY, iZ,iZ);						// cellDomain to check = current position
          	particleField.findParticles(cellDomain, foundParticles);	// find the particles (if any) in cellDomain
      		
        		if (foundParticles.size() != 0)
            {
          		plint effectiveAntiFnQty = 0;		// compute the effective nb of real particles at position
          		for (size_t i = 0 ; i < foundParticles.size() ; ++i)
          			effectiveAntiFnQty += foundParticles[i]->getTag();		// tag contains remaining real nb of particles for each super particle
          		plint effectiveFnQty = clotFlags.get(iX+offsetFlags.x, iY+offsetFlags.y, iZ+offsetFlags.z);		// get the effective nb of clot particles at position

          		plint effectiveFnToRemove = util::roundToInt(simParam.dt*simParam.k1*(simParam.dx*simParam.dx*simParam.dx)*effectiveFnQty*effectiveAntiFnQty); 		// compute reaction equation for Fn (divide by V_elem to convert concentration to qty)
          		plint effectiveAntiFnToRemove = util::roundToInt(simParam.dt*simParam.k2*(simParam.dx*simParam.dx*simParam.dx)*effectiveFnQty*effectiveAntiFnQty);		// compute reaction equation for antiFn
          		
          		// ------------- "t+1" --------------
          		// apply eq for Fn
          		if (effectiveFnQty > effectiveFnToRemove)
          			clotFlags.get(iX+offsetFlags.x, iY+offsetFlags.y, iZ+offsetFlags.z) = effectiveFnQty-effectiveFnToRemove;
          		else 
              {
          			clotFlags.get(iX+offsetFlags.x, iY+offsetFlags.y, iZ+offsetFlags.z) = flagNums::fibrinDestroyed;
    						clotSolidFraction.get(iX+offsetSF.x, iY+offsetSF.y, iZ+offsetSF.z) = 0;	// set to fluid node
    						clotSolidFractionPhys.get(iX+offsetSFPhys.x, iY+offsetSFPhys.y, iZ+offsetSFPhys.z) = 0;	// set to fluid node
          		}
          		// apply eq for anti-Fn
        			for (size_t i = 0 ; (i < foundParticles.size() && effectiveAntiFnToRemove > 0) ; ++i)		// check all the particles found
        			{	
        				if (effectiveAntiFnToRemove < foundParticles[i]->getTag())
                {
        					foundParticles[i]->setTag(foundParticles[i]->getTag() - effectiveAntiFnToRemove);
        					effectiveAntiFnToRemove = 0;
        				}

        				else	// if more effectiveAntiFnToRemove than antiFn available in i-th particle : remove that particle and update effAntToRemove
        				{

                        // update qtyAntiFnToRemove
                      effectiveAntiFnToRemove -= foundParticles[i]->getTag();
                      // set to 0 remaining nb of eff particles to the super-particle
                      foundParticles[i]->setTag(0);
                      	// destroy super-particles which tag (effNb) == 0
                      particleField.removeParticles(cellDomain, 0);
                      // update foundParticles' size for the loop
                      particleField.findParticles(cellDomain, foundParticles);
        				}
        			}
        		} // if particles found
        	} // if at a clot
  			} // for Z
    	} // for Y	
    } // for X

	} // processGenericBlocks

	void getTypeOfModification(std::vector<modif::ModifT>& modified) const
	{
		modified[0] = modif::dynamicVariables;	// particle field
		modified[1] = modif::staticVariables;	// clot flags
		modified[2] = modif::staticVariables;	// clot solid fraction
		modified[3] = modif::staticVariables;	// clot solid fraction phys
	}
	virtual clotParticleInteraction<T,Descriptor>* clone() const  {
		return new clotParticleInteraction<T,Descriptor>(*this);
	}
private:
	plint currentIter;
};
////////////////////////////////////////////////////////////////////////

// data processor that updates the solid fraction after reaction, based on the permeability model
template<typename T, template<typename U> class Descriptor>
class updatesF_R_L_x_t : public BoxProcessingFunctional3D
{
public:
	updatesF_R_L_x_t() {}
  virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields)
  {
    PLB_PRECONDITION( fields.size()==6 );
    ScalarField3D<double>& sF_x_t = *dynamic_cast<ScalarField3D<double>*>(fields[0]);
    ScalarField3D<double>& R_x_t = *dynamic_cast<ScalarField3D<double>*>(fields[1]);
    ScalarField3D<double>& L_x_t = *dynamic_cast<ScalarField3D<double>*>(fields[2]);
    ScalarField3D<plint>& clotFlagsIni = *dynamic_cast<ScalarField3D<plint>*>(fields[3]);
    ScalarField3D<plint>& clotFlagsCurr = *dynamic_cast<ScalarField3D<plint>*>(fields[4]);
    ScalarField3D<double>& sFPhys_x_t = *dynamic_cast<ScalarField3D<double>*>(fields[5]);

    Dot3D offsetR = computeRelativeDisplacement(sF_x_t, R_x_t);
    Dot3D offsetL = computeRelativeDisplacement(sF_x_t, L_x_t);
    Dot3D offsetFlagsIni = computeRelativeDisplacement(sF_x_t, clotFlagsIni);
    Dot3D offsetFlagsCurr = computeRelativeDisplacement(sF_x_t, clotFlagsCurr);
    Dot3D offsetsFPhys = computeRelativeDisplacement(sF_x_t, sFPhys_x_t);
    T sFTmp;
    
    for(plint iX=domain.x0; iX <= domain.x1; ++iX)
    {
    	for(plint iY=domain.y0; iY <= domain.y1; ++iY)
    	{
  			for(plint iZ=domain.z0; iZ<=domain.z1;++iZ)
  			{
  				// lysed quantity in this voxel, in µM
  				L_x_t.get(iX+offsetL.x,iY+offsetL.y,iZ+offsetL.z) = (clotFlagsIni.get(iX+offsetFlagsIni.x,iY+offsetFlagsIni.y,iZ+offsetFlagsIni.z) - clotFlagsCurr.get(iX+offsetFlagsCurr.x,iY+offsetFlagsCurr.y,iZ+offsetFlagsCurr.z))
  									  * 1e3 / double(simParam.dx*simParam.dx*simParam.dx*simParam.avogadro);		// convert from qty of matter to µM

  				// L_DV/DV = 0 if clot empty (clotFlagsCurr = 0), then sF = 0
  			  	double LDV_DV = clotFlagsIni.get(iX+offsetFlagsCurr.x,iY+offsetFlagsCurr.y,iZ+offsetFlagsCurr.z)
  			  			 / double(simParam.dx*simParam.dx*simParam.dx*1e6 * simParam.p0 * (2./(double)45) * simParam.pi * simParam.Rf0*simParam.Rf0);
  			  			 	// DV (cm^3)
  	  			 
  	  			// radius of fiber in clot voxel at pos x, current time
  			  	// check first if this was a clot position, to avoid division by 0 (LDV_DV == 0)
  			  	if (clotFlagsIni.get(iX+offsetFlagsIni.x,iY+offsetFlagsIni.y,iZ+offsetFlagsIni.z) == 0)
  			  		R_x_t.get(iX+offsetR.x,iY+offsetR.y,iZ+offsetR.z) = 0;
  		  		else
  					R_x_t.get(iX+offsetR.x,iY+offsetR.y,iZ+offsetR.z) = sqrt( simParam.Rf0*simParam.Rf0 - double(1./double(simParam.pi * simParam.p0)) *
  																		L_x_t.get(iX+offsetL.x,iY+offsetL.y,iZ+offsetL.z) / 
  																		((LDV_DV) * double(1./double(simParam.avogadro)) * ((double)2e6/(double)45e-3))
  																		);


    				// solid fraction
    				sFTmp = simParam.pi*pow(R_x_t.get(iX+offsetR.x,iY+offsetR.y,iZ+offsetR.z),2) * (LDV_DV) * 1e-21;
    				sFPhys_x_t.get(iX+offsetsFPhys.x,iY+offsetsFPhys.y,iZ+offsetsFPhys.z) = sFTmp;
    				// using Davies equation to compute ns(ns*_new)
    				if (sFTmp == 0.)
    					sF_x_t.get(iX,iY,iZ) = 0.;
    				else{
    					if (simParam.permeModel == "Davies")
    						sF_x_t.get(iX,iY,iZ) = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)/(16*pow(sFTmp,1.5)*(1+56*pow(sFTmp,3.)))/(simParam.nu_LB *pow(simParam.dx,2.)));
    					else if (simParam.permeModel == "Clague")
    						sF_x_t.get(iX,iY,iZ) = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)*0.50941*pow(pow(simParam.pi/T(4.*sFTmp), 0.5)-1., 2.)*exp(-1.8042*sFTmp) /(simParam.nu_LB *pow(simParam.dx,2.)));
    					else if (simParam.permeModel == "JJ")
    						sF_x_t.get(iX,iY,iZ) = T(1.)/T(1.+2.*pow(simParam.Rf0*1.e-9,2.)*(3./T(20*sFTmp))*(-log(sFTmp)-0.931)/(simParam.nu_LB *pow(simParam.dx,2.)));
    					else	// Walsh by default, no rescaling
    						sF_x_t.get(iX,iY,iZ) = sFTmp;
  				}
  		
  			}
    	}
    }
  }

  void getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
  	modified[0] = modif::staticVariables;
  	modified[1] = modif::staticVariables;
  	modified[2] = modif::staticVariables;
  	modified[3] = modif::nothing;
  	modified[4] = modif::nothing;
  	modified[5] = modif::nothing;
  }
  virtual updatesF_R_L_x_t<T,Descriptor>* clone() const  {
    return new updatesF_R_L_x_t<T,Descriptor>(*this);
  }
  };

// Class Clot
typedef double T;
class Clot
{
public:
	Clot(){}

	Clot(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_, plint zBegin_ = 0.0,
	 plint zEnd_ = 1.0, T porosity_ = 0.5){
    lattice = &lattice_;
    zBegin = zBegin_;
    zEnd = zEnd_;
    porosity = porosity_;
  }

	Clot(const Clot & clot){
    lattice = clot.lattice;
    zBegin = clot.zBegin;
    zEnd = clot.zEnd;
    porosity = clot.porosity;
  }

  // wrapper function that generates the clot
	void generateFibrin(MultiScalarField3D<plint> & flags, MultiScalarField3D<double> & clotSolidFraction, 
		MultiScalarField3D<double> & clotSolidFractionPhys, double frac){
  
    pcout << "Creating clot... "<< frac*100 << " % of final sF" << endl;                    
    
    // to prevent from making several times the call to getBoundingBox()
    Box3D latticeDomain = lattice->getBoundingBox();
    // define domain before data processor, so parallelization can be made properly (domain divided to cores)
    Box3D clotDomain(latticeDomain.x0,latticeDomain.x1, latticeDomain.y0,latticeDomain.y1, zBegin,zEnd);

    std::vector<MultiBlock3D*> clotFlagArg;
    clotFlagArg.push_back(&flags);
    clotFlagArg.push_back(&clotSolidFraction);   
    clotFlagArg.push_back(&clotSolidFractionPhys);    
    applyProcessingFunctional(new readClotAnaFromFile<double, DESCRIPTOR>(frac), clotDomain, clotFlagArg);
    
    std::vector<MultiBlock3D*> latticeClotArg;
    latticeClotArg.push_back(lattice);
    latticeClotArg.push_back(&clotSolidFraction);   
    applyProcessingFunctional(new ApplySolidFractionToRhoBar<double, DESCRIPTOR>(), clotDomain, latticeClotArg);
      
  }

  // wrapper function that calls the interaction between anti-fibrin and fibrin
	void interact(std::vector<MultiBlock3D*> particleFlagSFArg, plint iter){
    // to prevent from making several times the call to getBoundingBox()
    Box3D latticeDomain = lattice->getBoundingBox();
    // define domain before data processor, so parallelization can be made properly (domain divided to cores)
    Box3D clotDomain(latticeDomain.x0,latticeDomain.x1, latticeDomain.y0,latticeDomain.y1, zBegin,zEnd);  

    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR>>& particleField = *dynamic_cast<MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR>>*>(particleFlagSFArg[0]);
    pluint initNbParticles = countParticles(particleField, particleField.getBoundingBox());

    // do the particle/clot interaction
    pluint nbDestroyed = 0;
    applyProcessingFunctional(new clotParticleInteraction<double, DESCRIPTOR>(iter), clotDomain, particleFlagSFArg);

    nbDestroyed = initNbParticles - countParticles(particleField, particleField.getBoundingBox());
    
    if (nbDestroyed > 0)
      pcout << nbDestroyed << " super particle(s) destroyed at time " << iter*simParam.dt << "(it=" << iter << ")..." << endl;
    

    MultiScalarField3D<T>& clotSolidFraction = *dynamic_cast<MultiScalarField3D<T>*>(particleFlagSFArg[2]);
    std::vector<MultiBlock3D*> latticeClotArg;
    latticeClotArg.push_back(lattice);
    latticeClotArg.push_back(&clotSolidFraction);   
    applyProcessingFunctional(new ApplySolidFractionToRhoBar<double, DESCRIPTOR>(), clotDomain, latticeClotArg);
  }

private:
	MultiBlockLattice3D<T,DESCRIPTOR>* lattice;
	T porosity;
	plint zBegin;
	plint zEnd;

};

// count the length of the clot in z-direction based on the clot txt file
plint countClotZLength()
{
	std::fstream clotFile;
    clotFile.open(simParam.clotGeometryFile.c_str(), std::fstream::in);

    std::string line;

    size_t count=0;
    while(!clotFile.eof()){ // ss is used more like cin
    	getline(clotFile,line);
    	if (line.empty())
        	count++;
    }

    clotFile.close();
    return count-1;
}

#endif
