/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo                        | 
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef ANALYZER_FORCE_H
#define ANALYZER_FORCE_H
/**
* @file
*
* @class AnalyzerForce
*
* @brief Calculate force of chain monomers by counting propability of jumps up and
* down.
* 
* @details Calculate the force acting on a set of monomers in a predefined environment (ingredients) 
* using FeatureMoleculesIO and FeatureExcludedVolume.
* Moves are checked in z direction
*
* @tparam IngredientsType
**/

#include <iostream>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>



template<class IngredientsType>
class AnalyzerForce: public AbstractAnalyzer
{
public:
    
    AnalyzerForce(const IngredientsType& ing_, std::vector<uint32_t> monomers_, uint64_t begCal_);
    
    virtual void initialize();
    virtual bool execute();
    virtual void cleanup();
    
    //functions to get private variable for tests
    const std::vector<uint64_t> getCounterPlus() const { return counterPlus; }
    const std::vector<uint64_t> getCounterMinus() const { return counterMinus; }
    const uint64_t getCounterTries() const { return counterTries; }
  
private:
    
    //holds a reference of the complete system
    const IngredientsType& ingredients;

    IngredientsType forceIngredients;
    
    //! bool for initialize call
    bool isInitialized;

    //! number of mcs after which calculation begins for the first time
    uint64_t beginCalculation;

    //! container for monomer idx to calculate force
    const std::vector<uint32_t> idXSelectedMonomers;

    //counters for jumps up and down for both (all) strands and for every tried jump
    std::vector<uint64_t> counterPlus;
    std::vector<uint64_t> counterMinus;
    uint64_t counterTries;
    
};


/** 
* @brief Constructor handling the new systems parameters and initialize some variables
*
* @param ing_ a reference to the IngredientsType - mainly the system
* @param monomers_ set of monomers to calculate the force
* @param begCal_ age for starting to analyze
*/
template<class IngredientsType>
AnalyzerForce<IngredientsType>::AnalyzerForce(const IngredientsType& ing_, std::vector<uint32_t> monomers_, uint64_t begCal_)
 :ingredients(ing_),beginCalculation(begCal_),isInitialized(false),
 idXSelectedMonomers(monomers_),counterPlus(monomers_.size()),counterMinus(monomers_.size()),
 counterTries(0)
{}


/**
* The initialize function handles the new systems information.
* 
* @details Setup forceIngredients parameters
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template<class IngredientsType>
void AnalyzerForce<IngredientsType>::initialize(){
    if(!isInitialized){
	    std::cout << "AnalyzerForce: initialise" << std::endl;

        //setup forceIngredients without walls and with periodic boundary conditions
	    forceIngredients.modifyMolecules()=ingredients.getMolecules();

	    forceIngredients.setBoxX(ingredients.getBoxX());
	    forceIngredients.setBoxY(ingredients.getBoxY());
	    forceIngredients.setBoxZ(ingredients.getBoxZ());

	    forceIngredients.setPeriodicX(true);
	    forceIngredients.setPeriodicY(true);
	    forceIngredients.setPeriodicZ(true);

	    forceIngredients.modifyBondset().addBFMclassicBondset();

	    forceIngredients.synchronize();

        isInitialized=true;
    }

    execute();
}


/**
* execute()
* 
* @brief calculate number of possible moves
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template<class IngredientsType>
bool AnalyzerForce<IngredientsType>::execute(){
    
	//check if analyzer have been initialized. if not, exit and explain
	if(!isInitialized){
		initialize();
	}
	
	//begin calculation first time at age beginCalculation
	if(ingredients.getMolecules().getAge() >= beginCalculation)
    {
        // copy monomer positions to forceIngredients
        forceIngredients.modifyMolecules()=ingredients.getMolecules();
        forceIngredients.synchronize();
        
        for(uint32_t i=0; i<idXSelectedMonomers.size(); i++){
            //remove constraints on the monomer
            forceIngredients.modifyMolecules()[i].setMovableTag(true);

            MoveLocalSc movePlus;
            movePlus.init(forceIngredients, idXSelectedMonomers.at(i), VectorInt3(0,0,1));
            if(movePlus.check(forceIngredients)){
                counterPlus.at(i)++;
            }

            MoveLocalSc moveMinus;
            moveMinus.init(forceIngredients, idXSelectedMonomers.at(i), VectorInt3(0,0,-1));
            if(moveMinus.check(forceIngredients)){
                counterMinus.at(i)++;
            }
        }
        counterTries++;
    }
}


/**
* cleanup()
* 
* @brief Calculates force and writes output file.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template<class IngredientsType>
void AnalyzerForce<IngredientsType>::cleanup()
{
    // print results into a file
    std::string filename("force.dat");
    
    //construct a list
    std::vector<std::vector<double> > tmpResults(4,std::vector<double>());

    //fill tmpResults
    for(uint64_t i=0; i<idXSelectedMonomers.size(); i++){
        tmpResults[0].push_back(idXSelectedMonomers.at(i));
        tmpResults[1].push_back(counterMinus.at(i));
        tmpResults[2].push_back(counterPlus.at(i));
        tmpResults[3].push_back(log(double(counterMinus.at(i))/double(counterPlus.at(i))));
    }
    //write comments
    std::stringstream comment;
    comment << "# Analyzer force" << std::endl
            << "# idxMonomer\tn-\tn+\tlog(n-/n+)";

    //write file
    ResultFormattingTools::writeResultFile(filename, this->ingredients, tmpResults, comment.str());
}



#endif //ANALYZER_FORCE_H
