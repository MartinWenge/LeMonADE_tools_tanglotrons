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

#ifndef LEMONADE_UPDATERCREATE_CHAININSLIT
#define LEMONADE_UPDATERCREATE_CHAININSLIT
/**
 * @file
 *
 * @class UpdaterCreateChainInSlit
 *
 * @brief Updater setting up a simple system containing a linear chain in a slit with different monomer fixes
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>

template<class IngredientsType>
class UpdaterCreateChainInSlit: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
public:
  UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);

  enum FIX_TYPE{
    SINGLE_FIXPOINT_BOTTOM=0,
    DOUBLE_FIXED_AT_WALLS=1,
    FIXED_AT_WALL_AND_IN_SPACE=2
  };
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::addMonomerInsideConnectedPair;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;
  
  //! number of monomers in a chain
  uint32_t chainLength;
  
  //! size of box in z direction
  uint32_t slitSize;
  
  //! size of box in xy direction
  uint32_t boxXY;

  //! setup type using FIX_TYPE
  int fixType;

  //! in case of FIXED_AT_WALL_AND_IN_SPACE: distance of fixed monomer in space to the wall
  uint32_t distanceFixpointWall;
  
  //! bool for execution
  bool isInitialized;

  //! bool for execution
  bool isExecuted;

  //! helper function:
  int32_t pow2roundup(int32_t a);

};

/** 
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param chainLength_ number of bfm units ("monomers") in the chain
* @param slitSize_ boxsize in z direction
* @param boxXY_ boxsize in xy direction
* @param fixType_ type of system setup using FIX_TYPE
* @param distanceFixpointWall_ distance between fixpoint of chain end (FIXED_AT_WALL_AND_IN_SPACE) and wall
*/
template < class IngredientsType >
UpdaterCreateChainInSlit<IngredientsType>::UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0):
BaseClass(ingredients_), chainLength(chainLength_), slitSize(slitSize_), boxXY(boxXY_), fixType(fixType_),distanceFixpointWall(distanceFixpointWall_)
isInitialized(false), isExecuted(false)
{}

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateChainInSlit<IngredientsType>::initialize(){

  if(!isInitialized){
    std::cout << "initialize UpdaterCreateChainInSlit" << std::endl;

    // setup box
    ingredients.setBoxX(boxXY);
    ingredients.setBoxY(boxXY);
    // adjust the z boxsize to next power of two
    ingredients.setBoxZ(pow2roundup(slitSize));

    // add the flexible wall
    if(ingredients.getBoxZ()!=slitSize){
      Wall slitSizedWall;
      wallUp.setBase(0,0,slitSize);
      wallUp.setNormal(0,0,1);
      ingredients.addWall(slitSizedWall);
    }
    
    // set periodicity
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(false);

    // set bondset
    ingredients.modifyBondset().addBFMclassicBondset();

    // check parameters
    if ( (3*chainLength) < (slitSize-1) ){
      throw std::runtime_error("UpdaterCreateChainInSlit: chain length is too short for slit size");
    }

    isInitialized=true;
  }
  
  execute();
}

/**
* Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterCreateChainInSlit<IngredientsType>::execute(){
  if(isExecuted)
    return true;
  
  // start with first monomer at (0,0,0)
  addMonomerAtPosition(VectorInt3(0,0,0));

  // stack with (0,0,3) until wall
  int32_t sizeStack(slitSize/3);
  for(uint32_t i=0; i<sizeStack; i++){
    addMonomerAtPosition(VectorInt3(0,0,3*i));
    ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  }
  int32_t lastStep((slitSize-1)%3);
  switch (lastStep){
    case 0: break;
    case 1: {
              ingredients.modifyMolecules().[ingredients.getMolecules().size()-1] + VectorInt3(0,0,-1);
            }
    case 2: {
              addMonomerAtPosition(VectorInt3(0,0,3*slitSize+2));
              ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
              break;
            }
  }

  // add remaining monomers
  int32_t numOfRemainingMonos(chainLength-ingredients.getMolecules().size());
  for(int32_t i=0; i<numOfRemainingMonos; i++){
    int32_t mono( rng.r250_rand32() % (ingredients.getMolecules().size()-2)+1 );
    addMonomerInsideConnectedPair(mono,ingredients.getMolecules().getNeighborIndex(mono,1));
  }
  
  linearizeSystem();
}

/**
* Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateChainInSlit<IngredientsType>::cleanup(){
  
}

/**
* Helper function to use feature Lattice power of two
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
inline void UpdaterCreateChainInSlit<IngredientsType>::pow2roundup(int x){
    if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

#endif /* LEMONADE_UPDATERCREATE_CHAININSLIT */