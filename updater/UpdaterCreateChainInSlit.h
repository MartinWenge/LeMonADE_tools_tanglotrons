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

  //! getter for initialised bool
  const bool getIsInitialized() const { return isInitialized;}
  //! getter for number of executions
  const int32_t getIsExecuted() const { return isExecuted;}
  
private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::addMonomerInsideConnectedPair;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

  // static instance of rng
  RandomNumberGenerators rng;
  
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
UpdaterCreateChainInSlit<IngredientsType>::UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_):
BaseClass(ingredients_), chainLength(chainLength_), slitSize(slitSize_), boxXY(boxXY_), fixType(fixType_),distanceFixpointWall(distanceFixpointWall_),
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
      slitSizedWall.setBase(0,0,slitSize);
      slitSizedWall.setNormal(0,0,1);
      ingredients.addWall(slitSizedWall);
    }
    
    // set periodicity
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(false);

    // set bondset
    ingredients.modifyBondset().addBFMclassicBondset();

    // check parameters
    if ( (3*chainLength) < (slitSize-1) && (fixType != 0) ){
      throw std::runtime_error("UpdaterCreateChainInSlit: chain length is too short for slit size");
    }

    ingredients.synchronize();
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
  ingredients.modifyMolecules()[0].setMovableTag(false);

  // check if setup is completed
  if(ingredients.getMolecules().size() == chainLength){
    ingredients.synchronize();
    isExecuted=true;
    return true;
  }

  // stack with (0,0,2) until wall
  int32_t sizeStack((slitSize-2)/2);
  if(sizeStack > (chainLength-1) ){
    sizeStack = chainLength-1;
  }

  for(uint32_t i=0; i<sizeStack; i++){
    //std::cout << i<<":"<<2*(i+1)<<"("<<slitSize<<")"<<std::endl;
    addMonomerAtPosition(VectorInt3(0,0,2*(i+1)));
    ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
    // last monomer:
    if( (2*(i+1)) == (slitSize-2) && (fixType != 0) ){
      ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setMovableTag(false);
    }
    if( (2*(i+1)) == (slitSize-3) && (fixType != 0) ){
      ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setAllCoordinates(0,0,2*(i+1)+1);
      ingredients.modifyMolecules()[ingredients.getMolecules().size()-1].setMovableTag(false);
    }
  }

  // check if setup is completed
  if(ingredients.getMolecules().size() == chainLength){
    ingredients.synchronize();
    isExecuted=true;
    return true;
  }
  
  // add remaining monomers
  ingredients.synchronize();
  int32_t numOfRemainingMonos(chainLength-ingredients.getMolecules().size());
  for(int32_t i=0; i<numOfRemainingMonos; i++){
    int32_t mono( ingredients.getMolecules().size()/2 );
    
    if(mono>0){
      addMonomerInsideConnectedPair(mono,ingredients.getMolecules().getNeighborIdx(mono,0)) ? : throw std::runtime_error("UpdaterCreateChainInSlit: addMonomerInsideConnectedPair is not able to place a monomer!");
    }else{
      addMonomerInsideConnectedPair(mono,ingredients.getMolecules().getNeighborIdx(mono,1)) ? : throw std::runtime_error("UpdaterCreateChainInSlit: addMonomerInsideConnectedPair is not able to place a monomer!");
    }
  }
  
  // reorder monomer ids and synchronize
  linearizeSystem();
  ingredients.synchronize();

  isExecuted=true;
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
inline int32_t UpdaterCreateChainInSlit<IngredientsType>::pow2roundup(int32_t x){
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