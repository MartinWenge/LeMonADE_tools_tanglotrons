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

// use the catch file but do not add the #define CATCH_CONFIG_MAIN !!
#include "catch.hpp"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureWall.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>

#include "UpdaterCreateChainInSlit.h"
#include "AnalyzerForce.h"

typedef LOKI_TYPELIST_5(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <bool> >, FeatureWall, FeatureAttributes, FeatureFixedMonomers) Features;
typedef ConfigureSystem<VectorInt3,Features,4> Config;
typedef Ingredients<Config> IngredientsType;

TEST_CASE( "TestAnalyzerForce_constructor" ) {
    IngredientsType ingredients;
    // AnalyzerForce(const IngredientsType& ing_, std::vector<int32_t> monomers_, uint64_t begCal_);
    AnalyzerForce<IngredientsType> Anna(ingredients, std::vector<uint32_t> (4,0), 1234);
    CHECK(Anna.getCounterPlus().size() == 4);
    CHECK(Anna.getCounterMinus().size() == 4);
    CHECK(Anna.getCounterTries() == 0);
}

TEST_CASE( "TestAnalyzerForce_execute" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 1, 5, 16, UpdaterCreateChainInSlit<IngredientsType>::SINGLE_FIXPOINT_BOTTOM, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 1);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==8);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,5));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    CHECK(ingredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(ingredients.getMolecules()[0].getMovableTag() == false);

    // start using the analyzer
    AnalyzerForce<IngredientsType> Anna(ingredients, std::vector<uint32_t> (1,0), 0);
    CHECK(Anna.getCounterPlus().size() == 1);
    CHECK(Anna.getCounterMinus().size() == 1);
    CHECK(Anna.getCounterTries() == 0);
    Anna.initialize();
    CHECK(Anna.getCounterPlus().at(0) == 1);
    CHECK(Anna.getCounterMinus().at(0) == 1);
    CHECK(Anna.getCounterTries() == 1);
    Anna.execute();
    CHECK(Anna.getCounterPlus().at(0) == 2);
    CHECK(Anna.getCounterMinus().at(0) == 2);
    CHECK(Anna.getCounterTries() == 2);

    // set wall in front of monomer
    Wall blocker;
    blocker.setBase(0,0,2);
    blocker.setNormal(0,0,1);
    ingredients.addWall(blocker);

    // wall should be ignored
    AnalyzerForce<IngredientsType> Benno(ingredients, std::vector<uint32_t> (1,0), 0);
    CHECK(Benno.getCounterPlus().size() == 1);
    CHECK(Benno.getCounterMinus().size() == 1);
    CHECK(Benno.getCounterTries() == 0);
    Benno.initialize();
    CHECK(Benno.getCounterPlus().at(0) == 1);
    CHECK(Benno.getCounterMinus().at(0) == 1);
    CHECK(Benno.getCounterTries() == 1);

    //set monomer in front of grafted monomer
    ingredients.modifyMolecules().addMonomer(0,0,2);
    ingredients.modifyMolecules().connect(0,1);

    // monomer should block plus movement
    AnalyzerForce<IngredientsType> Calvin(ingredients, std::vector<uint32_t> (1,0), 0);
    CHECK(Calvin.getCounterPlus().size() == 1);
    CHECK(Calvin.getCounterMinus().size() == 1);
    CHECK(Calvin.getCounterTries() == 0);
    Calvin.initialize();
    CHECK(Calvin.getCounterPlus().at(0) == 0);
    CHECK(Calvin.getCounterMinus().at(0) == 1);
    CHECK(Calvin.getCounterTries() == 1);

    // move second monomer to block minus movement
    ingredients.modifyMolecules()[1].setAllCoordinates(0,0,3);

    AnalyzerForce<IngredientsType> David(ingredients, std::vector<uint32_t> (1,0), 0);
    CHECK(David.getCounterPlus().size() == 1);
    CHECK(David.getCounterMinus().size() == 1);
    CHECK(David.getCounterTries() == 0);
    David.initialize();
    CHECK(David.getCounterPlus().at(0) == 1);
    CHECK(David.getCounterMinus().at(0) == 0);
    CHECK(David.getCounterTries() == 1);

    // test the same with "old" analyzer Calvin
    Calvin.execute();
    CHECK(Calvin.getCounterPlus().at(0) == 1);
    CHECK(Calvin.getCounterMinus().at(0) == 1);
    CHECK(Calvin.getCounterTries() == 2);

    // try both monomers
    AnalyzerForce<IngredientsType> Emil(ingredients, std::vector<uint32_t> {0,1}, 0);
    CHECK(Emil.getCounterPlus().size() == 2);
    CHECK(Emil.getCounterMinus().size() == 2);
    CHECK(Emil.getCounterTries() == 0);
    Emil.initialize();
    CHECK(Emil.getCounterPlus().at(0) == 1);
    CHECK(Emil.getCounterPlus().at(1) == 0);
    CHECK(Emil.getCounterMinus().at(0) == 0);
    CHECK(Emil.getCounterMinus().at(1) == 1);
    CHECK(Emil.getCounterTries() == 1);

    // "move" the system
    ingredients.modifyMolecules()[1].setAllCoordinates(0,0,2);

    Emil.execute();
    CHECK(Emil.getCounterPlus().at(0) == 1);
    CHECK(Emil.getCounterPlus().at(1) == 1);
    CHECK(Emil.getCounterMinus().at(0) == 1);
    CHECK(Emil.getCounterMinus().at(1) == 1);
    CHECK(Emil.getCounterTries() == 2);

    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    Emil.execute();
    CHECK(Emil.getCounterPlus().at(0) == 2);
    CHECK(Emil.getCounterPlus().at(1) == 2);
    CHECK(Emil.getCounterMinus().at(0) == 2);
    CHECK(Emil.getCounterMinus().at(1) == 2);
    CHECK(Emil.getCounterTries() == 3);

    // test a small chain out of the box
    IngredientsType largeIngredients;

    //UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    UpdaterCreateChainInSlit<IngredientsType> Secundus(largeIngredients, 6, 16, 16, UpdaterCreateChainInSlit<IngredientsType>::FIXED_AT_WALL_AND_IN_SPACE, 12);
    Secundus.initialize();
    CHECK(largeIngredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(largeIngredients.getMolecules()[1]==VectorInt3(0,0,2));
    CHECK(largeIngredients.getMolecules()[4]==VectorInt3(0,0,8));
    CHECK(largeIngredients.getMolecules()[5]==VectorInt3(0,0,10));

    CHECK(largeIngredients.getMolecules()[0].getMovableTag()==false);
    CHECK(largeIngredients.getMolecules()[1].getMovableTag()==true);
    CHECK(largeIngredients.getMolecules()[4].getMovableTag()==true);
    CHECK(largeIngredients.getMolecules()[5].getMovableTag()==false);

    AnalyzerForce<IngredientsType> Fabian(largeIngredients, std::vector<uint32_t> {0,5}, 0);
    REQUIRE(Fabian.getCounterPlus().size() == 2);
    REQUIRE(Fabian.getCounterMinus().size() == 2);
    CHECK(Fabian.getCounterTries() == 0);
    Fabian.initialize();
    CHECK(Fabian.getCounterPlus().at(0) == 0);
    CHECK(Fabian.getCounterMinus().at(0) == 1);
    CHECK(Fabian.getCounterPlus().at(1) == 1);
    CHECK(Fabian.getCounterMinus().at(1) == 0);
    CHECK(Fabian.getCounterTries() == 1);
    Fabian.execute();
    CHECK(Fabian.getCounterPlus().at(0) == 0);
    CHECK(Fabian.getCounterMinus().at(0) == 2);
    CHECK(Fabian.getCounterPlus().at(1) == 2);
    CHECK(Fabian.getCounterMinus().at(1) == 0);
    CHECK(Fabian.getCounterTries() == 2);

}

