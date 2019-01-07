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

typedef LOKI_TYPELIST_5(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <bool> >, FeatureWall, FeatureAttributes, FeatureFixedMonomers) Features;
typedef ConfigureSystem<VectorInt3,Features,4> Config;
typedef Ingredients<Config> IngredientsType;

TEST_CASE( "UpdaterCreateChainInSlit_bottomFixedChain_oneMonomer" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 1, 14, 16, UpdaterCreateChainInSlit<IngredientsType>::SINGLE_FIXPOINT_BOTTOM, 0);

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
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,14));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    CHECK(ingredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(ingredients.getMolecules()[0].getMovableTag() == false);

}

TEST_CASE( "UpdaterCreateChainInSlit_bottomFixedChain_twoMonomer" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 2, 14, 16, UpdaterCreateChainInSlit<IngredientsType>::SINGLE_FIXPOINT_BOTTOM, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 2);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,14));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    CHECK(ingredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(ingredients.getMolecules()[1]==VectorInt3(0,0,2));
    CHECK(ingredients.getMolecules().areConnected(0,1));
    CHECK(ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK(ingredients.getMolecules()[1].getMovableTag() == true);

}

TEST_CASE( "UpdaterCreateChainInSlit_bottomFixedChain_chainShort" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 5, 8, 16, UpdaterCreateChainInSlit<IngredientsType>::SINGLE_FIXPOINT_BOTTOM, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 5);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==8);

    CHECK(ingredients.getWalls().size()==0);

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
    }
    CHECK(ingredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK(ingredients.getMolecules()[1].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[2].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[3].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[4].getMovableTag() == true);

}

TEST_CASE( "UpdaterCreateChainInSlit_bottomFixedChain_chainLong" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 23, 18, 16, UpdaterCreateChainInSlit<IngredientsType>::SINGLE_FIXPOINT_BOTTOM, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 23);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==32);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,18));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
    }

    CHECK(ingredients.getMolecules()[0]==VectorInt3(0,0,0));
    CHECK(ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK(ingredients.getMolecules()[1].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[2].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[3].getMovableTag() == true);
    CHECK(ingredients.getMolecules()[4].getMovableTag() == true);
}

TEST_CASE( "UpdaterCreateChainInSlit_doubleFixedStrand_chainShort" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 4, 9, 16, UpdaterCreateChainInSlit<IngredientsType>::DOUBLE_FIXED_AT_WALLS, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 4);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,9));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
    }
    CHECK( ingredients.getMolecules()[0] == VectorInt3(0,0,0) );
    CHECK( ingredients.getMolecules()[1] == VectorInt3(0,0,2) );
    CHECK( ingredients.getMolecules()[2] == VectorInt3(0,0,4) );
    CHECK( ingredients.getMolecules()[3] == VectorInt3(0,0,7) );
    CHECK( ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK( ingredients.getMolecules()[1].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[2].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[3].getMovableTag() == false);

}

TEST_CASE( "UpdaterCreateChainInSlit_doubleFixedStrand_chainLong" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 12, 9, 16, UpdaterCreateChainInSlit<IngredientsType>::DOUBLE_FIXED_AT_WALLS, 0);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 12);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==1);
    CHECK(ingredients.getWalls().at(0).getBase() == VectorInt3(0,0,9));
    CHECK(ingredients.getWalls().at(0).getNormal() == VectorInt3(0,0,1));

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
        CHECK(ingredients.getMolecules()[i].getZ() < 8 );
    }
    CHECK( ingredients.getMolecules()[0] == VectorInt3(0,0,0) );
    CHECK( ingredients.getMolecules()[11] == VectorInt3(0,0,7) );
    CHECK( ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK( ingredients.getMolecules()[1].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[2].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[3].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[4].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[5].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[6].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[7].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[8].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[9].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[10].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[11].getMovableTag() == false);
}

TEST_CASE( "UpdaterCreateChainInSlit_inSpaceFixedChain_chainShort" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 4, 16, 16, UpdaterCreateChainInSlit<IngredientsType>::FIXED_AT_WALL_AND_IN_SPACE, 9);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 4);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==0);

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
    }
    CHECK( ingredients.getMolecules()[0] == VectorInt3(0,0,0) );
    CHECK( ingredients.getMolecules()[1] == VectorInt3(0,0,2) );
    CHECK( ingredients.getMolecules()[2] == VectorInt3(0,0,4) );
    CHECK( ingredients.getMolecules()[3] == VectorInt3(0,0,7) );
    CHECK( ingredients.getMolecules()[0].getMovableTag() == false);
    CHECK( ingredients.getMolecules()[1].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[2].getMovableTag() == true);
    CHECK( ingredients.getMolecules()[3].getMovableTag() == false);

}


TEST_CASE( "UpdaterCreateChainInSlit_inSpaceFixedChain_chainLong" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;

    // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
    UpdaterCreateChainInSlit<IngredientsType> Primus(ingredients, 20, 16, 16, UpdaterCreateChainInSlit<IngredientsType>::FIXED_AT_WALL_AND_IN_SPACE, 9);

    CHECK(Primus.getIsInitialized()==false);
    CHECK(Primus.getIsExecuted()==false);

    Primus.initialize();
    CHECK(Primus.getIsInitialized()==true);
    CHECK(Primus.getIsExecuted()==true);

    CHECK(ingredients.getMolecules().size() == 20);
    CHECK(ingredients.isPeriodicX());
    CHECK(ingredients.isPeriodicY());
    CHECK(!ingredients.isPeriodicZ());

    CHECK(ingredients.getBoxX()==16);
    CHECK(ingredients.getBoxY()==16);
    CHECK(ingredients.getBoxZ()==16);

    CHECK(ingredients.getWalls().size()==0);

    for(uint32_t i=0;i<ingredients.getMolecules().size()-1;i++){
        CHECK(ingredients.getMolecules().areConnected(i, i+1));
        if(i==0 || i== 19){
            CHECK( ingredients.getMolecules()[i].getMovableTag() == false);
        }else{
            CHECK( ingredients.getMolecules()[i].getMovableTag() == true);
        }
    }

    CHECK( ingredients.getMolecules()[0] == VectorInt3(0,0,0) );
    CHECK( ingredients.getMolecules()[19] == VectorInt3(0,0,7) );
    
}
