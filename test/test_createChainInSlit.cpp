// use the catch file but do not add the #define CATCH_CONFIG_MAIN !!
#include "catch.hpp"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>

//#include "pdaterCreateChainInSlit.h"

typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <bool> >) Features;
typedef ConfigureSystem<VectorInt3,Features,4> Config;
typedef Ingredients<Config> IngredientsType;

TEST_CASE( "UpdaterCreateChainInSlit_freeChain" ) {
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    IngredientsType ingredients;
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(false);

    ingredients.modifyMolecules().addMonomer(0,0,0); // 0

    // REQUIRE(...); -> check if true, do not continue if false
    REQUIRE(ingredients.getMolecules().size() == 1);

    // CHECK(...); -> check if true, continue if false
    CHECK(ingredients.isPeriodicX());

}
