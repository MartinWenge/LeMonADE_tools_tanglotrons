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

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureWall.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include "UpdaterCreateChainInSlit.h"


// read in command line options
#include <boost/program_options.hpp>
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  /* read arguments
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  */
  std::string filename;
  uint32_t LinearChainLength, slitSize, box, mode, fixedPosition;

  
  try{
    options_description desc{"Create BFM-CoDendrimers and linear chains in simulation box\nAllowed options"};
    desc.add_options()
      ("help,h", "produce help message")
      ("filename,f", value<std::string>(&filename)->default_value("config.bfm"), "filename")
      ("chainlength,n", value<uint32_t>(&LinearChainLength)->default_value(1), "linear chain length")
      ("box,b", value<uint32_t>(&box)->default_value(128), "boxsize ( in x,y)")
      ("slit,s", value<uint32_t>(&slitSize)->default_value(0), "size of slit (in z)")
      ("positionZ,p", value<uint32_t>(&fixedPosition)->default_value(0), "fixed monomer position")
      ("mode,m", value<uint32_t>(&mode)->default_value(0), "mode: 0=grafted chain, 1=chain fixed between walls, 2=monomer fixed in space");
      
    variables_map options_map;
    store(parse_command_line(argc, argv, desc), options_map);
    notify(options_map); 
    
    // help option
    if (options_map.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }
  } catch (const error &ex){
    std::cerr << ex.what() << '\n';
  }
  
  // display all internal variables:
  std::cout << "file name = '" << filename <<"'"<<std::endl
  << "LinearChainLength = '" << LinearChainLength <<"'\t"
  << "mode = '" << mode <<"'"<<std::endl
  << "box size = '" << box <<"' ("<<slitSize<<")"<<std::endl;
  
  /* initialize system
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
  */

  typedef LOKI_TYPELIST_5(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <bool> >, FeatureWall, FeatureAttributes, FeatureFixedMonomers) Features;
  const uint max_bonds=4;
  // define maximal number of bonds
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  //defines ingredients type (3D Vektor, oben genannte features, maximale anzahl Bindungspartner)
  typedef Ingredients<Config> IngredientsType;
  //defines the ingedients type with the configuration above
  IngredientsType ingredients;
  
  /* set up random number generator (static object)
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
  RandomNumberGenerators rng;
  rng.seedAll();
  
  /* use TaskManager
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
  TaskManager taskManager;
   // UpdaterCreateChainInSlit(IngredientsType& ingredients_, uint32_t chainLength_, uint32_t slitSize_, uint32_t boxXY_, int fixType_, uint32_t distanceFixpointWall_=0);
    // SINGLE_FIXPOINT_BOTTOM=0,    DOUBLE_FIXED_AT_WALLS=1,    FIXED_AT_WALL_AND_IN_SPACE=2
  taskManager.addUpdater(new UpdaterCreateChainInSlit<IngredientsType>(ingredients, LinearChainLength, slitSize, box, mode, fixedPosition));

  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(filename,ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE ));

  taskManager.initialize();
  taskManager.run(1);
  taskManager.cleanup();
  /* */
  return true;
}
