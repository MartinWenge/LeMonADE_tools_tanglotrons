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

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include "AnalyzerForce.h"

// read in command line options
#include <boost/program_options.hpp>
using namespace boost::program_options;

int main(int argc, char* argv[])
{
  /* read arguments
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  */
  std::string ifilename,ofilename;
  int32_t max_mcs, save_interval, force_interval, relaxtime;
  std::vector<uint32_t> selectedMonomers;

  try{
    options_description desc{"Set up all paramters for SimulatorSlitChain with force measurement\nnummcs, nforce and nsave are requested to give useful values when dividing by each other\nAllowed options"};
    desc.add_options()
      ("help,h", "produce help message")
      ("ifilename,i", value<std::string>(&ifilename)->default_value("config.bfm"), "input filename")
      ("ofilename,o", value<std::string>(&ofilename)->default_value("configRun.bfm"), "output filename")
      ("nummcs,n", value<int32_t>(&max_mcs)->default_value(10000), "number of MCS")
      ("numsave,s", value<int32_t>(&save_interval)->default_value(1000), "mcs intervall to save the current config")
      ("nforce,f", value<int32_t>(&force_interval)->default_value(10), "mcs intervall to analyzer force")
      ("selection,v", value<vector<uint32_t> >(&selectedMonomers)->multitoken(), "vector of monomers to measure force {a,b,...}")
      ("relax,r", value<int32_t>(&relaxtime)->default_value(10), "num mcs before starting force calculation");
      
    variables_map options_map;
    store(parse_command_line(argc, argv, desc), options_map);
    notify(options_map); 
    
    // help option
    if (options_map.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    // show paramters
    
    for (const auto& it : options_map) {
      std::cout << it.first.c_str() << " ";
      auto& value = it.second.value();
      if (auto v = boost::any_cast<int32_t>(&value)){
	    std::cout << *v << std::endl;
      }else if (auto v = boost::any_cast<std::string>(&value)){
	    std::cout << *v << std::endl;
      }else if (auto v = boost::any_cast<std::vector<uint32_t> >(&value)){
          for( const auto& a:*v)
	         std::cout <<a << ", ";
          std::cout << std::endl;
      }else{
	    std::cout << "error"<< std::endl;
      }
    }
    

  }catch (const error &ex){
    std::cerr << ex.what() << '\n';
  }


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

  try{
    // prepare cycles
    int simulatorCycles(0), simulatorInterval(0), writePeriod(0);

    if(save_interval > force_interval){
        writePeriod=(save_interval/force_interval);
        simulatorCycles=(max_mcs/force_interval);
        simulatorInterval=force_interval;
    }else{
        throw std::runtime_error("force_intervall is smaller than save_interval");
    }

    TaskManager taskmanager;
    taskmanager.addUpdater(new UpdaterReadBfmFile<IngredientsType>(ifilename,ingredients,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE),0);
    taskmanager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,simulatorInterval));

    taskmanager.addAnalyzer(new AnalyzerForce<IngredientsType>(ingredients,selectedMonomers,relaxtime/force_interval));

    ofilename=(ofilename.substr(0,ofilename.find_last_of(".")));
    taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(ofilename+".bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND ),writePeriod);
    taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(ofilename+"_lastconfig.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::OVERWRITE ),writePeriod);

    taskmanager.initialize();
    taskmanager.run(simulatorCycles);
    taskmanager.cleanup();

  }catch(std::exception& err){
    std::cerr<<err.what();
  }
  
  return 0;

}
