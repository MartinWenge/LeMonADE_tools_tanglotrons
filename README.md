# LeMonADE-tools
The abbreviation LeMonADE stands for 
"**L**attice-based **e**xtensible **Mon**te-Carlo **A**lgorithm and **D**evelopment **E**nvironment".
The tools use the LeMonADE library to perform system setups, simulations and analysis of a bond-flutuation-model system [1] [BFM1] [2] [BFM2].

For further details of LeMonADE-project see: [LeMonADE-project]


[BFM1]: http://dx.doi.org/10.1021/ma00187a030  "I. Carmesin, K. Kremer; Macromolecules 21, 2819-2823 (1988)"
 
[BFM2]: http://dx.doi.org/10.1063/1.459901 "H. P. Deutsch, K. Binder; J. Chem. Phys. 94, 2294-2304 (1990)"

[LeMonADE-project]: https://github.com/LeMonADE-project/

## Installation

* Clone the LeMonADE library `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install the LeMonADE library - see [LeMonADE-project] 
* Clone the LeMonADE-tools application `git clone https://github.com/MartinWenge/LeMonADE_tools_tanglotrons.git`
* Install cmake (minimum version 2.6.2)
* You find there directories for analyzers, updaters, features, tests and executables
* To use the tools you only need the executable directory and add the path of the LeMonADE library [LeMonADE-project] to the CMakeLists.txt
* then on Linux use

````sh
    # generates the application in build-directory
    mkdir build
    cd build
    cmake ../
    make
````
The executables can be found in the bin directory.

## Getting Started

````sh
    # yust play around with the programs in the bin directory, for instance
    createFixedChainInSlit -h
````

## Tests
In the test directory unit tests for the analyzers, updaters and features are provided.
Tests are written using [Catch2](https://github.com/catchorg/Catch2).
To run the tests, add your Path to the [LeMonADE-library](https://github.com/LeMonADE-project/) in the CMakeLists.txt and do
````sh
    mkdir build
    cd build
    cmake ../
    make
    ./testTanglotron
````

## Troubleshooting

* if you are facing problems with the program just open an issue


## License

See the LICENSE in the LeMonADE library [LeMonADE-project]
