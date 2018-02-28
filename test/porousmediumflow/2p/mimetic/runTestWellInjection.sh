#!/bin/sh

make test_well2p

######################
###### Mimetic #######
inputArgs="-Newton.UseLineSearch 0 -Newton.EnableChop 0 -Grid.Refinement 1 -LinearSolver.Verbosity 1 -LinearSolver.MaxIterations 2000"
outputName="wellinjection_hybrid_ilu0bicgstab"


echo $outputName
rm $outputName".out"
./test_well2p -ParameterFile ./test_well2p.input -Problem.Name $outputName $inputArgs > $outputName".out"


######################
###### MPFA-O ########
make test_well2p_mpfao

inputArgs="-Newton.UseLineSearch 0 -Newton.EnableChop 0 -Grid.Refinement 1 -LinearSolver.Verbosity 1 -LinearSolver.MaxIterations 2000"
outputName="wellinjection_mpfao_ilu0bicgstab"


echo $outputName
rm $outputName".out"
./test_well2p_mpfao -ParameterFile ./test_well2p_mpfao.input -Problem.Name $outputName $inputArgs > $outputName".out"