// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
// ## The main file
// This is the main file for the biomineralization example. Here we can see the programme sequence and how the system is solved using newton's method.

// ### Includes
#include <config.h>

// Standard header file for C++, to get time and date information.
#include <ctime>

// Standard header file for C++, for in- and output.
#include <iostream>

// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers. So we need some includes from that.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

// In Dumux a property system is used to specify the model. For this, different properties are defined containing type definitions, values and methods. All properties are declared in the file properties.hh.
#include <dumux/common/properties.hh>
// The following file contains the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
#include <dumux/common/parameters.hh>
// The file dumuxmessage.hh contains the class defining the start and end message of the simulation.
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
// We include the linear solver to be used to solve the linear system
#include <dumux/linear/amgbackend.hh>
// We include the nonlinear newtons method
#include <dumux/nonlinear/newtonsolver.hh>
// Further we include assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation)
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

// We need the following class to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>

// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

// We include the problem file which defines initial and boundary conditions to describe our example problem
#include "problem.hh"

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-ParameterFile Parameter file (Input file) \n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

// ### Beginning of the main function
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv, usage);

    // ### Create the grid
    // A gridmanager tries to create the grid either from a grid file or the input file.
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // ### Setup and solving of the problem
    // #### Setup
    // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
    // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // In the problem, we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // We initialize the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
        const auto fileName = getParam<std::string>("Restart.File");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem, SolidSystem>();
        loadSolution(x, fileName, pvName, *gridGeometry);
    }
    else
        problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We intialize the vtk output module. Each model has a predefined model specific output with relevant parameters for that model.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    //we add permeability as a specific output
    vtkWriter.addField(problem->getPermeability(), "Permeability");
    //We update the output fields before write
    problem->updateVtkOutput(x);
    vtkWriter.write(restartTime);

    //We instantiate the time loop and define an episode index to keep track of the the actual episode number
    int episodeIdx_ = 0;
    //We  set the initial episodeIdx in the problem to zero
    problem->setEpisodeIdx(episodeIdx_);
    //We  set the time loop with episodes (check points)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //We use a time loop with episodes if an injection parameter file is specified in the input
    if (hasParam("Injection.InjectionParamFile"))
    {
        //We first read the times of the checkpoints from the injection file
        std::vector<Scalar> epiEnd_;
        std::string injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");
        std::ifstream injectionData;
        std::string row;
        //We open the injection data file
        injectionData.open( injectionParameters_);
        //We check if we could open the file
        if (not injectionData.is_open())
        {
            std::cerr << "\n\t -> Could not open file '"
                    << injectionParameters_
                    << "'. <- \n\n\n\n";
            exit(1) ;
        }
        Scalar tempTime = 0;

        //We read the episode end times data from the injection data file
        while (!injectionData.eof())
        {
            getline(injectionData, row);
            if (row == "EndTimes")
            {
                getline(injectionData, row);
                while (row != "#")
                {
                    if (row != "#")
                    {
                        std::istringstream ist(row);
                        ist >> tempTime;
                        //The injection data file contains the time in hours, we convert that to seconds!
                        epiEnd_.push_back(tempTime*3600);
                    }
                    getline(injectionData, row);
                }
            }
        }
        injectionData.close();
        //We check the injection data against the number of injections specified in the parameter file and print an error message it this check fails
        int numInjections_ = getParam<Scalar>("Injection.numInjections");
        if (epiEnd_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection end times specified in the injection data file do not match!"
                        <<"\n numInjections from parameter file = "<<numInjections_
                        <<"\n numEpisodes from injection data file = "<<epiEnd_.size()
                        <<"\n Abort!\n";
            exit(1) ;
        }
        //We set the time loop with episodes of various lengths as specified in the injection file
        // set the episode ends /check points:
        for (int epiIdx = 0; epiIdx < epiEnd_.size(); ++epiIdx)
        {
            timeLoop->setCheckPoint(epiEnd_[epiIdx]);
        }
        //We also set the initial episodeIdx in the problem to zero, as in the problem the boundary conditions depend on the episode.
        problem->setEpisodeIdx(episodeIdx_);
    }

    //If nothing specified, we use no episodes
    else
    {
        //In this case, we set the time loop with one big episodes
        timeLoop->setCheckPoint(tEnd);
    }

    //We set the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop);

    //We set the linear solver
    using LinearSolver = Dumux::ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    //We set the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // We start the time loop.
    timeLoop->start(); do
    {
        //We set the time and time step size to be used in the problem
        problem->setTime( timeLoop->time() + timeLoop->timeStepSize() );
        problem->setTimeStepSize( timeLoop->timeStepSize() );

        //We set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        //We solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        //We make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        //We advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        //We update the output fields before write
        problem->updateVtkOutput(x);

        //We write vtk output
        vtkWriter.write(timeLoop->time());

        //We report statistics of this time step
        timeLoop->reportTimeStep();

        //If using episodes/check points, we count episodes and update the episode indices in the main file and the problem.
        if (hasParam("Injection.InjectionParamFile") || hasParam("TimeLoop.EpisodeLength"))
        {
            if(timeLoop->isCheckPoint() && hasParam("Injection.InjectionParamFile"))
            {

                episodeIdx_++;
                problem->setEpisodeIdx(episodeIdx_);
            }
        }

        //We set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());
    timeLoop->finalize(leafGridView.comm());

    // ### Final Output
    // We print dumux end message.
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
