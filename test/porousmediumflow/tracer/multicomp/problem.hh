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
/**
 * \file
 * \ingroup TracerTest
 * \brief Definition of a problem for the MaxwellStefan problem:
 * A rotating velocity field mixes a MaxwellStefan band in a porous groundwater reservoir.
 */

#ifndef DUMUX_MAXWELL_STEFAN_TEST_PROBLEM_HH
#define DUMUX_MAXWELL_STEFAN_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"
#include <dumux/flux/maxwellstefanslaw.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux {
template <class TypeTag>
class MaxwellStefanTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct MaxwellStefanTest { using InheritsFrom = std::tuple<Tracer>; };
struct MaxwellStefanTestCC { using InheritsFrom = std::tuple<MaxwellStefanTest, CCTpfaModel>; };
struct MaxwellStefanTestBox { using InheritsFrom = std::tuple<MaxwellStefanTest, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanTest> { using type = MaxwellStefanTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MaxwellStefanTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = MaxwellStefanTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::MaxwellStefanTest> { static constexpr bool value = true; };

//! Here we set FicksLaw or MaxwellStefansLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::MaxwellStefanTest> { using type = MaxwellStefansLaw<TypeTag>; };

//! A simple fluid system with one MaxwellStefan component
template<class TypeTag>
class MaxwellStefanTracerFluidSystem
: public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, MaxwellStefanTracerFluidSystem<TypeTag>>

{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    static constexpr bool isTracerFluidSystem()
    { return true; }
    //! The number of components
    static constexpr int numComponents = 3;

    static constexpr int compOneIdx = 0;//first major component
    static constexpr int compTwoIdx = 1;//second major component
    static constexpr int compThreeIdx = 2;//secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
        case compOneIdx: return "CompOne";
        case compTwoIdx: return "CompTwo";
        case compThreeIdx:return "CompThree";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Gas"; }

    //! Molar mass in kg/mol of the component with index compIdx.
    static Scalar molarMass(unsigned int compIdx)
    { return 0.02896; /*air*/ }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
      if (compIdx == compOneIdx)
          return 0;
      if (compIdx == compTwoIdx)
          return 83.3e-6;
      if (compIdx == compThreeIdx)
          return 68.0e-6;
       DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of component "
                       << compIdx <<" is undefined!\n");
    }

    //! Binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIIdx,
                                             unsigned int compJIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        if (compIIdx == compOneIdx && compJIdx == compTwoIdx)
            return 83.3e-6;
        if (compIIdx == compOneIdx && compJIdx == compThreeIdx)
            return 68.0e-6;
        if (compIIdx == compTwoIdx && compJIdx == compThreeIdx)
            return 16.8e-6;
        DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx << " is undefined!\n");
    }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component \f$M_\kappa\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\kappa} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        return density(fluidState, phaseIdx)/molarMass(0);
    }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MaxwellStefanTest> { using type = MaxwellStefanTracerFluidSystem<TypeTag>; };

} // end namespace Properties

/*!
 * \ingroup TracerTest
 * \brief Definition of a problem for the MaxwellStefan problem.
 *
 * This problem uses the MaxwellStefan equations.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxmaxwellstefan -ParameterFile ./test_boxmaxwellstefan.input</tt> or
 * <tt>./test_ccmaxwellstefan -ParameterFile ./test_ccMaxwellstefan.input</tt>
 */
template <class TypeTag>
class MaxwellStefanTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    MaxwellStefanTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';

        plotOutput_ = false;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    //! Called after every time step
    void plotComponentsOverTime(const SolutionVector& curSol, Scalar time)
    {
        if (plotOutput_)
        {
            Scalar x_CompThree_left = 0.0;
            Scalar x_CompTwo_left = 0.0;
            Scalar x_CompThree_right = 0.0;
            Scalar x_CompTwo_right = 0.0;
            Scalar x_CompOne_left = 0.0;
            Scalar x_CompOne_right = 0.0;
            Scalar i = 0.0;
            Scalar j = 0.0;
            if (!(time < 0.0))
            {
                for (const auto& element : elements(this->gridGeometry().gridView()))
                {
                    auto fvGeometry = localView(this->gridGeometry());
                    fvGeometry.bindElement(element);

                    const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto& globalPos = scv.dofPosition();
                        VolumeVariables volVars;
                        volVars.update(elemSol, *this, element, scv);

                        if (globalPos[0] < 0.5)
                        {
                            x_CompThree_left += volVars.moleFraction(0,2);

                            x_CompTwo_left += volVars.moleFraction(0,1);
                            x_CompOne_left += volVars.moleFraction(0,0);
                            i +=1;
                        }
                        else
                        {
                            x_CompThree_right += volVars.moleFraction(0,2);
                            x_CompTwo_right += volVars.moleFraction(0,1);
                            x_CompOne_right += volVars.moleFraction(0,0);
                            j +=1;
                        }

                    }
                }
                x_CompThree_left /= i;
                x_CompTwo_left /= i;
                x_CompOne_left /= i;
                x_CompThree_right /= j;
                x_CompTwo_right /= j;
                x_CompOne_right /= j;

                //do a gnuplot
                x_.push_back(time); // in seconds
                y_.push_back(x_CompTwo_left);
                y2_.push_back(x_CompTwo_right);
                y3_.push_back(x_CompThree_left);
                y4_.push_back(x_CompThree_right);
                y5_.push_back(x_CompOne_left);
                y6_.push_back(x_CompOne_right);

                gnuplot_.resetPlot();
                gnuplot_.setXRange(0, std::min(time, 72000.0));
                gnuplot_.setYRange(0.4, 0.6);
                gnuplot_.setXlabel("time [s]");
                gnuplot_.setYlabel("mole fraction mol/mol");
                gnuplot_.addDataSetToPlot(x_, y_, "CompTwo_left.dat", "w l t 'CompTwo left'");
                gnuplot_.addDataSetToPlot(x_, y2_, "CompTwo_right.dat", "w l t 'CompTwo right'");
                gnuplot_.plot("mole_fraction_CompTwo");

                gnuplot2_.resetPlot();
                gnuplot2_.setXRange(0, std::min(time, 72000.0));
                gnuplot2_.setYRange(0.0, 0.6);
                gnuplot2_.setXlabel("time [s]");
                gnuplot2_.setYlabel("mole fraction mol/mol");
                gnuplot2_.addDataSetToPlot(x_, y3_, "CompThree_left.dat", "w l t 'CompThree left'");
                gnuplot2_.addDataSetToPlot(x_, y4_, "CompThree_right.dat", "w l t 'CompThree right");
                gnuplot2_.plot("mole_fraction_CompThree");

                gnuplot3_.resetPlot();
                gnuplot3_.setXRange(0, std::min(time, 72000.0));
                gnuplot3_.setYRange(0.0, 0.6);
                gnuplot3_.setXlabel("time [s]");
                gnuplot3_.setYlabel("mole fraction mol/mol");
                gnuplot3_.addDataSetToPlot(x_, y5_, "CompOne_left.dat", "w l t 'CompOne left'");
                gnuplot3_.addDataSetToPlot(x_, y6_, "CompOne_right.dat", "w l t 'CompOne right'");
                gnuplot3_.plot("mole_fraction_CompOne");
           }
        }
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     * The units must be according to either using mole or mass fractions (mole/(m^2*s) or kg/(m^2*s)).
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        if (globalPos[0] < 0.5)
        {
           initialValues[FluidSystem::compOneIdx] = 0.0;
           initialValues[FluidSystem::compTwoIdx] = 0.50086;
           initialValues[FluidSystem::compThreeIdx] = 0.49914;
        }
        else
        {
           initialValues[FluidSystem::compOneIdx] = 0.50121;
           initialValues[FluidSystem::compTwoIdx] = 0.49879;
           initialValues[FluidSystem::compThreeIdx] = 0.0;
        }
        return initialValues;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;

    Dumux::GnuplotInterface<double> gnuplot_;
    Dumux::GnuplotInterface<double> gnuplot2_;
    Dumux::GnuplotInterface<double> gnuplot3_;

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> y2_;
    std::vector<double> y3_;
    std::vector<double> y4_;
    std::vector<double> y5_;
    std::vector<double> y6_;

    bool plotOutput_;
};

} // end namespace Dumux

#endif
