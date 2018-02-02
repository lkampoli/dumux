// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 *
 * \brief Tutorial problem for a fully coupled twophase box model.
 */
#ifndef DUMUX_EXERCISE_THREE_A_PROBLEM_HH
#define DUMUX_EXERCISE_THREE_A_PROBLEM_HH

// The numerical model
#include <dumux/porousmediumflow/2p/model.hh>

//The box discretization
#include <dumux/discretization/box/properties.hh>

// The base porous media box problem
#include <dumux/porousmediumflow/problem.hh>

// Spatially dependent parameters
#include "spatialparams.hh"

// The water component
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>

// The components that will be created in this exercise
#include "components/myincompressiblecomponent.hh"
// #include "components/mycompressiblecomponent.hh"

// We will only have liquid phases Here
#include <dumux/material/fluidsystems/liquidphase.hh>

// The two-phase immiscible fluid system
#include <dumux/material/fluidsystems/2pimmiscible.hh>

namespace Dumux{
// Forward declaration of the problem class
template <class TypeTag> class ExerciseThreeProblemTwoP;

namespace Properties {
// Create a new type tag for the problem
NEW_TYPE_TAG(ExerciseThreeTwoPTypeTag, INHERITS_FROM(TwoP, ExerciseThreeSpatialParams));
NEW_TYPE_TAG(ExerciseThreeBoxTwoPTypeTag, INHERITS_FROM(BoxModel, ExerciseThreeTwoPTypeTag));

// Set the "Problem" property
SET_TYPE_PROP(ExerciseThreeTwoPTypeTag, Problem, ExerciseThreeProblemTwoP<TypeTag>);

// Set grid and the grid creator to be used
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(ExerciseThreeTwoPTypeTag, Grid, Dune::ALUGrid</*dim=*/2, 2, Dune::cube, Dune::nonconforming>);
#elif HAVE_UG
SET_TYPE_PROP(ExerciseThreeTwoPTypeTag, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(ExerciseThreeTwoPTypeTag, Grid, Dune::YaspGrid<2>);
#endif // HAVE_DUNE_ALUGRID

// we use the immiscible fluid system here
SET_PROP(ExerciseThreeTwoPTypeTag, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using TabulatedH2O = TabulatedComponent<Scalar, H2O<Scalar>>;
    using WettingPhase = typename FluidSystems::LiquidPhase<Scalar, TabulatedH2O>;
    /*!
     * Uncomment first line and comment second line for using the incompressible component
     * Uncomment second line and comment first line for using the compressible component
     */
    using NonWettingPhase = typename FluidSystems::LiquidPhase<Scalar, MyIncompressibleComponent<Scalar> >;
    // using NonWettingPhase = typename FluidSystems::LiquidPhase<Scalar, MyCompressibleComponent<Scalar> >;

public:
    using type = typename FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonWettingPhase>;
};

}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief  Tutorial problem for a fully coupled twophase box model.
 */
template <class TypeTag>
class ExerciseThreeProblemTwoP : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // Grid dimension
    enum { dim = GridView::dimension,
           dimWorld = GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

    // Dumux specific types
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    ExerciseThreeProblemTwoP(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
        , eps_(3e-6)
    {
#if !(HAVE_DUNE_ALUGRID || HAVE_UG)
      std::cout << "If you want to use simplices instead of cubes, install and use dune-ALUGrid or UGGrid." << std::endl;
#endif // !(HAVE_DUNE_ALUGRID || HAVE_UG)

      // initialize the tables for the water properties
      std::cout << "Initializing the tables for the water properties" << std::endl;
      TabulatedComponent<Scalar, H2O<Scalar>>::init(/*tempMin=*/273.15,
                                                    /*tempMax=*/623.15,
                                                    /*numTemp=*/100,
                                                    /*pMin=*/0.0,
                                                    /*pMax=*/20e6,
                                                    /*numP=*/200);

      // set the depth of the bottom of the reservoir
      depthBOR_ = this->fvGridGeometry().bBoxMax()[dimWorld-1];

        // name of the problem and output file
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }


    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar temperature() const
    { return 283.15; }

     /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
         BoundaryTypes bcTypes;

        if (globalPos[0] < eps_ || globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
           bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;

        priVars[Indices::pwIdx] = 200.0e3 + 9.81*1000*(depthBOR_ - globalPos[dimWorld-1]); // 200 kPa = 2 bar
        priVars[Indices::snIdx] = 0.0; // 0 % oil saturation on left boundary
       return priVars;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        // initialize values to zero, i.e. no-flow Neumann boundary conditions
        PrimaryVariables values(0.0);

        Scalar up = this->fvGridGeometry().bBoxMax()[dimWorld-1];
        // extraction of oil on the right boundary for approx. 1.e6 seconds
        if (globalPos[dimWorld-1] > up - eps_ && globalPos[0] > 20 && globalPos[0] < 40) {
            // oil outflux of 30 g/(m * s) on the right boundary.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = -3e-2;
        } else {
            // no-flow on the remaining Neumann-boundaries.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 0;
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const

    {
        PrimaryVariables values(0.0);

        values[Indices::pwIdx] = 200.0e3 + 9.81*1000*(depthBOR_ - globalPos[dimWorld-1]); // 200 kPa = 2 bar
        values[Indices::snIdx] = 0.0;

        return values;
    }
    // \}

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        return values;
    }

private:
    // small epsilon value
    Scalar eps_;
    std::string name_; //! Problem name

    // depth at the bottom of the reservoir
    Scalar depthBOR_;
};
}

#endif