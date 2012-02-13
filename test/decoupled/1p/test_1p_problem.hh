// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
*   Copyright (C) 2007-2008 by Markus Wolff                                  *
*   Institute for Modelling Hydraulic and Environmental Systems              *
*   University of Stuttgart, Germany                                         *
*   email: <givenname>.<name>@iws.uni-stuttgart.de                           *
*                                                                            *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/
/*!
 * \file
 *
 * \brief test problem for the decoupled one-phase model.
 */
#ifndef DUMUX_TEST_1P_PROBLEM_HH
#define DUMUX_TEST_1P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/unit.hh>

#include <dumux/decoupled/1p/diffusion/fv/fvpressureproperties1p.hh>
#include <dumux/decoupled/1p/diffusion/diffusionproblem1p.hh>
#include <dumux/decoupled/common/fv/fvvelocity.hh>

#include "test_1p_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestProblemOneP;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TestProblemOneP, INHERITS_FROM(FVPressureOneP));

// set the GridCreator property
SET_TYPE_PROP(TestProblemOneP, GridCreator, CubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(TestProblemOneP, Grid)
{
        typedef Dune::YaspGrid<2> type;
//    typedef Dune::SGrid<2, 2> type;
};

// Set the wetting phase
SET_PROP(TestProblemOneP, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the spatial parameters
SET_TYPE_PROP(TestProblemOneP, SpatialParameters, Dumux::TestOnePSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TestProblemOneP, EnableGravity, false);

//Set the problem
SET_TYPE_PROP(TestProblemOneP, Problem, Dumux::TestProblemOneP<TypeTag>);


SET_INT_PROP(TestProblemOneP, LinearSolverVerbosity, 1);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the decoupled one-phase model.
 */
template<class TypeTag>
class TestProblemOneP: public DiffusionProblem1P<TypeTag >
{
    typedef DiffusionProblem1P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;


public:
    TestProblemOneP(TimeManager &timeManager, const GridView &gridView) :
        ParentType(gridView), velocity_(*this)
    {
        delta_ = 1e-3 ;

        try
        {
            if (ParameterTree::tree().hasKey("delta"))
                delta_       = GET_RUNTIME_PARAM(TypeTag, Scalar, delta_);
            int numRefine;
            numRefine = GET_RUNTIME_PARAM(TypeTag, int, numRefine);
            GridCreator::grid().globalRefine(numRefine);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

        this->spatialParameters().setDelta(delta_);
    }

    /*!
    * \name Problem parameters
    */
    // \{

    /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
    const char *name() const
    {
        return "test_1p";
    }

    bool shouldWriteRestartFile() const
    { return false; }

    void addOutputVtkFields()
    {
        velocity_.calculateVelocity();
        velocity_.addOutputVtkFields(this->resultWriter());
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperature(const Element& element) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressure(const Element& element) const
    {
        return 1e5; // -> 10°C
    }

    //!source term [kg/(m^3 s)]
    void sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
        {
        double pi = 4.0*atan(1.0);
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        double ux = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        double uy = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);
        double kxx = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        double kxy = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        double kyy = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
        double f0 = sin(pi*globalPos[0])*sin(pi*globalPos[1])*pi*pi*(1.0 + delta_)*(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])
        + cos(pi*globalPos[0])*sin(pi*globalPos[1])*pi*(1.0 - 3.0*delta_)*globalPos[0]
                                                                                    + cos(pi*globalPos[1])*sin(pi*globalPos[0])*pi*(1.0 - 3.0*delta_)*globalPos[1]
                                                                                                                                                                + cos(pi*globalPos[1])*cos(pi*globalPos[0])*2.0*pi*pi*(1.0 - delta_)*globalPos[0]*globalPos[1];
        values = (f0 + 2.0*(globalPos[0]*(kxx*ux + kxy*uy) + globalPos[1]*(kxy*ux + kyy*uy)))/rt;
        }

    /*!
    * \brief Returns the type of boundary condition.
    *
    * BC can be dirichlet (pressure) or neumann (flux).
    */
    void boundaryTypes(BoundaryTypes &bcType,
            const Intersection& intersection) const
    {
        bcType.setAllDirichlet();
    }

    //! return dirichlet condition  (pressure, [Pa])
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = exact(globalPos);
    }


    //! return neumann condition  (flux, [kg/(m^2 s)])
    void neumann(PrimaryVariables &values, const Intersection& intersection) const
        {
        values = 0;
        }

private:
    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
        {
        Dune::FieldVector<Scalar,dim> grad(0);
        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
        }

    double delta_;
    Dumux::FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity) > velocity_;
};
} //end namespace

#endif
