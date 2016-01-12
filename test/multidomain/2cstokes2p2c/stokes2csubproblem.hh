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
/**
 * \file
 * \brief Isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 */
#ifndef DUMUX_STOKES2C_SUBPROBLEM_HH
#define DUMUX_STOKES2C_SUBPROBLEM_HH

#include <dumux/freeflow/stokesnc/stokesncmodel.hh>
#include <dumux/multidomain/2cstokes2p2c/stokesnccouplinglocalresidual.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cSubProblem;

namespace Properties
{
NEW_TYPE_TAG(Stokes2cSubProblem,
             INHERITS_FROM(BoxStokesnc, SubDomain));

// Set the problem property
SET_TYPE_PROP(Stokes2cSubProblem, Problem, Dumux::Stokes2cSubProblem<TypeTag>);

// Set the property for the material parameters by extracting it from the material law.
SET_TYPE_PROP(Stokes2cSubProblem,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

// Use the local residual extended for the coupling the local residual extended for the coupling
SET_TYPE_PROP(Stokes2cSubProblem, LocalResidual, StokesncCouplingLocalResidual<TypeTag>);

// Used the fluid system from the coupled problem
SET_TYPE_PROP(Stokes2cSubProblem,
              FluidSystem,
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));

// Use formulation based on mass fractions
SET_BOOL_PROP(Stokes2cSubProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(Stokes2cSubProblem, ProblemEnableGravity, false);

// Switch inertia term off
SET_BOOL_PROP(Stokes2cSubProblem, EnableNavierStokes, false);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief Isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 *
 * The Stokes subdomain is sized 0.25m times 0.25m. The boundary conditions
 * for the momentum balances are all set to Dirichlet, except on the right
 * boundary, where outflow conditions are set. The mass balance receives
 * outflow BCs, which are replaced in the localresidual by the sum
 * of the two momentum balances. On the right boundary Dirichlet BCs are
 * set for the pressure.
 *
 * This sub problem uses the \ref StokesncModel. It is part of the
 * 2cstokes2p2c model and is combined with the 2p2csubproblem for
 * the Darcy domain.
 */
template <class TypeTag>
class Stokes2cSubProblem : public StokesProblem<TypeTag>
{
    typedef Stokes2cSubProblem<TypeTag> ThisType;
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,

        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance

        transportEqIdx = Indices::transportEqIdx //!< Index of the transport equation (massfraction)
    };
    enum { // primary variable indices
            pressureIdx = Indices::pressureIdx,
            velocityXIdx = Indices::velocityXIdx,
            velocityYIdx = Indices::velocityYIdx,
            velocityZIdx = Indices::velocityZIdx,
            massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };
    enum { phaseIdx = Indices::phaseIdx };
    enum { numComponents = Indices::numComponents };
    enum {
        transportCompIdx = Indices::transportCompIdx, //!< water component index
        phaseCompIdx = Indices::phaseCompIdx          //!< air component index
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    /*!
     * \brief The sub-problem for the Stokes subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    Stokes2cSubProblem(TimeManager &timeManager, const GridView gridView)
        : ParentType(timeManager, gridView)
    {
        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, LowerLeftX);
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightX);
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightY);
        runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX); // first part of the interface without coupling

        refVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);
        refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefPressure);
        refMassfrac_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
        refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);

        sinusVAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelAmplitude);
        sinusVPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelPeriod);
        sinusPAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressureAmplitude);
        sinusPPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressurePeriod);
        sinusXAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationAmplitude);
        sinusXPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationPeriod);
        sinusTAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperatureAmplitude);
        sinusTPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperaturePeriod);

        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc ImplicitProblem::bBoxMax()
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NameFF); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature, which can
     * be set in the parameter file.
     */
    Scalar temperature() const
    {
        return refTemperature_;
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        const Scalar time = this->timeManager().time();

        values.setAllDirichlet();

        if (onUpperBoundary_(globalPos))
            values.setNeumann(transportEqIdx);

        // Left inflow boundaries should be Neumann, otherwise the
        // evaporative fluxes are much more grid dependent
        if (onLeftBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onRightBoundary_(globalPos))
        {
            values.setAllOutflow();

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onLowerBoundary_(globalPos))
        {
            values.setAllDirichlet();
            values.setNeumann(transportEqIdx);

            if (globalPos[0] > runUpDistanceX_-eps_ && time > initializationTime_)
            {
                values.setAllCouplingDirichlet();
                values.setCouplingNeumann(momentumXIdx);
                values.setCouplingNeumann(momentumYIdx);
            }
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for stokes
        values.setOutflow(massBalanceIdx);

        // set pressure at one point, do NOT specify this
        // if the Darcy domain has a Dirichlet condition for pressure
        if (onRightBoundary_(globalPos))
        {
            if (time > initializationTime_)
                values.setDirichlet(pressureIdx);
            else
                if (!onLowerBoundary_(globalPos) && !onUpperBoundary_(globalPos))
                    values.setDirichlet(pressureIdx);
        }
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = refPressure()
                + density*this->gravity()[1]*(globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param values The Neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^{\textrm{dim}-1} \cdot s )] \f$
     * \param globalPos The global position
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.;

        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);
        const Scalar xVelocity = xVelocity_(globalPos);

        if (onLeftBoundary_(globalPos)
                && globalPos[1] > bBoxMin_[1] && globalPos[1] < bBoxMax_[1])
        {
            // rho*v*X at inflow
            values[transportEqIdx] = -xVelocity * density * refMassfrac();
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        // The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }
    // \}

    //! \brief Returns the reference velocity.
    const Scalar refVelocity() const
    { return refVelocity_ + variation_(sinusVAmplitude_, sinusVPeriod_); }

    //! \brief Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_ + variation_(sinusPAmplitude_, sinusPPeriod_); }

    //! \brief Returns the reference mass fraction.
    const Scalar refMassfrac() const
    { return refMassfrac_ + variation_(sinusXAmplitude_, sinusXPeriod_); }

    //! \brief Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_+ variation_(sinusTAmplitude_, sinusTPeriod_); }

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.;

        values[pressureIdx] = refPressure()
                + density*this->gravity()[1]*(globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
    }

    //! \brief set the profile of the inflow velocity (horizontal direction)
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        const Scalar vmax = refVelocity();
        return  4*vmax*(globalPos[1] - bBoxMin_[1])*(bBoxMax_[1] - globalPos[1])
                / (height_()*height_()) + 0.00134;
    }

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState) const
    {
        fluidState.setTemperature(refTemperature());
        fluidState.setPressure(phaseIdx, refPressure());
        // setMassFraction() has only to be called 1-numComponents times
        fluidState.setMassFraction(phaseIdx, transportCompIdx, refMassfrac());
    }

    // can be used for the variation of a boundary condition
    const Scalar variation_(const Scalar amplitude, const Scalar period) const
    { return sin(2*M_PI*this->timeManager().time()/period) * amplitude; }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    // the height of the free-flow domain
    const Scalar height_() const
    { return bBoxMax_[1] - bBoxMin_[1]; }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMassfrac_;
    Scalar refTemperature_;

    Scalar sinusVAmplitude_;
    Scalar sinusVPeriod_;
    Scalar sinusPAmplitude_;
    Scalar sinusPPeriod_;
    Scalar sinusXAmplitude_;
    Scalar sinusXPeriod_;
    Scalar sinusTAmplitude_;
    Scalar sinusTPeriod_;

    Scalar runUpDistanceX_;
    Scalar initializationTime_;
};
} //end namespace Dumux

#endif // DUMUX_STOKES2C_SUBPROBLEM_HH
