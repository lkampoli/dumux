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
 * \brief A test problem for the one-phase blood flow model:
 * Blood is flowing through a 1d network grid.
 */
#ifndef DUMUX_BLOOD_FLOW_PROBLEM_HH
#define DUMUX_BLOOD_FLOW_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>
#include "bloodflowspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class BloodFlowProblem;

namespace Properties
{
NEW_TYPE_TAG(BloodFlowProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(BloodFlowCCProblem, INHERITS_FROM(CCTpfaModel, BloodFlowProblem));

// Set the grid type
SET_TYPE_PROP(BloodFlowProblem, Grid, Dune::FoamGrid<1, 3>);

// Set the grid parameter group
SET_STRING_PROP(BloodFlowProblem, GridParameterGroup, "VesselGrid");

// Set the problem property
SET_TYPE_PROP(BloodFlowProblem, Problem, BloodFlowProblem<TypeTag>);

SET_PROP(BloodFlowProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, Constant<TypeTag, Scalar>>;
};

// Set the spatial parameters
SET_TYPE_PROP(BloodFlowProblem, SpatialParams, BloodFlowSpatialParams<TypeTag>);

// evaluate analytical permeability field at scvf center
SET_BOOL_PROP(BloodFlowProblem, EvaluatePermeabilityAtScvfCenter, true);

// disable gravity
SET_BOOL_PROP(BloodFlowProblem, ProblemEnableGravity, false);

}

/*!
 * \ingroup OneDStokes
 * \ingroup ImplicitTestProblems
 * \brief Exact solution 1D-3D
 */
template <class TypeTag>
class BloodFlowProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // copy some indices for convenience
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

public:
    BloodFlowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //read parameters from input file
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_1d";
        p_in_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryConditions1D, PressureInput);
        delta_p_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BoundaryConditions1D, DeltaPressure);
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        ParentType::init();
        exactPressure_.resize(this->model().numDofs());
        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
                exactPressure_[scv.dofIndex()] = exactSolution(scv.dofPosition());
        }
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolutionVector& elemSol) const
    {
        const auto eIdx = this->elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \brief Returns true if a restart file should be written to disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 37.0; } // Body temperature

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        if(globalPos[2] > 0.5)
            values[pressureIdx] = p_in_;
        else
            values[pressureIdx] = p_in_-delta_p_;
        return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const auto& bulkVolVars = this->couplingManager().bulkVolVars(source.id());
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pressureIdx];

        // calculate the source
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        const Scalar sourceValue = beta*(bulkVolVars.pressure() - pressure1D)*bulkVolVars.density();

        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    //! The exact pressure solution
    Scalar exactSolution(const GlobalPosition &globalPos) const
    { return 1 + globalPos[2]; }

    //! Called after every time step
    //! Output the total global exchange term
    void postTimeStep()
    {
        ParentType::postTimeStep();

        PrimaryVariables source(0.0);

        if (!(this->timeManager().time() < 0.0))
        {
            for (const auto& element : elements(this->gridView()))
            {
                auto fvGeometry = localView(this->model().globalFvGeometry());
                fvGeometry.bindElement(element);

                auto elemVolVars = localView(this->model().curGlobalVolVars());
                elemVolVars.bindElement(element, fvGeometry, this->model().curSol());

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                    pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                    source += pointSources;
                }
            }
        }

        std::cout << "Global integrated source (1D): " << source << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void addVtkOutputFields(VtkOutputModule<TypeTag>& outputModule) const
    {
        auto& p = outputModule.createScalarField("exact pressure", dofCodim);
        p = exactPressure_;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::vector<Scalar> exactPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    Scalar p_in_;
    Scalar delta_p_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
