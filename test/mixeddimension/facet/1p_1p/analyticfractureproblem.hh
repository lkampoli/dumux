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
 * \brief A test problem for the one-dimensional single-phase fracture model
 */
#ifndef DUMUX_1P_ANALYTICAL_FRACTURE_PROBLEM_HH
#define DUMUX_1P_ANALYTICAL_FRACTURE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "analyticfracturespatialparams.hh"

namespace Dumux
{

//! Forward declaration of the problem class
template <class TypeTag>
class OnePFractureProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePFractureProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePCCFractureProblem, INHERITS_FROM(CCTpfaModel, OnePFractureProblem));

SET_PROP(OnePFractureProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Unit<Scalar> > type;
};

// Set the problem property
SET_TYPE_PROP(OnePFractureProblem, Problem, OnePFractureProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(OnePFractureProblem, SpatialParams, OnePFractureSpatialParams<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(OnePFractureProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(OnePFractureProblem, ProblemEnableGravity, false);
SET_BOOL_PROP(OnePCCFractureProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(OnePCCFractureProblem, EnableGlobalFluxVariablesCache, true);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 */
template <class TypeTag>
class OnePFractureProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

    enum
    {
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    OnePFractureProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_fracture";
        eps_ = 1e-6;
        aperture_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture);
        integrateAnalyticalSolution_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, IntegrateFracSol);
    }

    /*!
     * \brief The problem name.
     *        This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    { return aperture_; }


    /*!
     * \brief Return the sources within the domain.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        // we have only sources coming from the bulk domain
        auto sources = couplingManager().evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        sources /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return sources;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be used.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    { return PrimaryVariables(exact(element, scvf.ipGlobal())); }

    Scalar exact(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::cosh;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        return cos(x)*cosh(y);
    }

    Scalar exact(const Element& element, const GlobalPosition& globalPos) const
    {
        if (integrateAnalyticalSolution_)
        {
            // integrate and average the exact pressure
            auto upperEdge = globalPos;
            upperEdge[1] += aperture_/2.0;

            auto lowerEdge = globalPos;
            lowerEdge[1] -= aperture_/2.0;

            static const int order = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, FracSolQuadOrder);
            const auto& rule = Dune::QuadratureRules<Scalar, 1>::rule(Dune::GeometryType(1), order);

            std::vector<GlobalPosition> corners(2);
            corners[0] = lowerEdge;
            corners[1] = upperEdge;
            auto g = Dune::AffineGeometry<Scalar, 1, 2>(Dune::GeometryType(1), corners);

            Scalar result = 0.0;
            for (auto&& qp : rule)
            {
                const auto u = exact(g.global(qp.position()));
                result += u*qp.weight()*g.integrationElement(qp.position());
            }
            result /= aperture_;

            return result;
        }
        else
            return PrimaryVariables(exact(globalPos));
    }

    GlobalPosition exactGradient(const GlobalPosition& globalPos) const
    {
        using std::cos;
        using std::sin;
        using std::cosh;
        using std::sinh;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        GlobalPosition gradU;
        gradU[0] = -sin(x)*cosh(y);
        gradU[1] = cos(x)*sinh(y);

        return gradU;
    }

    Scalar exactFlux(const Element& element, const SubControlVolumeFace& scvf) const
    {
        static const Scalar kf = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FracturePermeability);

        if (integrateAnalyticalSolution_)
        {
            // integrate the exact flux over the aperture
            auto upperEdge = scvf.ipGlobal();
            upperEdge[1] += aperture_/2.0;

            auto lowerEdge = scvf.ipGlobal();
            lowerEdge[1] -= aperture_/2.0;

            static const int order = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, FracSolQuadOrder);
            const auto& rule = Dune::QuadratureRules<Scalar, 1>::rule(Dune::GeometryType(1), order);

            std::vector<GlobalPosition> corners(2);
            corners[0] = lowerEdge;
            corners[1] = upperEdge;
            auto g = Dune::AffineGeometry<Scalar, 1, 2>(Dune::GeometryType(1), corners);

            Scalar result = 0.0;
            for (auto&& qp : rule)
            {
                const auto gradU = exactGradient(g.global(qp.position()));
                const auto flux = -1.0*kf*(gradU*scvf.unitOuterNormal());
                result += flux*qp.weight()*g.integrationElement(qp.position());
            }
            return result;
        }
        else
        {
            const auto gradU = exactGradient(scvf.ipGlobal());
            return -1.0*kf*(gradU*scvf.unitOuterNormal());
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvf) const
    { return exactFlux(element, scvf); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); };

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class OutputModule>
    void addVtkOutputFields(OutputModule& outputModule) const
    {
        // create the required scalar fields
        auto& exactSol = outputModule.createScalarField("p_exact [N/m^2]", 0);

        for (const auto& element : elements(this->gridView()))
            exactSol[this->elementMapper().index(element)] = exact(element.geometry().center());
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    CouplingManager& couplingManager()
    { return *couplingManager_; }

private:
    std::string name_;
    Scalar eps_;
    Scalar aperture_;
    bool integrateAnalyticalSolution_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
