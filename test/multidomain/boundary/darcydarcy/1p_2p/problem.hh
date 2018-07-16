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
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_ONEP_SUB_TEST_PROBLEM_HH
#define DUMUX_ONEP_SUB_TEST_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>
#include "spatialparams.hh"

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief  Multidomain test problem for the incompressible one-phase model
 * \todo doc me!
 */
template<class TypeTag, std::size_t tag>
class OnePTestProblem
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto domainIdx = Dune::index_constant<tag>{};

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        injectionRate_ = getParam<double>("Problem.InjectionRate");
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        if (couplingManager_->isCoupled(domainIdx, scvf))
            values.setAllCouplingNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto bcTypes = boundaryTypes(element, scvf);
        if (bcTypes.hasCouplingNeumann())
            values[Indices::conti0EqIdx] = couplingManager_->advectiveFluxCoupling(domainIdx, element, fvGeometry, elemVolVars, scvf, FluidSystem::phase0Idx);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e6;
        if (NumEqVector::dimension > 1)
            values[Indices::pressureIdx + 1] = 1e-10; // fully saturated
        return values;
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the values method of the point source
     * has to return the absolute rate values in units
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    template<class PointSource>
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        NumEqVector values(0.0);
        // if this is the 2p problem inject methane
        if (NumEqVector::dimension > 1)
            values[Indices::conti0EqIdx + 1] = injectionRate_; // kg/s
        pointSources.emplace_back(GlobalPosition(0.0), values);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos); }


private:
    Scalar injectionRate_;
    std::shared_ptr<CouplingManager> couplingManager_;
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif