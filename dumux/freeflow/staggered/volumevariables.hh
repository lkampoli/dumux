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
 * \brief Quantities required by the one-phase fully implicit model defined on a vertex.
 */
#ifndef DUMUX_1P_VOLUME_VARIABLES_HH
#define DUMUX_1P_VOLUME_VARIABLES_HH

#include "properties.hh"
#include <dumux/discretization/staggered/volumevariables.hh>

#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class NavierStokesVolumeVariables : public StaggeredVolumeVariables<TypeTag>
{
    using ParentType = StaggeredVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(priVars, problem, element, scv);

        completeFluidState(priVars, problem, element, scv, fluidState_);
    };

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(priVars, problem, element, scv);
        fluidState.setTemperature(t);
        fluidState.setSaturation(/*phaseIdx=*/0, 1.);

        fluidState.setPressure(/*phaseIdx=*/0, priVars[Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setViscosity(/*phaseIdx=*/0, value);

        // compute and set the enthalpy
        value = ParentType::enthalpy(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setEnthalpy(/*phaseIdx=*/0, value);
    }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return the saturation
     */
    Scalar saturation(int phaseIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Return the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the mobility \f$\mathrm{[1/(Pa s)]}\f$.
     *
     * This function enables the use of ImplicitDarcyFluxVariables
     * with the 1p fully implicit model, ALTHOUGH the term mobility is
     * usually not employed in the one phase context.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx = 0) const
    { return 1.0/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the fluid state of the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

protected:
    FluidState fluidState_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
