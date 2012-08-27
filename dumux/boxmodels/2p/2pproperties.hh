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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup TwoPBoxModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the twophase BOX model.
 */

#ifndef DUMUX_2P_PROPERTIES_HH
#define DUMUX_2P_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>

namespace Dumux
{



////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
NEW_TYPE_TAG(BoxTwoP, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableGravity); //!< DEPRECATED Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the weight of the upwind direction in the mass conservation equations
NEW_PROP_TAG(MassUpwindWeight); //!< DEPRECATED The value of the weight of the upwind direction in the mass conservation equations
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(MobilityUpwindWeight); //!< DEPRECATED Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(Formulation);   //!< The formulation of the model
NEW_PROP_TAG(TwoPIndices); //!< DEPRECATED Enumerations for the 2p models
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(SpatialParameters); //!< DEPRECATED The type of the spatial parameters
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(WettingPhase); //!< The wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase); //!< The non-wetting phase for two-phase models
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about the phases
NEW_PROP_TAG(FluidState); //!<The phases state
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
NEW_PROP_TAG(EnableVelocityOutput); //!< DEPRECATED Returns whether vertex velocity vectors are written into the vtk output
}

}

#endif
