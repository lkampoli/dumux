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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup MultidomainModel
 *
 * \brief Defines the properties required for the coupled 2cstokes2p2c model.
 */

#ifndef DUMUX_TWOCSTOKESTWOPTWOC_PROPERTIES_HH
#define DUMUX_TWOCSTOKESTWOPTWOC_PROPERTIES_HH

#include <dumux/multidomain/propertydefaults.hh>

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

//! The type tags for the coupled 2cstokes2p2c model
NEW_TYPE_TAG(TwoCStokesTwoPTwoC, INHERITS_FROM(MultiDomain));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(BoundaryLayerModel); //!< Type of the used boundary layer model
NEW_PROP_TAG(MassTransferModel); //!< Type of the used mass transfer model

} // end namespace properties

} // end namespace Dumux


#endif