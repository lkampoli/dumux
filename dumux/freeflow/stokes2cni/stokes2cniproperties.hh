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
 * \ingroup BoxStokes2cniModel
 *
 * \file
 *
 * \brief Defines the additional properties required for the non-isothermal compositional
 * Stokes box model.
 */
#ifndef DUMUX_STOKES2CNI_PROPERTIES_HH
#define DUMUX_STOKES2CNI_PROPERTIES_HH

#include <dumux/freeflow/stokes2c/stokes2cproperties.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the non-isothermal compositional Stokes problems
NEW_TYPE_TAG(BoxStokes2cni, INHERITS_FROM(BoxStokes2c));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Stokes2cniIndices); //!< DEPRECATED Enumerations for the compositional Stokes models
NEW_PROP_TAG(NumComponents); //!< Number of components
}

}
#endif
