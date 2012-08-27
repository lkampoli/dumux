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
 * \ingroup OnePhase
 * \ingroup IMPETProperties
 * \file
 *
 * \brief Defines the properties required for the single phase sequential model.
 */

#ifndef DUMUX_1PPROPERTIES_HH
#define DUMUX_1PPROPERTIES_HH

//Dumux-includes
#include <dumux/decoupled/common/decoupledproperties.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////
template <class TypeTag>
class CellData1P;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase problem
NEW_TYPE_TAG(DecoupledOneP, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( SpatialParameters ); //!< The type of the spatial parameters object
NEW_PROP_TAG( SpatialParams ); //!< The type of the spatial parameters object
NEW_PROP_TAG( ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( EnableGravity); //!< DEPRECATED Returns whether gravity is considered in the problem
NEW_PROP_TAG( Fluid ); //!< The fluid for one-phase models
NEW_PROP_TAG( Indices ); //!< Set of indices for the one-phase model
NEW_PROP_TAG( CellData ); //!< The cell data storage class
}
}

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/1p/cellData1p.hh>
#include <dumux/decoupled/1p/1pindices.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Set number of equations to 1 for isothermal one-phase models
SET_INT_PROP(DecoupledOneP, NumEq, 1);

//! Set number of phases to 1 for one-phase models
SET_INT_PROP(DecoupledOneP, NumPhases, 1)
;
//!< Each phase consists of 1 pure component
SET_INT_PROP(DecoupledOneP, NumComponents, 1);

//! Chose the set of indices for the one-phase formulation
SET_TYPE_PROP(DecoupledOneP, Indices, DecoupledOnePCommonIndices);

//! Set general decoupled VariableClass as default
SET_TYPE_PROP(DecoupledOneP, Variables, VariableClass<TypeTag>);

//! Set standart CellData of immiscible one-phase models as default
SET_TYPE_PROP(DecoupledOneP, CellData, CellData1P<TypeTag>);

//! DEPRECATED SpatialParameters property
SET_TYPE_PROP(DecoupledOneP, SpatialParameters, typename GET_PROP_TYPE(TypeTag, SpatialParams));

//Has to be removed if DEPRECATED EnableGravity is removed!
SET_BOOL_PROP(DecoupledOneP, ProblemEnableGravity, GET_PROP_VALUE(TypeTag, EnableGravity));
}
}
#endif
