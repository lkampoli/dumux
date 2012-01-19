/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \file
 *
 * \brief Defines the properties required for the non-isothermal compositional
 * stokes BOX model.
 */
#ifndef DUMUX_STOKES2CNIPROPERTIES_HH
#define DUMUX_STOKES2CNIPROPERTIES_HH

#include "stokes2cniindices.hh"

#include <dumux/freeflow/stokes2c/stokes2cproperties.hh>
#include "stokes2cnivolumevariables.hh"
#include "stokes2cnifluxvariables.hh"

namespace Dumux
{
/*!
 * \addtogroup Stokes2cniModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class Stokes2cniModel;

template<class TypeTag>
class Stokes2cniLocalResidual;

template <class TypeTag>
class Stokes2cniVolumeVariables;

template <class TypeTag>
class Stokes2cniFluxVariables;

template <class TypeTag>
class Stokes2cniNewtonController;


////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the compositional stokes problems
NEW_TYPE_TAG(BoxStokes2cni, INHERITS_FROM(BoxStokes2c));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Stokes2cniIndices); //!< Enumerations for the compositional stokes models
NEW_PROP_TAG(NumComponents); //!< Number of components

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(BoxStokes2cni, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;
 public:
    static constexpr int value = 3 + dim;
};

//! Use the stokes2cni local jacobian operator for the compositional stokes model
SET_TYPE_PROP(BoxStokes2cni,
              LocalResidual,
              Stokes2cniLocalResidual<TypeTag>);

//! Use the stokes2cni specific newton controller for the compositional stokes model
SET_PROP(BoxStokes2cni, NewtonController)
{public:
    typedef Stokes2cniNewtonController<TypeTag> type;
};

//! the Model property
SET_TYPE_PROP(BoxStokes2cni, Model, Stokes2cniModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokes2cni, VolumeVariables, Stokes2cniVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokes2cni, FluxVariables, Stokes2cniFluxVariables<TypeTag>);

// the indices for the Stokes2cni model
SET_TYPE_PROP(BoxStokes2cni,
              Stokes2cniIndices,
              Stokes2cniCommonIndices<TypeTag>);

//! Set the number of components to 2
SET_INT_PROP(BoxStokes2cni, NumComponents, 2);
}

}
#endif
