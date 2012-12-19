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
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxModel
 * \file
 *
 * \brief Default properties for box models
 */
#ifndef DUMUX_BOX_PROPERTY_DEFAULTS_HH
#define DUMUX_BOX_PROPERTY_DEFAULTS_HH

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/backend/novlpistlsolverbackend.hh>
#include <dumux/linear/amgbackend.hh>
#endif

#include <dumux/implicit/common/implicitpropertydefaults.hh>
#include "boxassembler.hh"
#include "boxfvelementgeometry.hh"
#include "boxelementboundarytypes.hh"
#include "boxlocalresidual.hh"
#include "boxelementvolumevariables.hh"
#include "boxproperties.hh"

namespace Dumux {

// forward declarations
template<class TypeTag> class BoxModel;
template<class TypeTag> class BoxLocalResidual;
template<class TypeTag> class BoxElementBoundaryTypes;
template<class TypeTag> class BoxElementVolumeVariables;
template<class TypeTag> class BoxFVElementGeometry;
template<class TypeTag> class AMGBackend;

namespace Properties {
//! Set the default for the FVElementGeometry
SET_TYPE_PROP(BoxModel, FVElementGeometry, Dumux::BoxFVElementGeometry<TypeTag>);

//! Disable evaluation of shape function gradients at the sub-control volume center by default
// The shape function gradients at the sub-control volume center are currently only
// needed for the Stokes and the linear elastic models
SET_BOOL_PROP(BoxModel, EvalGradientsAtSCVCenter, false);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(BoxModel, ElementBoundaryTypes, Dumux::BoxElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(BoxModel, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, Dumux::BoxLocalResidual<TypeTag>);

//! An array of secondary variable containers
SET_TYPE_PROP(BoxModel, ElementVolumeVariables, Dumux::BoxElementVolumeVariables<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(BoxModel, JacobianAssembler, Dumux::BoxAssembler<TypeTag>);

//! disable two-point-flux by default
SET_BOOL_PROP(BoxModel, ImplicitUseTwoPointFlux, false);

//! indicate that this is a box discretization
SET_BOOL_PROP(BoxModel, ImplicitIsBox, true);

#if HAVE_DUNE_PDELAB
//! use the local FEM space associated with cubes by default
SET_PROP(BoxModel, ImplicitLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type;
};

SET_PROP(BoxModel, ImplicitPDELabBackend)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GridOperator> type;
};

/*SET_PROP(BoxModel, SolutionVector)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};*/
#endif

} // namespace Properties
} // namespace Dumux

#endif
