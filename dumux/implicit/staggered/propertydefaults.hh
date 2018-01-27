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
 * \ingroup CCTpfaProperties
 * \ingroup StaggeredModel
 * \file
 *
 * \brief Default properties for cell centered models
 */
#ifndef DUMUX_STAGGERED_PROPERTY_DEFAULTS_HH
#define DUMUX_STAGGERED_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/discretization/staggered/globalfvgeometry.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/implicit/staggered/properties.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/discretization/staggered/globalfluxvariablescache.hh>
#include <dumux/discretization/staggered/elementfluxvariablescache.hh>

#include <dumux/discretization/staggered/elementvolumevariables.hh>
#include <dumux/discretization/staggered/globalvolumevariables.hh>
#include <dumux/discretization/staggered/globalfacevariables.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/common/intersectionmapper.hh>

#include "assembler.hh"
#include "localresidual.hh"
#include "localjacobian.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "newtonconvergencewriter.hh"
#include "model.hh"
#include "primaryvariables.hh"

namespace Dumux {

// forward declarations
template<class TypeTag> class CCElementBoundaryTypes;

namespace Properties {
//! Set the corresponding discretization method property
SET_PROP(StaggeredModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Staggered;
};


SET_TYPE_PROP(StaggeredModel, BaseModel, StaggeredBaseModel<TypeTag>);


//! Set the default for the global finite volume geometry
SET_TYPE_PROP(StaggeredModel, GlobalFVGeometry, StaggeredGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(StaggeredModel, FVElementGeometry, StaggeredFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>);

//! The sub control volume
SET_PROP(StaggeredModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::CCSubControlVolume<ScvGeometry, IndexType> type;
};

SET_TYPE_PROP(StaggeredModel, GlobalFaceVars, Dumux::StaggeredGlobalFaceVariables<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(StaggeredModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(StaggeredModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The global current volume variables vector class
SET_TYPE_PROP(StaggeredModel, GlobalVolumeVariables, Dumux::StaggeredGlobalVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(StaggeredModel, GlobalFluxVariablesCache, Dumux::StaggeredGlobalFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The local jacobian operator
SET_TYPE_PROP(StaggeredModel, LocalJacobian, Dumux::StaggeredLocalJacobian<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(StaggeredModel, JacobianAssembler, Dumux::StaggeredAssembler<TypeTag>);

//! The local flux variables cache vector class
SET_TYPE_PROP(StaggeredModel, ElementFluxVariablesCache, Dumux::StaggeredElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(StaggeredModel, ElementVolumeVariables, Dumux::StaggeredElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! Set the BaseLocalResidual to StaggeredLocalResidual
SET_TYPE_PROP(StaggeredModel, BaseLocalResidual, Dumux::StaggeredLocalResidual<TypeTag>);

//! Set the BaseLocalResidual to StaggeredLocalResidual
SET_TYPE_PROP(StaggeredModel, IntersectionMapper, Dumux::ConformingGridIntersectionMapper<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(StaggeredModel, ImplicitIsBox, false);

SET_TYPE_PROP(StaggeredModel, NewtonController, StaggeredNewtonController<TypeTag>);

SET_INT_PROP(StaggeredModel, NumEqCellCenter, 1);
SET_INT_PROP(StaggeredModel, NumEqFace, 1);

SET_PROP(StaggeredModel, NumEq)
{
private:
    static constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
    static constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
public:
    static constexpr auto value = numEqCellCenter + numEqFace;
};

//! A vector of primary variables
SET_TYPE_PROP(StaggeredModel,
              CellCenterPrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEqCellCenter)>);
//! A vector of primary variables
SET_TYPE_PROP(StaggeredModel,
              FacePrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEqFace)>);

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(StaggeredModel,
              CellCenterSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables)>);

//! The type of a solution for the whole grid at a fixed time
SET_TYPE_PROP(StaggeredModel,
              FaceSolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables)>);

//! default property value for the solution vector only used for monolithic solver
SET_PROP(StaggeredModel, SolutionVector)
{
private:
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
public:
    typedef typename Dune::MultiTypeBlockVector<CellCenterSolutionVector, FaceSolutionVector> type;
};

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(StaggeredModel, JacobianMatrix)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };

public:
    // the sub-blocks
    using MatrixLittleBlockCCToCC = typename Dune::FieldMatrix<Scalar, numEqCellCenter, numEqCellCenter>; // 2x2
    using MatrixLittleBlockCCToFace = typename Dune::FieldMatrix<Scalar, numEqCellCenter, numEqFace>; // 2x1

    using MatrixLittleBlockFaceToFace = typename Dune::FieldMatrix<Scalar, numEqFace, numEqFace>; // 1x1
    using MatrixLittleBlockFaceToCC = typename Dune::FieldMatrix<Scalar, numEqFace, numEqCellCenter>; // 1x2

    // the BCRS matrices of the subproblems as big blocks
    using MatrixBlockCCToCC = typename Dune::BCRSMatrix<MatrixLittleBlockCCToCC>;
    using MatrixBlockCCToFace = typename Dune::BCRSMatrix<MatrixLittleBlockCCToFace>;


    using MatrixBlockFaceToFace = typename Dune::BCRSMatrix<MatrixLittleBlockFaceToFace>;
    using MatrixBlockFaceToCC = typename Dune::BCRSMatrix<MatrixLittleBlockFaceToCC>;

    // the row types
    using RowCellCenter = typename Dune::MultiTypeBlockVector<MatrixBlockCCToCC, MatrixBlockCCToFace>;
    using RowFace = typename Dune::MultiTypeBlockVector<MatrixBlockFaceToCC, MatrixBlockFaceToFace>;

    // the jacobian matrix
    using type = typename Dune::MultiTypeBlockMatrix<RowCellCenter, RowFace>;
};

//! Definition of the indices for cell center and face dofs in the global solution vector
SET_PROP(StaggeredModel, DofTypeIndices)
{
    using CellCenterIdx = Dune::index_constant<0>;
    using FaceIdx = Dune::index_constant<1>;
};

//! set default solver
// SET_TYPE_PROP(StaggeredModel, LinearSolver, Dumux::GSBiCGSTABBackend<TypeTag>);
SET_TYPE_PROP(StaggeredModel, LinearSolver, Dumux::UMFPackBackend<TypeTag>);

//! set the block level to 2, suitable for e.g. the Dune::MultiTypeBlockMatrix
SET_INT_PROP(StaggeredModel, LinearSolverPreconditionerBlockLevel, 2);

//! Boundary types at a single degree of freedom
SET_PROP(StaggeredModel, BoundaryTypes)
{
private:
    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };
public:
    using type = BoundaryTypes<numEqCellCenter + numEqFace>;
};

SET_PROP(StaggeredModel, PrimaryVariables)
{
private:
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterPrimaryVariables, FacePrimaryVariables>;
};

//! use the plain newton convergence writer by default
SET_TYPE_PROP(StaggeredModel, NewtonConvergenceWriter, StaggeredNewtonConvergenceWriter<TypeTag>);

//! Write separate vtp files for face variables by default
SET_BOOL_PROP(StaggeredModel, VtkWriteFaceData, true);

//! For compatibility
SET_BOOL_PROP(StaggeredModel, EnableInteriorBoundaries, false);

//! Set staggered vtk output module
SET_TYPE_PROP(StaggeredModel, VtkOutputModule, StaggeredVtkOutputModule<TypeTag>);

//! Set one or different base epsilons for the calculations of the localJacobian's derivatives
SET_PROP(StaggeredModel, BaseEpsilon)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    static constexpr Scalar dCCdCC = 1e-8;
    static constexpr Scalar dCCdFace = 1e-8;
    static constexpr Scalar dFacedCC = 1e-8;
    static constexpr Scalar dFacedFace = 1e-8;

public:
    static constexpr auto getEps()
    {
        return std::array<std::array<Scalar, 2>, 2>{{{dCCdCC, dCCdFace},
                                                     {dFacedCC, dFacedFace}}};
    }
};

} // namespace Properties

} // namespace Dumux

#endif