// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Common
 *
 * \brief _Declares_ all properties used in Dumux.
 * \note Include this to forward declare properties in your headers.
 */

#ifndef DUMUX_PROPERTIES_HH
#define DUMUX_PROPERTIES_HH

// explicitly guard the include so that the property system
// header doesn't need to be opened and checked all the time
#ifndef DUMUX_PROPERTY_SYSTEM_HH
#include <dumux/common/properties/propertysystem.hh>

// per default, we do not allow the old property macros
// remove this after release 3.2
#ifndef DUMUX_ENABLE_OLD_PROPERTY_MACROS
#define DUMUX_ENABLE_OLD_PROPERTY_MACROS 0
#endif

// remove this after release 3.2 to remove macros completely
#if DUMUX_ENABLE_OLD_PROPERTY_MACROS
#include <dumux/common/properties/propertysystemmacros.hh>
#endif // DUMUX_ENABLE_OLD_PROPERTY_MACROS
#endif // DUMUX_PROPERTY_SYSTEM_HH

namespace Dumux {
namespace Properties {

///////////////////////////////////////
// Basic properties of numeric models:
///////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct Scalar { using type = UndefinedProperty; };                 //!< Property to specify the type of scalar values.
template<class TypeTag, class MyTypeTag>
struct ModelDefaultParameters { using type = UndefinedProperty; }; //!< Property which defines the group that is queried for parameters by default
template<class TypeTag, class MyTypeTag>
struct Grid { using type = UndefinedProperty; };                   //!< The DUNE grid type
template<class TypeTag, class MyTypeTag>
struct PrimaryVariables { using type = UndefinedProperty; };       //!< A vector of primary variables
template<class TypeTag, class MyTypeTag>
struct NumEqVector { using type = UndefinedProperty; };            //!< A vector of size number equations that can be used for Neumann fluxes, sources, residuals, ...
template<class TypeTag, class MyTypeTag>
struct GridView { using type = UndefinedProperty; };               //!< The type of the grid view according to the grid type
template<class TypeTag, class MyTypeTag>
struct ModelTraits { using type = UndefinedProperty; };            //!< Traits class encapsulating model specifications
template<class TypeTag, class MyTypeTag>
struct BaseModelTraits { using type = UndefinedProperty; };        //!< Model traits to be used as a base for nonisothermal, mineralization ... models
template<class TypeTag, class MyTypeTag>
struct Problem { using type = UndefinedProperty; };                //!< Property to specify the type of a problem which has to be solved
template<class TypeTag, class MyTypeTag>
struct PointSource { using type = UndefinedProperty; };            //!< Property defining the type of point source used
template<class TypeTag, class MyTypeTag>
struct PointSourceHelper { using type = UndefinedProperty; };      //!< Property defining the class that computes which sub control volume point sources belong to
// TODO: Remove deprecated property VtkOutputFields
template<class TypeTag, class MyTypeTag>
struct VtkOutputFields { using type = UndefinedProperty; };        //!< A class helping models to define default vtk output parameters
template<class TypeTag, class MyTypeTag>
struct IOFields { using type = UndefinedProperty; };               //!< A class helping models to define input and output fields
template<class TypeTag, class MyTypeTag>
struct BaseLocalResidual { using type = UndefinedProperty; };      //!< The type of the base class of the local residual (specific to a discretization scheme)
template<class TypeTag, class MyTypeTag>
struct JacobianMatrix { using type = UndefinedProperty; };         //!< Type of the global jacobian matrix
template<class TypeTag, class MyTypeTag>
struct SolutionVector { using type = UndefinedProperty; };         //!< Vector containing all primary variable vector of the grid
template<class TypeTag, class MyTypeTag>
struct BoundaryTypes { using type = UndefinedProperty; };          //!< Stores the boundary types of a single degree of freedom

//! The type of the local residual function, i.e. the equation to be solved. Must inherit
//! from the BaseLocalResidual property and fulfill its interfaces.
template<class TypeTag, class MyTypeTag>
struct LocalResidual { using type = UndefinedProperty; };

//! TODO: Remove this property as soon as the decoupled models are integrated
template<class TypeTag, class MyTypeTag>
struct LinearSolver { using type = UndefinedProperty; };

////////////////////////////////////////////////
// Basic properties regarding balance equations
/////////////////////////////////////////////////
// TODO: Integrate UseMoles into BalanceEqOpts
template<class TypeTag, class MyTypeTag>
struct UseMoles { using type = UndefinedProperty; };               //!< Property whether to use moles or kg as amount unit for balance equations
template<class TypeTag, class MyTypeTag>
struct ReplaceCompEqIdx { using type = UndefinedProperty; };       //!< The component balance index that should be replaced by the total mass/mole balance
template<class TypeTag, class MyTypeTag>
struct BalanceEqOpts { using type = UndefinedProperty; };          //!< A class that collects options for the evaluation of the balance equations

/////////////////////////////////////////////
// Properties used by finite volume schemes:
/////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct ElementBoundaryTypes { using type = UndefinedProperty; };                //!< Stores the boundary types on an element

#if defined(__clang__) && !defined(DONT_EMIT_CLANG_GRIDGEOMETRY_WARNING)
#warning "The properties `FVGridGeometry` and `EnableFVGridGeometryCache` \
are deprecated in favor of `GridGeometry` and `EnableGridGeometryCache`. \
The old properties will be removed after release 3.1. \
If clang is used, no deprecation warnings are emitted. \
We recommend to use gcc for getting rid of the warnings. \
You can suppress this message by defining the preprocessor variable \
DONT_EMIT_CLANG_GRIDGEOMETRY_WARNING."
#endif
// TODO: Remove deprecated property FVGridGeometry after 3.1
template<class TypeTag, class MyTypeTag>
struct [[deprecated("Use GridGeometry instead.")]] FVGridGeometry { using type = UndefinedProperty; }; //!< The type of the global finite volume geometry

// TODO: Remove deprecated property EnableFVGridGeometryCache after 3.1
template<class TypeTag, class MyTypeTag>
struct [[deprecated("Use EnableGridGeometryCache instead.")]] EnableFVGridGeometryCache { using type = UndefinedProperty; };           //!< specifies if geometric data is saved (faster, but more memory consuming)

// Dumux 3.1 changes the property `FVGridGeometry` to `GridGeometry`.
// For ensuring backward compatibility, it is necessary to set the default value
// of the new property to the old one, see the discussion in MR 1647.
// Use diagnostic pragmas to prevent the emission of a warning message.
// TODO after 3.1: change default vale to `UndefinedProperty`, remove pragmas
// and comment.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

template<class TypeTag, class T>
struct GridGeometryHelper
{ using type = GetPropType<TypeTag, Properties::FVGridGeometry>; };

template<class TypeTag>
struct GridGeometryHelper<TypeTag, UndefinedProperty>
{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct GridGeometry
{
    using type = typename GridGeometryHelper<TypeTag, typename FVGridGeometry<TypeTag, MyTypeTag>::type>::type;
};

template<class TypeTag, bool hasParentTypeTag>
struct EnableGridGeometryCacheHelper
{ using type = UndefinedProperty; };

template<class TypeTag>
struct EnableGridGeometryCacheHelper<TypeTag, false>
{
    // fallback
    static constexpr bool value = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
};

// only use the fallback (EnableFVGridGeometryCache) if none
// of the TypeTags define EnableGridGeometryCache
template <class TypeTag, class MyTypeTag>
struct EnableGridGeometryCache : public EnableGridGeometryCacheHelper<TypeTag, Detail::hasParentTypeTag<MyTypeTag>(int{})>
{};

#pragma GCC diagnostic pop

template<class TypeTag, class MyTypeTag>
struct VolumeVariables { using type = UndefinedProperty; };                     //!< The secondary variables within a sub-control volume
template<class TypeTag, class MyTypeTag>
struct GridVolumeVariables { using type = UndefinedProperty; };                 //!< The type for a global container for the volume variables
template<class TypeTag, class MyTypeTag>
struct EnableGridVolumeVariablesCache { using type = UndefinedProperty; };      //!< If disabled, the volume variables are not stored (reduces memory, but is slower)
template<class TypeTag, class MyTypeTag>
struct FluxVariables { using type = UndefinedProperty; };                       //!< Container storing the different types of flux variables
template<class TypeTag, class MyTypeTag>
struct FluxVariablesCache { using type = UndefinedProperty; };                  //!< Stores data associated with flux vars
template<class TypeTag, class MyTypeTag>
struct FluxVariablesCacheFiller { using type = UndefinedProperty; };            //!< The engine behind the global flux cache (how to fill caches for the stencil)
template<class TypeTag, class MyTypeTag>
struct GridFluxVariablesCache { using type = UndefinedProperty; };              //!< The global vector of flux variable containers
template<class TypeTag, class MyTypeTag>
struct EnableGridFluxVariablesCache { using type = UndefinedProperty; };        //!< specifies if data on flux vars should be saved (faster, but more memory consuming)
template<class TypeTag, class MyTypeTag>
struct GridVariables { using type = UndefinedProperty; };                       //!< The grid variables object managing variable data on the grid (volvars/fluxvars cache)

/////////////////////////////////////////////////////////////////
// Additional properties used by the cell-centered mpfa schemes:
/////////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct PrimaryInteractionVolume { using type = UndefinedProperty; };            //!< The primary interaction volume type
template<class TypeTag, class MyTypeTag>
struct SecondaryInteractionVolume { using type = UndefinedProperty; };          //!< The secondary interaction volume type used e.g. on the boundaries
template<class TypeTag, class MyTypeTag>
struct DualGridNodalIndexSet { using type = UndefinedProperty; };               //!< The type used for the nodal index sets of the dual grid

/////////////////////////////////////////////////////////////
// Properties used by models involving flow in porous media:
/////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct EnergyLocalResidual { using type = UndefinedProperty; };                 //!< The local residual of the energy equation
template<class TypeTag, class MyTypeTag>
struct AdvectionType { using type = UndefinedProperty; };                       //!< The type for the calculation the advective fluxes
template<class TypeTag, class MyTypeTag>
struct SolutionDependentAdvection { using type = UndefinedProperty; };          //!< specifies if the parameters for the advective fluxes depend on the solution
template<class TypeTag, class MyTypeTag>
struct MolecularDiffusionType { using type = UndefinedProperty; };              //!< The type for the calculation of the molecular diffusion fluxes
template<class TypeTag, class MyTypeTag>
struct SolutionDependentMolecularDiffusion { using type = UndefinedProperty; }; //!< specifies if the parameters for the diffusive fluxes depend on the solution
template<class TypeTag, class MyTypeTag>
struct HeatConductionType { using type = UndefinedProperty; };                  //!< The type for the calculation of the heat conduction fluxes
template<class TypeTag, class MyTypeTag>
struct SolutionDependentHeatConduction { using type = UndefinedProperty; };     //!< specifies if the parameters for the heat conduction fluxes depend on the solution

template<class TypeTag, class MyTypeTag>
struct SpatialParams { using type = UndefinedProperty; };                       //!< The type of the spatial parameters object
template<class TypeTag, class MyTypeTag>
struct FluidSystem { using type = UndefinedProperty; };                         //!< The type of the fluid system to use
template<class TypeTag, class MyTypeTag>
struct FluidState { using type = UndefinedProperty; };                          //!< The type of the fluid state to use
template<class TypeTag, class MyTypeTag>
struct SolidSystem { using type = UndefinedProperty; };                         //!< The type of the solid system to use
template<class TypeTag, class MyTypeTag>
struct SolidState { using type = UndefinedProperty; };                           //!< The type of the solid state to use
template<class TypeTag, class MyTypeTag>
struct EffectiveDiffusivityModel { using type = UndefinedProperty; };           //!< The employed model for the computation of the effective diffusivity
template<class TypeTag, class MyTypeTag>
struct ThermalConductivityModel { using type = UndefinedProperty; };            //!< Model to be used for the calculation of the effective conductivity
template<class TypeTag, class MyTypeTag>
struct VelocityOutput { using type = UndefinedProperty; };                      //!< specifies the velocity calculation module to be used
template<class TypeTag, class MyTypeTag>
struct Formulation { using type = UndefinedProperty; };                         //!< The formulation of the model
// TODO: is this useful? -> everything is a constraint solver just a different type
template<class TypeTag, class MyTypeTag>
struct UseConstraintSolver { using type = UndefinedProperty; };                 //!< Whether to use a contraint solver for computing the secondary variables

// When using the box method in a multi-phase context, an interface solver might be necessary
template<class TypeTag, class MyTypeTag>
struct EnableBoxInterfaceSolver { using type = UndefinedProperty; };

//////////////////////////////////////////////////////////////
// Additional properties used by the 2pnc and 2pncmin models:
//////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct Chemistry { using type = UndefinedProperty; };                           //!< The chemistry class with which solves equlibrium reactions
template<class TypeTag, class MyTypeTag>
struct SetMoleFractionsForFirstPhase { using type = UndefinedProperty; };       //!< Set the mole fraction in the wetting or non-wetting phase

//////////////////////////////////////////////////////////////
// Additional properties used by the richards model
//////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct EnableWaterDiffusionInAir { using type = UndefinedProperty; }; //!< Property for turning Richards into extended Richards

//////////////////////////////////////////////////////////////
// Additional properties used by the 3pwateroil model:
//////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct OnlyGasPhaseCanDisappear { using type = UndefinedProperty; }; //!< reduces the phasestates to threePhases and wnPhaseOnly

/////////////////////////////////////////////////////////////
// Properties used by geomechanical models:
/////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct StressType { using type = UndefinedProperty; };       //!< The type used for the evaluation of stress tensors and forces

/////////////////////////////////////////////////////////////
// Properties used by the staggered-grid discretization method
/////////////////////////////////////////////////////////////

template<class TypeTag, class MyTypeTag>
struct NumEqCellCenter { using type = UndefinedProperty; };                     //!< The number of equations for cell-centered dofs
template<class TypeTag, class MyTypeTag>
struct NumEqFace { using type = UndefinedProperty; };                           //!< The number of equations for face dofs
template<class TypeTag, class MyTypeTag>
struct CellCenterSolutionVector { using type = UndefinedProperty; };            //!< The solution vector type for cell-centered dofs
template<class TypeTag, class MyTypeTag>
struct FaceSolutionVector { using type = UndefinedProperty; };                  //!< The solution vector type for face dofs
template<class TypeTag, class MyTypeTag>
struct GridFaceVariables { using type = UndefinedProperty; };                   //!< Global vector containing face-related data
template<class TypeTag, class MyTypeTag>
struct CellCenterPrimaryVariables { using type = UndefinedProperty; };          //!< The primary variables container type for cell-centered dofs
template<class TypeTag, class MyTypeTag>
struct FacePrimaryVariables { using type = UndefinedProperty; };                //!< The primary variables container type for face dofs
template<class TypeTag, class MyTypeTag>
struct IntersectionMapper { using type = UndefinedProperty; };                  //!< Specifies the intersection mapper
template<class TypeTag, class MyTypeTag>
struct StaggeredPrimaryVariables { using type = UndefinedProperty; };           //!< The hybrid primary variables container type
template<class TypeTag, class MyTypeTag>
struct BaseEpsilon { using type = UndefinedProperty; };                         //!< A base epsilon for numerical differentiation, can contain multiple values
template<class TypeTag, class MyTypeTag>
struct FaceVariables { using type = UndefinedProperty; };                       //!< Class containing local face-related data
template<class TypeTag, class MyTypeTag>
struct BoundaryValues { using type = UndefinedProperty; };                      //!< Class containing local boundary data
template<class TypeTag, class MyTypeTag>
struct StaggeredFaceSolution { using type = UndefinedProperty; };               //!< A vector containing the solution for a face (similar to ElementSolution)
template<class TypeTag, class MyTypeTag>
struct EnableGridFaceVariablesCache { using type = UndefinedProperty; };        //!< Switch on/off caching of face variables
template<class TypeTag, class MyTypeTag>
struct UpwindSchemeOrder { using type = UndefinedProperty; };                   //!< Specifies the order of the upwinding scheme (1 == first order, 2 == second order(tvd methods))

/////////////////////////////////////////////////////////////
// Properties used by the mpnc model
/////////////////////////////////////////////////////////////

template<class TypeTag, class MyTypeTag>
struct PressureFormulation { using type = UndefinedProperty; }; //! the formulation of the pressure e.g most wetting first

/////////////////////////////////////////////////////////////
// Properties used by the nonequilibrium model
/////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct EquilibriumModelTraits { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct EquilibriumLocalResidual { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct EquilibriumIndices { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct EquilibriumIOFields { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct NumEqBalance { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct EnableThermalNonEquilibrium { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct EnableChemicalNonEquilibrium { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct NumEnergyEqFluid { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct NumEnergyEqSolid { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct NusseltFormulation { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct SherwoodFormulation { using type = UndefinedProperty; };

/////////////////////////////////////////////////////////////
// Properties used by free flow models
/////////////////////////////////////////////////////////////

template<class TypeTag, class MyTypeTag>
struct NormalizePressure { using type = UndefinedProperty; }; //!<  Returns whether to normalize the pressure term in the momentum balance or not

/////////////////////////////////////////////////////////////
// Properties used by multidomain simulations
/////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct CouplingManager { using type = UndefinedProperty; };

///////////////////////////////////////
// Basic properties of sequential models:
///////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct TimeManager { using type = UndefinedProperty; };

} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_PROPERTIES_HH
