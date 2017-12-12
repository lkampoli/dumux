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
 * \file
 *
 * \brief _Declares_ all properties used in Dumux.
 * \note Include this to forward declare properties in your headers.
 */

#ifndef DUMUX_PROPERTIES_HH
#define DUMUX_PROPERTIES_HH

#ifndef DUMUX_PROPERTY_SYSTEM_HH
#include <dumux/common/properties/propertysystem.hh>
#endif

namespace Dumux
{
namespace Properties
{
///////////////////////////////////////
// Basic properties of numeric models:
///////////////////////////////////////
NEW_PROP_TAG(Scalar);                 //! Property to specify the type of scalar values.
NEW_PROP_TAG(ModelParameterGroup);    //! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(ModelDefaultParameters); //! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(GridCreator);            //! Property which provides a GridCreator (manages grids)
NEW_PROP_TAG(Grid);                   //! The DUNE grid type
NEW_PROP_TAG(NumEq);                  //! The number of equations to solve (equal to number of primary variables)
NEW_PROP_TAG(Indices);                //! Enumerations for the numeric model
NEW_PROP_TAG(PrimaryVariables);       //! A vector of primary variables
NEW_PROP_TAG(NumEqVector);            //! A vector of size number equations that can be used for Neumann fluxes, sources, residuals, ...
NEW_PROP_TAG(GridView);               //! The type of the grid view according to the grid type
NEW_PROP_TAG(Problem);                //! Property to specify the type of a problem which has to be solved
NEW_PROP_TAG(PointSource);            //! Property defining the type of point source used
NEW_PROP_TAG(PointSourceHelper);      //! Property defining the class that computes which sub control volume point sources belong to
NEW_PROP_TAG(VtkOutputFields);        //! A class helping models to define default vtk output parameters
NEW_PROP_TAG(BaseLocalResidual);      //! The type of the base class of the local residual (specific to a discretization scheme)
NEW_PROP_TAG(JacobianMatrix);         //! Type of the global jacobian matrix
NEW_PROP_TAG(SolutionVector);         //! Vector containing all primary variable vector of the grid
NEW_PROP_TAG(BoundaryTypes);          //! Stores the boundary types of a single degree of freedom
NEW_PROP_TAG(DiscretizationMethod);   //! Property for the used discretization method
NEW_PROP_TAG(VertexMapper);           //! mapper for vertices
NEW_PROP_TAG(ElementMapper);          //! mapper for elements

//! The type of the local residual function, i.e. the equation to be solved. Must inherit
//! from the BaseLocalResidual property and fulfill its interfaces.
NEW_PROP_TAG(LocalResidual);

//! TODO: Remove this property as soon as the decoupled models are integrated
NEW_PROP_TAG(LinearSolver);
NEW_PROP_TAG(LinearSolverPreconditionerBlockLevel); //! Block level depth for the preconditioner

////////////////////////////////////////////////
// Basic properties regarding balance equations
/////////////////////////////////////////////////
// TODO: Integrate UseMoles into BalanceEqOpts
NEW_PROP_TAG(UseMoles);               //! Property whether to use moles or kg as amount unit for balance equations
NEW_PROP_TAG(ReplaceCompEqIdx);       //! The component balance index that should be replaced by the total mass/mole balance
NEW_PROP_TAG(BalanceEqOpts);          //! A class that collects options for the evaluation of the balance equations

/////////////////////////////////////////////
// Properties used by finite volume schemes:
/////////////////////////////////////////////
NEW_PROP_TAG(ElementBoundaryTypes);                //! Stores the boundary types on an element
NEW_PROP_TAG(ElementSolutionVector);               //! A vector of primary variables within an element
NEW_PROP_TAG(AssemblyMap);                         //! Connectivity map (transposed) used for assembling the Jacobian matrix entries

NEW_PROP_TAG(SubControlVolume);                    //! The type of the sub control volume
NEW_PROP_TAG(SubControlVolumeFace);                //! The type of the sub control volume face
NEW_PROP_TAG(FVElementGeometry);                   //! The type of the local finite volume geometry (iterators over scvs, scvfs)
NEW_PROP_TAG(FVGridGeometry);                      //! The type of the global finite volume geometry
NEW_PROP_TAG(EnableFVGridGeometryCache);           //! specifies if geometric data is saved (faster, but more memory consuming)

NEW_PROP_TAG(VolumeVariables);                     //! The secondary variables within a sub-control volume
NEW_PROP_TAG(ElementVolumeVariables);              //! The type for a local (element/stencil) container for the volume variables
NEW_PROP_TAG(GlobalVolumeVariables);               //! The type for a global container for the volume variables
NEW_PROP_TAG(EnableGlobalVolumeVariablesCache);    //! If disabled, the volume variables are not stored (reduces memory, but is slower)
NEW_PROP_TAG(FluxVariables);                       //! Container storing the different types of flux variables
NEW_PROP_TAG(FluxVariablesCache);                  //! Stores data associated with flux vars
NEW_PROP_TAG(ElementFluxVariablesCache);           //! A local vector of flux variable caches per element
NEW_PROP_TAG(GlobalFluxVariablesCache);            //! The global vector of flux variable containers
NEW_PROP_TAG(EnableGlobalFluxVariablesCache);      //! specifies if data on flux vars should be saved (faster, but more memory consuming)
NEW_PROP_TAG(GridVariables);                       //! The grid variables object managing variable data on the grid (volvars/fluxvars cache)

/////////////////////////////////////////////////////////////////
// Additional properties used by the cell-centered mpfa schemes:
/////////////////////////////////////////////////////////////////
NEW_PROP_TAG(MpfaMethod);                          //! Specifies the mpfa method to be used
NEW_PROP_TAG(MpfaHelper);                          //! A Helper class depending on the mpfa method and grid dimension
NEW_PROP_TAG(PrimaryInteractionVolume);            //! The primary interaction volume type
NEW_PROP_TAG(SecondaryInteractionVolume);          //! The secondary interaction volume type used e.g. on the boundaries


/////////////////////////////////////////////////////////////
// Properties used by models involving flow in porous media:
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(EnergyLocalResidual);                 //! The local residual of the energy equation
NEW_PROP_TAG(EnableAdvection);                     //! specifies if advection is considered in the model
NEW_PROP_TAG(AdvectionType);                       //! The type for the calculation the advective fluxes
NEW_PROP_TAG(SolutionDependentAdvection);          //! specifies if the parameters for the advective fluxes depend on the solution
NEW_PROP_TAG(EnableMolecularDiffusion);            //! specifies if molecular diffusive fluxes are considered in the model
NEW_PROP_TAG(MolecularDiffusionType);              //! The type for the calculation of the molecular diffusion fluxes
NEW_PROP_TAG(SolutionDependentMolecularDiffusion); //! specifies if the parameters for the diffusive fluxes depend on the solution
NEW_PROP_TAG(EnableEnergyBalance);                 //! Specifies if the model solves an energy equation
NEW_PROP_TAG(HeatConductionType);                  //! The type for the calculation of the heat conduction fluxes
NEW_PROP_TAG(SolutionDependentHeatConduction);     //! specifies if the parameters for the heat conduction fluxes depend on the solution

NEW_PROP_TAG(NumPhases);                           //! Number of fluid phases in the system
NEW_PROP_TAG(PhaseIdx);                            //! A phase index to allow using a two-phase fluidsystem for one-phase models
NEW_PROP_TAG(NumComponents);                       //! Number of fluid phases in the system
NEW_PROP_TAG(SpatialParams);                       //! The type of the spatial parameters object
NEW_PROP_TAG(FluidSystem);                         //! The type of the fluid system to use
NEW_PROP_TAG(FluidState);                          //! The type of the fluid state to use
NEW_PROP_TAG(PrimaryVariableSwitch);               //! The primary variable switch needed for compositional models
NEW_PROP_TAG(EffectiveDiffusivityModel);           //! The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(ThermalConductivityModel);            //! Model to be used for the calculation of the effective conductivity
NEW_PROP_TAG(VelocityOutput);                      //! specifies the velocity calculation module to be used

NEW_PROP_TAG(MaterialLaw);                         //! The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(Formulation);                         //! The formulation of the model
// TODO: is this useful? -> everything is a constraint solver just a different type
NEW_PROP_TAG(UseConstraintSolver);                 //! Whether to use a contraint solver for computing the secondary variables
NEW_PROP_TAG(UseKelvinEquation);                   //! If we use Kelvin equation to lower the vapor pressure as a function of capillary pressure, temperature

////////////////////////////////////////////////////////////////////////////////
// Properties used by models involving mineralization:
////////////////////////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumSPhases);
NEW_PROP_TAG(NonMineralizationVtkOutputFields);
NEW_PROP_TAG(NonMineralizationVolumeVariables);

NEW_PROP_TAG(UseConstraintSolver);                 //! Determines whether the constraint solver should be used#

/////////////////////////////////////////////////////////////
// non-isothermal porous medium flow models
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(IsothermalVtkOutputFields);
NEW_PROP_TAG(IsothermalVolumeVariables);
NEW_PROP_TAG(IsothermalLocalResidual);
NEW_PROP_TAG(IsothermalIndices);
NEW_PROP_TAG(IsothermalNumEq);

// specify if we evaluate the permeability in the volume (for discontinuous fields)
// or at the scvf center for analytical permeability fields (e.g. convergence studies)
NEW_PROP_TAG(EvaluatePermeabilityAtScvfIP);

//////////////////////////////////////////////////////////////
// Additional properties used by the 2pnc and 2pncmin models:
//////////////////////////////////////////////////////////////
NEW_PROP_TAG(Chemistry);                           //!< The chemistry class with which solves equlibrium reactions
NEW_PROP_TAG(NumMajorComponents);                  //!< Number of major fluid components which are considered in the calculation of the phase density
NEW_PROP_TAG(SetMoleFractionsForWettingPhase);     //!< Set the mole fraction in the wetting or non-wetting phase

//////////////////////////////////////////////////////////////
// Additional properties used by the richards model
//////////////////////////////////////////////////////////////
NEW_PROP_TAG(EnableWaterDiffusionInAir); //!< Property for turning Richards into extended Richards

/////////////////////////////////////////////////////////////
// Properties used by the staggered-grid discretization method
/////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumEqCellCenter);                     //! The number of equations for cell-centered dofs
NEW_PROP_TAG(NumEqFace);                           //! The number of equations for face dofs
NEW_PROP_TAG(CellCenterSolutionVector);            //! The solution vector type for cell-centered dofs
NEW_PROP_TAG(FaceSolutionVector);                  //! The solution vector type for face dofs
NEW_PROP_TAG(GlobalFaceVars);                      //! Class containing face-related data
NEW_PROP_TAG(CellCenterPrimaryVariables);          //! The primary variables container type for cell-centered dofs
NEW_PROP_TAG(FacePrimaryVariables);                //! The primary variables container type for face dofs
NEW_PROP_TAG(IntersectionMapper);                  //! Specifies the intersection mapper
NEW_PROP_TAG(DofTypeIndices);                      //! Specifies index types for accessing the multi type block vectors/matrices
NEW_PROP_TAG(StaggeredGeometryHelper);             //! Specifies a helper class for the staggered grid geometry
NEW_PROP_TAG(StaggeredPrimaryVariables);           //! The hybrid primary variables container type
NEW_PROP_TAG(BaseEpsilon);                         //! A base epsilon for numerical differentiation, can contain multiple values
NEW_PROP_TAG(FaceVariables);                       //! Class containing local face-related data
NEW_PROP_TAG(BoundaryValues);                      //! Class containing local boundary data

} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_PROPERTIES_HH