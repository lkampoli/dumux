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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableInertiaTerms);
NEW_PROP_TAG(ReplaceCompEqIdx);
}

/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */



// forward declaration
template<class TypeTag, bool enableComponentTransport>
class StaggeredNavierStokesResidualImpl;

template<class TypeTag>
using StaggeredNavierStokesResidual = StaggeredNavierStokesResidualImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport)>;


template<class TypeTag>
class StaggeredNavierStokesResidualImpl<TypeTag, false> : public Dumux::StaggeredLocalResidual<TypeTag>
{
    using ParentType = StaggeredLocalResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using EnergyFluxVariables = typename GET_PROP_TYPE(TypeTag, EnergyFluxVariables);


    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CellCenterResidual = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidual = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using FaceResidualVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
    static constexpr bool normalizePressure = GET_PROP_VALUE(TypeTag, NormalizePressure);

public:

    using ParentType::ParentType;


    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const GlobalFaceVars& globalFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        CellCenterPrimaryVariables flux = fluxVars.computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars,
                                                 globalFaceVars, scvf, elemFluxVarsCache[scvf]);

        EnergyFluxVariables::energyFlux(flux, problem, element, fvGeometry, elemVolVars, globalFaceVars, scvf, elemFluxVarsCache[scvf]);

        return flux;
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Problem& problem,
                                                          const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const GlobalFaceVars& globalFaceVars,
                                                          const SubControlVolume &scv) const
    {
        return CellCenterPrimaryVariables(0.0);
        // TODO sources
    }


     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv,
                                    const VolumeVariables& volVars)
    {
        CellCenterPrimaryVariables storage;
        storage[0] = volVars.density();
        EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars);
        return storage;
    }

     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scvf The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    FacePrimaryVariables computeStorageForFace(const SubControlVolumeFace& scvf,
                                               const VolumeVariables& volVars,
                                               const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables storage(0.0);
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();
        storage[0] = volVars.density() * velocity;
        return storage;
    }

    FacePrimaryVariables computeSourceForFace(const Problem& problem,
                                              const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const GlobalFaceVars& globalFaceVars) const
    {
        FacePrimaryVariables source(0.0);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        source += problem.gravity()[scvf.directionIndex()] * insideVolVars.density();

        source += problem.sourceAtPos(scvf.center())[faceIdx][scvf.directionIndex()];

        return source;
    }

     /*!
     * \brief Returns the complete momentum flux for a face
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computeFluxForFace(const Problem& problem,
                                            const Element& element,
                                            const SubControlVolumeFace& scvf,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const GlobalFaceVars& globalFaceVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FacePrimaryVariables flux(0.0);
        FluxVariables fluxVars;
        flux += fluxVars.computeNormalMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, globalFaceVars);
        flux += fluxVars.computeTangetialMomentumFlux(problem, element, scvf, fvGeometry, elemVolVars, globalFaceVars);
        flux += computePressureTerm_(problem, element, scvf, fvGeometry, elemVolVars, globalFaceVars);
        return flux;
    }

protected:

     /*!
     * \brief Evaluate boundary conditions
     */
    void evalBoundary_(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const GlobalFaceVars& faceVars,
                       const ElementBoundaryTypes& elemBcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        evalBoundaryForCellCenter_(element, fvGeometry, elemVolVars, faceVars, elemBcTypes, elemFluxVarsCache);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            evalBoundaryForFace_(element, fvGeometry, scvf, elemVolVars, faceVars, elemBcTypes, elemFluxVarsCache);
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a cell center dof
     */
    void evalBoundaryForCellCenter_(CellCenterResidual& residual,
                                    const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const GlobalFaceVars& faceVars,
                                    const ElementBoundaryTypes& elemBcTypes,
                                    const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                auto boundaryFlux = computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, faceVars, scvf, elemFluxVarsCache);

                // handle the actual boundary conditions:
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                if(bcTypes.hasNeumann())
                {
                    // handle Neumann BCs, i.e. overwrite certain fluxes by user-specified values
                    for(int eqIdx = 0; eqIdx < GET_PROP_VALUE(TypeTag, NumEqCellCenter); ++eqIdx)
                        if(bcTypes.isNeumann(eqIdx))
                        {
                            const auto extrusionFactor = 1.0; //TODO: get correct extrusion factor
                            boundaryFlux[eqIdx] = problem.neumann(element, scvf)[cellCenterIdx][eqIdx]
                                                   * extrusionFactor * scvf.area();
                        }
                }

                residual += boundaryFlux;

                asImp_().setFixedCell_(residual, problem, fvGeometry.scv(scvf.insideScvIdx()), elemVolVars, bcTypes);
            }
        }
    }

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     *
     * \param insideScv The sub control volume
     * \param elemVolVars The current or previous element volVars
     * \param bcTypes The boundary types
     */
    void setFixedCell_(CellCenterResidual& residual,
                       const Problem& problem,
                       const SubControlVolume& insideScv,
                       const ElementVolumeVariables& elemVolVars,
                       const BoundaryTypes& bcTypes) const
    {
        // set a fixed pressure for cells adjacent to a wall
        if(bcTypes.isDirichletCell(massBalanceIdx))
        {
            const auto& insideVolVars = elemVolVars[insideScv];
            residual[pressureIdx] = insideVolVars.pressure() - problem.dirichletAtPos(insideScv.center())[cellCenterIdx][pressureIdx];
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a face dof
     */
    void evalBoundaryForFace_(FaceResidual& residual,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf,
                              const ElementVolumeVariables& elemVolVars,
                              const GlobalFaceVars& faceVars,
                              const ElementBoundaryTypes& elemBcTypes,
                              const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if (scvf.boundary())
        {
            // handle the actual boundary conditions:
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            // set a fixed value for the velocity for Dirichlet boundary conditions
            if(bcTypes.isDirichlet(momentumBalanceIdx))
            {
                const Scalar velocity = faceVars.faceVars(scvf.dofIndex()).velocity();
                const Scalar dirichletValue = problem.dirichlet(element, scvf)[faceIdx][scvf.directionIndex()];
                residual = velocity - dirichletValue;
            }

            // For symmetry boundary conditions, there is no flow accross the boundary and
            // we therefore treat it like a Dirichlet boundary conditions with zero velocity
            if(bcTypes.isSymmetry())
            {
                const Scalar velocity = faceVars.faceVars(scvf.dofIndex()).velocity();
                const Scalar fixedValue = 0.0;
                residual = velocity - fixedValue;
            }

            // outflow condition for the momentum balance equation
            if(bcTypes.isOutflow(momentumBalanceIdx))
            {
                if(bcTypes.isDirichlet(massBalanceIdx))
                    residual += computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, faceVars, elemFluxVarsCache);
                else
                    DUNE_THROW(Dune::InvalidStateException, "Face at " << scvf.center()  << " has an outflow BC for the momentum balance but no Dirichlet BC for the pressure!");
            }
        }
//        std::cout << "** staggered/localresidual: massBalanceIdx = " << massBalanceIdx
//        		<< ", momBalanceIdx = " << momentumBalanceIdx << std::endl;
    }


private:

     /*!
     * \brief Returns the pressure term
     * \param scvf The sub control volume face
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     */
    FacePrimaryVariables computePressureTerm_(const Problem& problem,
                                              const Element& element,
                                              const SubControlVolumeFace& scvf,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const GlobalFaceVars& globalFaceVars) const
    {
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        const Scalar deltaP = normalizePressure ? problem.initialAtPos(scvf.center())[cellCenterIdx][pressureIdx] : 0.0;

        Scalar result = (insideVolVars.pressure() - deltaP) * scvf.area() * -1.0 * sign(scvf.outerNormalScalar());

        // treat outflow BCs
        if(scvf.boundary())
        {
            const Scalar pressure = problem.dirichlet(element, scvf)[cellCenterIdx][pressureIdx] - deltaP;
            result += pressure * scvf.area() * sign(scvf.outerNormalScalar());
        }
        return result;
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};
}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
