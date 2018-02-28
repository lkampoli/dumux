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
#ifndef DUMUX_WELL_MIMETIC_LOCAL_RESIDUAL_HH
#define DUMUX_WELL_MIMETIC_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

namespace Dumux
{

/*!
 * \ingroup CCModel
 * \ingroup WellImmiscibleMimeticLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */


template<class TypeTag>
class WellImmiscibleMimeticLocalResidual : public Dumux::StaggeredLocalResidual<TypeTag>
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

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

    enum { conti0EqIdx = Indices::conti0EqIdx};

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:
    // copying the local residual class is not a good idea
    WellImmiscibleMimeticLocalResidual(const WellImmiscibleMimeticLocalResidual &) = delete;

    WellImmiscibleMimeticLocalResidual() = default;


    CellCenterPrimaryVariables computeFluxForCellCenter(const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& globalFaceVars,
                                  const SubControlVolumeFace &scvf,
                                  const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        FluxVariables fluxVars;
        fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, globalFaceVars, scvf, elemFluxVarsCache);

        PrimaryVariables flux(0.0);

        auto globalPos = element.geometry().center();
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        if(!(std::abs(x-7.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8) &&
           !(std::abs(x-15.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8))
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                // the physical quantities for which we perform upwinding
                auto upwindTerm = [phaseIdx](const VolumeVariables& volVars)
                                  { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };

                auto eqIdx = conti0EqIdx + phaseIdx;
                flux[eqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
                EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
            }

            //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
        }

        return flux[cellCenterIdx];
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Element &element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const GlobalFaceVars& globalFaceVars,
                                     const SubControlVolume &scv)
    {
        CellCenterPrimaryVariables source(0.0);

        // add contributions from volume flux sources
        source += this->problem().source(element, fvGeometry, elemVolVars, scv)[cellCenterIdx];

        // add contribution from possible point sources
        //source += this->problem().scvPointSources(element, fvGeometry, elemVolVars, scv);

        return source;
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
        // partial time derivative of the phase mass
        PrimaryVariables storage(0.0);

        auto globalPos = scv.center();
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        if(!(std::abs(x-7.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8) &&
           !(std::abs(x-15.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8))
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                auto eqIdx = conti0EqIdx + phaseIdx;
                storage[eqIdx] = volVars.porosity()
                                 * volVars.density(phaseIdx)
                                 * volVars.saturation(phaseIdx);

                //! The energy storage in the fluid phase with index phaseIdx
                EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
            }

            //! The energy storage in the solid matrix
            EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);
        }
        return storage[cellCenterIdx];
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
        return storage;
    }

    FacePrimaryVariables computeSourceForFace(const SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables zero(0.0);
        return zero;
    }


    FacePrimaryVariables computeFluxForFace(const Element& element,
                                            const SubControlVolumeFace& scvf,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const GlobalFaceVars& globalFaceVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        FluxVariables fluxVars;
        fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, globalFaceVars, scvf, elemFluxVarsCache);

        PrimaryVariables flux(0.0);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // the physical quantities for which we perform upwinding
            auto upwindTerm = [phaseIdx](const VolumeVariables& volVars)
                              { return volVars.density(phaseIdx)/volVars.viscosity(phaseIdx); };

            auto eqIdx = conti0EqIdx + phaseIdx;
            flux[eqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
        }

        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux[cellCenterIdx];
    }

protected:
     /*!
     * \brief Evaluate boundary conditions for a cell center dof
     */
    void evalBoundaryForCellCenter_(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const GlobalFaceVars& globalFaceVars,
                                    const ElementBoundaryTypes& elemBcTypes,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                // handle the actual boundary conditions:
                const auto bcTypes = this->problem().boundaryTypes(element, scvf);

                // set a fixed value for the velocity
                if(bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
                {
                    auto neumannFluxes = this->problem().neumann(element, fvGeometry, elemVolVars, scvf)[cellCenterIdx];

                    // multiply neumann fluxes with the area and the extrusion factor
                    auto&& scv = fvGeometry.scv(scvf.insideScvIdx());
                    neumannFluxes *= scvf.area()*elemVolVars[scv].extrusionFactor();

                    this->ccResidual_ += neumannFluxes;
                }
                else if(!bcTypes.hasNeumann() && bcTypes.hasDirichlet())
                {
                    this->ccResidual_ += computeFluxForCellCenter(element, fvGeometry, elemVolVars, globalFaceVars, scvf, elemFluxVarsCache);
                }
            }
        }
    }

     /*!
     * \brief Evaluate boundary conditions for a face dof
     */
    void evalBoundaryForFace_(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf,
                              const ElementVolumeVariables& elemVolVars,
                              const GlobalFaceVars& globalFaceVars,
                              const ElementBoundaryTypes& elemBcTypes,
                              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (scvf.boundary())
        {
            // handle the actual boundary conditions:
            const auto bcTypes = this->problem().boundaryTypes(element, scvf);

            // set a fixed value for the velocity
            if(bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVolVars = elemVolVars[insideScv];

                auto neumannFlux = this->evalNeumannSegment_(element, fvGeometry, elemVolVars, scvf, bcTypes);

                FluxVariables fluxVars;
                fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, globalFaceVars, scvf, elemFluxVarsCache);

//                Scalar pdiff = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[nPhaseIdx] - globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[wPhaseIdx];
//                Scalar Sw = MaterialLaw::sw(this->problem().spatialParams().materialLawParamsAtPos(scvf.center()), pdiff);
//                std::vector<Scalar> kr(numPhases);
//                kr[wPhaseIdx] =  MaterialLaw::krw(this->problem().spatialParams().materialLawParamsAtPos(scvf.center()), Sw);
//                kr[nPhaseIdx] =  MaterialLaw::krn(this->problem().spatialParams().materialLawParamsAtPos(scvf.center()), Sw);

                PrimaryVariables flux;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    auto upwindTerm = [phaseIdx](const VolumeVariables& volVars)
                                      { return volVars.density(phaseIdx)/volVars.viscosity(phaseIdx); };

                    auto eqIdx = conti0EqIdx + phaseIdx;
                    flux[eqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                }

                EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

                this->faceResiduals_[scvf.localFaceIdx()] = flux[cellCenterIdx];
                this->faceResiduals_[scvf.localFaceIdx()] -= neumannFlux;

//                this->faceResiduals_[scvf.localFaceIdx()] = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()
//                                                          - insideVolVars.pressure(0);

//               this->faceResiduals_[scvf.localFaceIdx()] = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars();

            }
            else if(!bcTypes.hasNeumann() && bcTypes.hasDirichlet())
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVolVars = elemVolVars[insideScv];
                this->faceResiduals_[scvf.localFaceIdx()] = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()
                                                          - this->problem().dirichlet(element, scvf)[faceIdx];
            }
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single scv face
     */
    FacePrimaryVariables evalNeumannSegment_(const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace &scvf,
                                         const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the neumann boundary fluxes
        FacePrimaryVariables flux(0);

        auto neumannFluxes = this->problem().neumann(element, fvGeometry, elemVolVars, scvf)[faceIdx];

        // multiply neumann fluxes with the area and the extrusion factor
        auto&& scv = fvGeometry.scv(scvf.insideScvIdx());
        neumannFluxes *= scvf.area()*elemVolVars[scv].extrusionFactor();

        // add fluxes to the residual
//        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//            if (bcTypes.isNeumann(eqIdx))
//                flux[eqIdx] += neumannFluxes[faceIdx][eqIdx];

        return neumannFluxes;
    }


};
}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH