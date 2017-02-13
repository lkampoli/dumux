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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_FREELOW_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
}

// forward declaration
template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
class FreeFlowFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using FreeFlowFluxVariables = FreeFlowFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport),
                                                                 GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
// specialization for immiscible, isothermal flow
template<class TypeTag>
class FreeFlowFluxVariablesImpl<TypeTag, false, false>
: public FluxVariablesBase<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const GlobalFaceVars& globalFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        CellCenterPrimaryVariables flux(0.0);
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();

        const bool insideIsUpstream = sign(scvf.outerNormalScalar()) == sign(velocity) ? true : false;
        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;

        const Scalar upWindWeight = GET_PROP_VALUE(TypeTag, ImplicitUpwindWeight);

        flux = (upWindWeight * upstreamVolVars.density() +
               (1.0 - upWindWeight) * downstreamVolVars.density()) * velocity;

        return flux * scvf.area() * sign(scvf.outerNormalScalar());
    }

    void computeCellCenterToCellCenterStencil(Stencil& stencil,
                                              const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf)
    {
        // the first entry is always the cc dofIdx itself
        if(stencil.empty())
            stencil.push_back(scvf.insideScvIdx());
        if(!scvf.boundary())
            stencil.push_back(scvf.outsideScvIdx());
    }

    void computeCellCenterToFaceStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    void computeFaceToCellCenterStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        const int eIdx = scvf.insideScvIdx();
        stencil.push_back(scvf.insideScvIdx());

        for(const auto& data : scvf.pairData())
        {
            auto& normalFace = fvGeometry.scvf(eIdx, data.localNormalFaceIdx);
            const auto outerParallelElementDofIdx = normalFace.outsideScvIdx();
            if(!normalFace.boundary())
                stencil.push_back(outerParallelElementDofIdx);
        }
    }

    void computeFaceToFaceStencil(Stencil& stencil,
                                  const Problem& problem,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        // the first entries are always the face dofIdx itself and the one of the opposing face
        if(stencil.empty())
        {
            stencil.push_back(scvf.dofIndex());
            stencil.push_back(scvf.dofIndexOpposingFace());
        }

        for(const auto& data : scvf.pairData())
        {
            stencil.push_back(data.normalPair.first);
            const auto outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
            if(outerParallelFaceDofIdx >= 0)
                stencil.push_back(outerParallelFaceDofIdx);
            if(!scvf.boundary())
                stencil.push_back(data.normalPair.second);
        }
    }

    /*!
    * \brief Returns the normal part of the momentum flux
    * \param scvf The sub control volume face
    * \param fvGeometry The finite-volume geometry
    * \param elemVolVars All volume variables for the element
    * \param globalFaceVars The face variables
    */
   FacePrimaryVariables computeNormalMomentumFlux(const Problem& problem,
                                                  const SubControlVolumeFace& scvf,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const GlobalFaceVars& globalFaceVars)
   {
       const auto insideScvIdx = scvf.insideScvIdx();
       const auto& insideVolVars = elemVolVars[insideScvIdx];
       const Scalar velocitySelf = globalFaceVars.faceVars(scvf.dofIndex()).velocity() ;
       const Scalar velocityOpposite = globalFaceVars.faceVars(scvf.dofIndexOpposingFace()).velocity();
       FacePrimaryVariables normalFlux(0.0);

       if(navierStokes)
       {
           // advective part
           const Scalar vAvg = (velocitySelf + velocityOpposite) * 0.5;
           const Scalar vUp = (sign(scvf.outerNormalScalar()) == sign(vAvg)) ? velocityOpposite : velocitySelf;
           normalFlux += vAvg * vUp * insideVolVars.density();
       }

       // diffusive part
       const Scalar deltaV = scvf.normalInPosCoordDir() ?
                             (velocitySelf - velocityOpposite) :
                             (velocityOpposite - velocitySelf);

       const Scalar deltaX = scvf.selfToOppositeDistance();
       normalFlux -= insideVolVars.viscosity() * 2.0 * deltaV/deltaX;

       // account for the orientation of the face
       const Scalar sgn = -1.0 * sign(scvf.outerNormalScalar());

       Scalar result = normalFlux * sgn * scvf.area();

       // treat outflow conditions
       if(navierStokes && scvf.boundary())
       {
           const auto& upVolVars = (sign(scvf.outerNormalScalar()) == sign(velocitySelf)) ?
                                   elemVolVars[insideScvIdx] : elemVolVars[scvf.outsideScvIdx()] ;

           result += velocitySelf * velocitySelf * upVolVars.density() * sign(scvf.outerNormalScalar()) * scvf.area() ;
       }
       return result;
   }

   /*!
   * \brief Returns the tangential part of the momentum flux
   * \param scvf The sub control volume face
   * \param fvGeometry The finite-volume geometry
   * \param elemVolVars All volume variables for the element
   * \param globalFaceVars The face variables
   */
  FacePrimaryVariables computeTangetialMomentumFlux(const Problem& problem,
                                                    const SubControlVolumeFace& scvf,
                                                    const FVElementGeometry& fvGeometry,
                                                    const ElementVolumeVariables& elemVolVars,
                                                    const GlobalFaceVars& globalFaceVars)
  {
      FacePrimaryVariables tangentialFlux(0.0);

      // convenience function to get the velocity on a face
      auto velocity = [&globalFaceVars](const int dofIdx)
      {
          return globalFaceVars.faceVars(dofIdx).velocity();
      };

      // account for all sub-faces
      for(auto subFaceData : scvf.pairData())
      {
          const int eIdx = scvf.insideScvIdx();
          const auto& normalFace = fvGeometry.scvf(eIdx, subFaceData.localNormalFaceIdx);

          if(navierStokes)
              tangentialFlux += computeAdvectivePartOfTangentialMomentumFlux_(problem, scvf, normalFace, subFaceData, elemVolVars, velocity);

          tangentialFlux += computeDiffusivePartOfTangentialMomentumFlux_(problem, scvf, normalFace, subFaceData, elemVolVars, velocity);
      }
      return tangentialFlux;
  }

private:

  template<class SubFaceData, class VelocityHelper>
  FacePrimaryVariables computeAdvectivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const SubFaceData& subFaceData,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     VelocityHelper velocity)
  {
      const Scalar transportingVelocity = velocity(subFaceData.normalPair.first);
      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      const bool innerElementIsUpstream = ( sign(normalFace.outerNormalScalar()) == sign(transportingVelocity) );

      const auto& upVolVars = innerElementIsUpstream ? elemVolVars[insideScvIdx] : elemVolVars[outsideScvIdx];

      Scalar transportedVelocity(0.0);

      if(innerElementIsUpstream)
          transportedVelocity = velocity(scvf.dofIndex());
      else
      {
          const int outerDofIdx = subFaceData.outerParallelFaceDofIdx;
          if(outerDofIdx >= 0)
              transportedVelocity = velocity(outerDofIdx);
          else // this is the case when the outer parallal dof would lie outside the domain
              transportedVelocity = problem.dirichletAtPos(scvf.center())[faceIdx][scvf.directionIndex()];
      }

      const Scalar momentum = upVolVars.density() * transportedVelocity;
      const int sgn = sign(normalFace.outerNormalScalar());

      return transportingVelocity * momentum * sgn * normalFace.area() * 0.5;
  }

  template<class SubFaceData, class VelocityHelper>
  FacePrimaryVariables computeDiffusivePartOfTangentialMomentumFlux_(const Problem& problem,
                                                                     const SubControlVolumeFace& scvf,
                                                                     const SubControlVolumeFace& normalFace,
                                                                     const SubFaceData& subFaceData,
                                                                     const ElementVolumeVariables& elemVolVars,
                                                                     VelocityHelper velocity)
  {
      FacePrimaryVariables tangentialDiffusiveFlux(0.0);

      const auto normalDirIdx = normalFace.directionIndex();
      const auto insideScvIdx = normalFace.insideScvIdx();
      const auto outsideScvIdx = normalFace.outsideScvIdx();

      const auto& insideVolVars = elemVolVars[insideScvIdx];
      const auto& outsideVolVars = elemVolVars[outsideScvIdx];

      // the averaged viscosity at the face normal to our face of interest (where we assemble the face residual)
      const Scalar muAvg = (insideVolVars.viscosity() + outsideVolVars.viscosity()) * 0.5;

      // the normal derivative
      const int innerNormalVelocityIdx = subFaceData.normalPair.first;
      const int outerNormalVelocityIdx = subFaceData.normalPair.second;

      const Scalar innerNormalVelocity = velocity(innerNormalVelocityIdx);

      const Scalar outerNormalVelocity = outerNormalVelocityIdx >= 0 ?
                                  velocity(outerNormalVelocityIdx) :
                                  problem.dirichletAtPos(subFaceData.virtualOuterNormalFaceDofPos)[faceIdx][normalDirIdx];

      const Scalar normalDeltaV = scvf.normalInPosCoordDir() ?
                                    (outerNormalVelocity - innerNormalVelocity) :
                                    (innerNormalVelocity - outerNormalVelocity);

      const Scalar normalDerivative = normalDeltaV / subFaceData.normalDistance;
      tangentialDiffusiveFlux -= muAvg * normalDerivative;

      // the parallel derivative
      const Scalar innerParallelVelocity = velocity(scvf.dofIndex());

      const int outerParallelFaceDofIdx = subFaceData.outerParallelFaceDofIdx;
      const Scalar outerParallelVelocity = outerParallelFaceDofIdx >= 0 ?
                                           velocity(outerParallelFaceDofIdx) :
                                           problem.dirichletAtPos(subFaceData.virtualOuterParallelFaceDofPos)[faceIdx][scvf.directionIndex()];

      const Scalar parallelDeltaV = normalFace.normalInPosCoordDir() ?
                                   (outerParallelVelocity - innerParallelVelocity) :
                                   (innerParallelVelocity - outerParallelVelocity);

      const Scalar parallelDerivative = parallelDeltaV / subFaceData.parallelDistance;
      tangentialDiffusiveFlux -= muAvg * parallelDerivative;

      const Scalar sgn = sign(normalFace.outerNormalScalar());
      return tangentialDiffusiveFlux * sgn * normalFace.area() * 0.5;
  }
};


// specialization for miscible, isothermal flow
template<class TypeTag>
class FreeFlowFluxVariablesImpl<TypeTag, true, false> : public FreeFlowFluxVariablesImpl<TypeTag, false, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
    static constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    static constexpr bool useMoles = true;

    //! The index of the component balance equation that gets replaced with the total mass balance
    static const int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);

    using ParentType = FreeFlowFluxVariablesImpl<TypeTag, false, false>;

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        conti0EqIdx = Indices::conti0EqIdx
    };

public:
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const GlobalFaceVars& globalFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux(0.0);

        flux += advectiveFluxForCellCenter_(problem, fvGeometry, elemVolVars, globalFaceVars, scvf);
        flux += diffusiveFluxForCellCenter_(problem, fvGeometry, elemVolVars, scvf);
        return flux;
    }

private:

    CellCenterPrimaryVariables advectiveFluxForCellCenter_(const Problem& problem,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const GlobalFaceVars& globalFaceVars,
                                                          const SubControlVolumeFace &scvf)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();

        const bool insideIsUpstream = sign(scvf.outerNormalScalar()) == sign(velocity) ? true : false;
        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;

        const Scalar upWindWeight = GET_PROP_VALUE(TypeTag, ImplicitUpwindWeight);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            auto eqIdx = conti0EqIdx + compIdx;

            const Scalar upstreamDensity = useMoles ? upstreamVolVars.molarDensity() : upstreamVolVars.density();
            const Scalar downstreamDensity = useMoles ? downstreamVolVars.molarDensity() : downstreamVolVars.density();
            const Scalar upstreamFraction = useMoles ? upstreamVolVars.moleFraction(0, compIdx) : upstreamVolVars.massFraction(0, compIdx);
            const Scalar downstreamFraction = useMoles ? downstreamVolVars.moleFraction(0, compIdx) : downstreamVolVars.massFraction(0, compIdx);

            Scalar advFlux = 0.0;

            if(scvf.boundary() && eqIdx > conti0EqIdx)
            {
                const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isDirichlet(eqIdx))
                    advFlux = upstreamDensity * problem.dirichletAtPos(scvf.center())[eqIdx] * velocity;
                if(bcTypes.isOutflow(eqIdx))
                    advFlux = upstreamDensity * upstreamFraction * velocity;
            }
            else
                advFlux = (upWindWeight * upstreamDensity * upstreamFraction +
                          (1.0 - upWindWeight) * downstreamDensity * downstreamFraction) * velocity;

            if (eqIdx != replaceCompEqIdx)
                flux[eqIdx] += advFlux;

            // in case one balance is substituted by the total mass balance
            if (replaceCompEqIdx < numComponents)
                flux[replaceCompEqIdx] += advFlux;
        }

        flux *= scvf.area() * sign(scvf.outerNormalScalar());
        return flux;
    }


    CellCenterPrimaryVariables diffusiveFluxForCellCenter_(const Problem& problem,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace &scvf)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideDensity = useMoles ? insideVolVars.molarDensity() : insideVolVars.density();

        for(int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            auto eqIdx = conti0EqIdx + compIdx;

            const Scalar tij = transmissibility_(problem, fvGeometry, elemVolVars, scvf, compIdx);
            const Scalar insideFraction = useMoles ? insideVolVars.moleFraction(0, compIdx) : insideVolVars.massFraction(0, compIdx);

            if(scvf.boundary())
            {
                const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isOutflow(eqIdx))
                    return flux;
                else if(bcTypes.isNeumann(eqIdx))
                    return flux; // TODO: implement neumann
                else
                {
                    const Scalar dirichletFraction = problem.dirichletAtPos(scvf.center())[eqIdx];
                    flux[eqIdx] = insideDensity * tij * (insideFraction - dirichletFraction);
                }
            }
            else
            {
                const Scalar outsideDensity = useMoles ? outsideVolVars.molarDensity() : outsideVolVars.density();
                const Scalar avgDensity = 0.5*(insideDensity + outsideDensity);
                const Scalar outsideFraction = useMoles ? outsideVolVars.moleFraction(0, compIdx) : outsideVolVars.massFraction(0, compIdx);
                flux[eqIdx] = avgDensity * tij * (insideFraction - outsideFraction);
            }
        }

        if (replaceCompEqIdx >= numComponents)
        {
            const Scalar cumulativeFlux = std::accumulate(flux.begin(), flux.end(), 0.0);
            flux[0] = - cumulativeFlux;
        }

        return flux;
    }

    static Scalar transmissibility_(const Problem& problem,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int compIdx)
    {
        Scalar tij = 0.0;
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
        const Scalar insideD = insideVolVars.diffusionCoefficient(0, compIdx);
        const Scalar ti = calculateOmega_(insideDistance, insideD, 1.0);

        if(scvf.boundary())
            tij = scvf.area() * ti;
        else
        {
            const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const Scalar outsideD = outsideVolVars.diffusionCoefficient(0, compIdx);
            const Scalar tj = calculateOmega_(outsideDistance, outsideD, 1.0);

            tij = scvf.area()*(ti * tj)/(ti + tj);
        }
        return tij;
    }

    static Scalar calculateOmega_(const Scalar distance,
                                  const Scalar D,
                                  const Scalar extrusionFactor)
    {
        Scalar omega = D / distance;
        omega *= extrusionFactor;

        return omega;
    }
};

} // end namespace

#endif
