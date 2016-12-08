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
 * \brief This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCTpfa method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvFace,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvFace];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvFace.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvFace.outsideScvIdx()];

        auto hInside = insideVolVars.pressure(phaseIdx);
        auto hOutside = scvFace.numOutsideScvs() <= 1 ? outsideVolVars.pressure(phaseIdx)
                        : branchingFacetPressure_(phaseIdx, problem, element, fvGeometry,
                                                  elemVolVars, scvFace,
                                                  elemFluxVarsCache, hInside);

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density over all neighboring elements
            const auto rho = scvFace.numOutsideScvs() <= 1 ? (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5
                             : branchingFacetDensity_(phaseIdx, elemVolVars, scvFace, insideVolVars.density(phaseIdx));

            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);

            hInside -= rho*(gInside*xInside);

            // and the outside neighbor
            if (scvFace.boundary() || scvFace.numOutsideScvs() > 1)
            {
                const auto xOutside = scvFace.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
            else
            {
                const auto outsideScvIdx = scvFace.outsideScvIdx();
                // as we assemble fluxes from the neighbor to our element the outside index
                // refers to the scv of our element, so we use the scv method
                const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
                const auto xOutside = outsideScv.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
        }

        return fluxVarsCache.tij()*(hInside - hOutside);
    }

    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvFace)
    {
        if (scvFace.boundary())
            return Stencil({scvFace.insideScvIdx()});
        else if (scvFace.numOutsideScvs() > 1)
        {
            Stencil stencil({scvFace.insideScvIdx()});
            for (unsigned int i = 0; i < scvFace.numOutsideScvs(); ++i)
                stencil.push_back(scvFace.outsideScvIdx(i));
            return stencil;
        }
        else
            return Stencil({scvFace.insideScvIdx(), scvFace.outsideScvIdx()});
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static Scalar calculateTransmissibilities(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const SubControlVolumeFace& scvFace)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto insideK = problem.spatialParams().intrinsicPermeability(insideScv, insideVolVars);
        Scalar ti = calculateOmega_(problem, scvFace, insideK, element, insideScv);

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvFace.boundary() || scvFace.numOutsideScvs() > 1)
        {
            tij = scvFace.area()*ti;
        }
        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto outsideElement = fvGeometry.globalFvGeometry().element(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto outsideK = problem.spatialParams().intrinsicPermeability(outsideScv, outsideVolVars);
            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*calculateOmega_(problem, scvFace, outsideK, outsideElement, outsideScv);
            else
                tj = calculateOmega_(problem, fvGeometry.flipScvf(scvFace.index()), outsideK, outsideElement, outsideScv);

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvFace.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    //! compute the transmissibility ti, overload for tensor permeabilites
    static Scalar calculateOmega_(const Problem& problem,
                                  const SubControlVolumeFace& scvFace,
                                  const DimWorldMatrix &K,
                                  const Element& element,
                                  const SubControlVolume &scv)
    {
        GlobalPosition Knormal;
        K.mv(scvFace.unitOuterNormal(), Knormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }

    //! compute the transmissibility ti, overload for scalar permeabilites
    static Scalar calculateOmega_(const Problem& problem,
                                  const SubControlVolumeFace& scvFace,
                                  const Scalar K,
                                  const Element& element,
                                  const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }

    //! compute the pressure at branching facets for network grids
    static Scalar branchingFacetPressure_(int phaseIdx,
                                          const Problem& problem,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars,
                                          const SubControlVolumeFace& scvFace,
                                          const ElementFluxVarsCache& elemFluxVarsCache,
                                          Scalar insideP)
    {
        const auto& insideFluxVarsCache = elemFluxVarsCache[scvFace];
        Scalar sumTi(insideFluxVarsCache.tij());
        Scalar sumPTi(insideFluxVarsCache.tij()*insideP);

        for (unsigned int i = 0; i < scvFace.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx(i);
            const auto& flippedScvf = fvGeometry.flipScvf(scvFace.index(), i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
            sumTi += outsideFluxVarsCache.tij();
            sumPTi += outsideFluxVarsCache.tij()*outsideVolVars.pressure(phaseIdx);
        }
        return sumPTi/sumTi;
    }

    //! compute the density at branching facets for network grids as arithmetic mean
    static Scalar branchingFacetDensity_(int phaseIdx,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace& scvFace,
                                         Scalar insideRho)
    {
        Scalar rho(insideRho);
        for (unsigned int i = 0; i < scvFace.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            rho += outsideVolVars.density(phaseIdx);
        }
        return rho/(scvFace.numOutsideScvs()+1);
    }
};

} // end namespace

#endif
