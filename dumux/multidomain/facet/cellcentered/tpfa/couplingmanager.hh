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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \copydoc Dumux::CCTpfaFacetCouplingManager
 */

#ifndef DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH
#define DUMUX_CCTPFA_FACETCOUPLING_MANAGER_HH

#include <algorithm>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

//! Forward declaration of the manager coupling three domains
template<class MDTraits, class CouplingMapper>
class CCTpfaFacetCouplingThreeDomainManager;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation
 *        is to be used in conjunction with models using the cell-centered tpfa scheme.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam idOffset The offset to be used for the domain ids. This is used to specify
 *                  which of the present grids on the hierarchy are to be managed by this
 *                  class. For instance, if a bulk, a facet and an edge grid is present and
 *                  you want to use this coupling manager for the coupling between the facet
 *                  and the edge grid, you should provide an offet of 1.
 */
template<class MDTraits, class CouplingMapper, std::size_t idOffset = 0>
class CCTpfaFacetCouplingManager : public CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // convenience aliases and instances of the two domain ids
    using BulkIdType = typename MDTraits::template DomainIdx<idOffset+0>;
    using LowDimIdType = typename MDTraits::template DomainIdx<idOffset+1>;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto lowDimId = LowDimIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, NumEqVector);
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using VolumeVariables = typename ElementVolumeVariables<id>::VolumeVariables;
    template<std::size_t id> using GridFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridFluxVariablesCache<id>::LocalView;

    /*!
     * \brief The coupling context of the bulk domain. Contains all data of the lower-
     *        dimensional domain which is required for the computation of a bulk element
     *        residual. This boils down to the geometries and volume variables of all
     *        lower-dimensional elements connected to a given bulk element.
     */
    struct BulkCouplingContext
    {
        bool isSet;
        IndexType< bulkId > elementIdx;
        std::vector< FVElementGeometry<lowDimId> > lowDimFvGeometries;
        std::vector< VolumeVariables<lowDimId> > lowDimVolVars;

        void reset()
        {
            lowDimFvGeometries.clear();
            lowDimVolVars.clear();
            isSet = false;
        }
    };

    /*!
     * \brief The coupling context of the lower-dimensional (codim 1) domain. Contains
     *        all data of the bulk domain which is required for computation of element
     *        residuals in the lower-dimensional domain. This is essentially everything
     *        that is necessary to compute the fluxes in bulk domain entering a given
     *        lower-dimensional element. Thus, we store and bind the local views of one
     *        of the neighboring elements, which will be enough to compute the fluxes
     *        stemming from all neighboring elements.
     *
     * \note We need unique ptrs here because the local views have no default constructor.
     */
    struct LowDimCouplingContext
    {
        bool isSet;
        IndexType< lowDimId > elementIdx;
        std::unique_ptr< FVElementGeometry<bulkId> > bulkFvGeometry;
        std::unique_ptr< ElementVolumeVariables<bulkId> > bulkElemVolVars;
        std::unique_ptr< ElementFluxVariablesCache<bulkId> > bulkElemFluxVarsCache;
        std::unique_ptr< LocalResidual<bulkId> > bulkLocalResidual;

        void reset()
        {
            bulkFvGeometry.reset(nullptr);
            bulkElemVolVars.reset(nullptr);
            bulkElemFluxVarsCache.reset(nullptr);
            isSet = false;
        }
    };

public:

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j = (i == bulkId) ? lowDimId : bulkId>
    using CouplingStencilType = typename CouplingMapper::template Stencil<j>;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the bulk domain
     * \param lowDimProblem The problem to be solved on the lower-dimensional domain
     * \param couplingMapper The mapper object containing the connectivity between the domains
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<lowDimId> > lowDimProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        couplingMapperPtr_ = couplingMapper;

        // set up tuple containing the sub-problems
        problemTuple_ = std::make_tuple(bulkProblem, lowDimProblem);

        // copy the solution vector
        ParentType::updateSolution(curSol);

        // determine all bulk elements/scvfs that couple to low dim elements
        bulkElemIsCoupled_.resize(bulkProblem->fvGridGeometry().gridView().size(0), false);
        bulkScvfIsCoupled_.resize(bulkProblem->fvGridGeometry().numScvf(), false);

        const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        for (const auto& entry : bulkMap)
        {
            bulkElemIsCoupled_[entry.first] = true;
            for (const auto& scvfs : entry.second.couplingScvfs)
                bulkScvfIsCoupled_[scvfs[0]] = true;
        }
    }

    /*!
     * \brief The coupling stencil of a given bulk domain element.
     */
    const CouplingStencilType<bulkId>& couplingStencil(BulkIdType domainI,
                                                       const Element<bulkId>& element,
                                                       LowDimIdType domainJ) const
    {
        const auto eIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(element);

        if (bulkElemIsCoupled_[eIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            auto it = map.find(eIdx);
            assert(it != map.end());
            if (it == map.end()) DUNE_THROW(Dune::InvalidStateException, "not found!");
            return it->second.couplingStencil;
        }

        return getEmptyStencil(lowDimId);
    }

    /*!
     * \brief The coupling stencil of the lower-dimensional domain with the bulk domain.
     */
    const CouplingStencilType<lowDimId>& couplingStencil(LowDimIdType domainI,
                                                         const Element<lowDimId>& element,
                                                         BulkIdType domainJ) const
    {
        const auto eIdx = problem<lowDimId>().fvGridGeometry().elementMapper().index(element);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);
        auto it = map.find(eIdx);
        if (it != map.end()) return it->second.couplingStencil;
        else return getEmptyStencil(bulkId);
    }

    /*!
     * \brief returns true if a bulk scvf coincides with a facet element.
     */
    bool isCoupled(const Element<bulkId>& element,
                   const SubControlVolumeFace<bulkId>& scvf) const
    { return bulkScvfIsCoupled_[scvf.index()]; }

    /*!
     * \brief returns the vol vars of a lower-dimensional element coinciding with a bulk scvf.
     */
    const VolumeVariables<lowDimId>& getLowDimVolVars(const Element<bulkId>& element,
                                                      const SubControlVolumeFace<bulkId>& scvf) const
    {
        assert(bulkContext_.isSet);
        assert(bulkScvfIsCoupled_[scvf.index()]);
        assert(scvf.insideScvIdx() == problem<bulkId>().fvGridGeometry().elementMapper().index(element));

        const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        const auto& couplingData = map.find(scvf.insideScvIdx())->second;
        const auto lowDimElemIdx = couplingData.getCoupledFacetElementIdx(scvf.index());

        const auto& s = map.find(bulkContext_.elementIdx)->second.couplingStencil;
        const auto& idxInContext = std::distance( s.begin(), std::find(s.begin(), s.end(), lowDimElemIdx) );
        assert(std::find(s.begin(), s.end(), lowDimElemIdx) != s.end()); if (std::find(s.begin(), s.end(), lowDimElemIdx) == s.end()) DUNE_THROW(Dune::InvalidStateException, "hure element!");
        return bulkContext_.lowDimVolVars[idxInContext];
    }

    /*!
     * \brief Evaluates the coupling element residual of a bulk domain element with respect
     *        to a dof in the lower-dimensional domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the bulk element facets that coincide with the lower-dimensional
     *        element whose dof idx is dofIdxGlobalJ.
     */
    template< class BulkLocalAssembler >
    typename LocalResidual<bulkId>::ElementResidualVector
    evalCouplingResidual(BulkIdType,
                         const BulkLocalAssembler& bulkLocalAssembler,
                         LowDimIdType,
                         IndexType<lowDimId> dofIdxGlobalJ)
    {
        typename LocalResidual<bulkId>::ElementResidualVector res(1);
        res = 0.0;
        res[0] = evalFluxToFacetElement_(bulkLocalAssembler.element(),
                                         bulkLocalAssembler.fvGeometry(),
                                         bulkLocalAssembler.curElemVolVars(),
                                         bulkLocalAssembler.elemFluxVarsCache(),
                                         bulkLocalAssembler.localResidual(),
                                         dofIdxGlobalJ);
        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of a lower-dimensional domain element
     *        with respect to a dof in the bulk domain (dofIdxGlobalJ). This is essentially
     *        the fluxes across the facets of the neighboring bulk element that coincide with
     *        the given element.
     */
    template< class LowDimLocalAssembler >
    typename LocalResidual<lowDimId>::ElementResidualVector
    evalCouplingResidual(LowDimIdType,
                         const LowDimLocalAssembler& lowDimLocalAssembler,
                         BulkIdType,
                         IndexType<bulkId> dofIdxGlobalJ)
    {
        // make sure this is called for the element for which the context was set
        assert(lowDimContext_.isSet);
        assert(problem<lowDimId>().fvGridGeometry().elementMapper().index(lowDimLocalAssembler.element()) == lowDimContext_.elementIdx);

        // since we use cc schemes: dof index = element index
        const auto elementJ = problem<bulkId>().fvGridGeometry().element(dofIdxGlobalJ);
        typename LocalResidual<lowDimId>::ElementResidualVector res(1);
        res = 0.0;
        res[0] = evalFluxToFacetElement_(elementJ,
                                         *lowDimContext_.bulkFvGeometry,
                                         *lowDimContext_.bulkElemVolVars,
                                         *lowDimContext_.bulkElemFluxVarsCache,
                                         *lowDimContext_.bulkLocalResidual,
                                         lowDimContext_.elementIdx);
        res[0] *= -1.0;
        return res;
    }

    /*!
     * \brief Computes the sources in a lower-dimensional element stemming from the bulk domain.
     */
    NumEqVector<lowDimId> evalSourcesFromBulk(const Element<lowDimId>& element,
                                              const FVElementGeometry<lowDimId>& fvGeometry,
                                              const ElementVolumeVariables<lowDimId>& elemVolVars,
                                              const SubControlVolume<lowDimId>& scv)
    {
        NumEqVector<lowDimId> sources(0.0);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);
        auto it = map.find(lowDimContext_.elementIdx);
        if (it == map.end())
            return sources;

        // make sure this is called for the element for which the context was set
        assert(lowDimContext_.isSet);
        assert(problem<lowDimId>().fvGridGeometry().elementMapper().index(element) == lowDimContext_.elementIdx);

        for (const auto& embedment : it->second.embedments)
            sources += evalFluxToFacetElement_(problem<bulkId>().fvGridGeometry().element(embedment.first),
                                               *lowDimContext_.bulkFvGeometry,
                                               *lowDimContext_.bulkElemVolVars,
                                               *lowDimContext_.bulkElemFluxVarsCache,
                                               *lowDimContext_.bulkLocalResidual,
                                               lowDimContext_.elementIdx);

        return sources;
    }

    /*!
     * \brief For the assembly of the element residual of a bulk domain element
     *        we need to prepare all variables of lower-dimensional domain elements
     *        that are coupled to the given bulk element
     */
    template< class Assembler >
    void bindCouplingContext(BulkIdType, const Element<bulkId>& element, const Assembler& assembler)
    {
        // clear context
        bulkContext_.reset();

        // set index in context in any case
        const auto bulkElemIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(element);
        bulkContext_.elementIdx = bulkElemIdx;

        // if element is coupled, actually set the context
        if (bulkElemIsCoupled_[bulkElemIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);

            auto it = map.find(bulkElemIdx); assert(it != map.end()); if (it == map.end()) DUNE_THROW(Dune::InvalidStateException, "as9d");
            const auto stencilSize = it->second.couplingStencil.size();
            bulkContext_.lowDimFvGeometries.reserve(stencilSize);
            bulkContext_.lowDimVolVars.reserve(stencilSize);

            for (const auto lowDimIdx : it->second.couplingStencil)
            {
                const auto& ldGridGeometry = problem<lowDimId>().fvGridGeometry();

                const auto elemJ = ldGridGeometry.element(lowDimIdx);
                auto fvGeom = localView(ldGridGeometry);
                fvGeom.bindElement(elemJ);

                const auto elemSol = elementSolution(elemJ, this->curSol()[lowDimId], ldGridGeometry);
                VolumeVariables<lowDimId> volVars;
                volVars.update(elemSol, problem<lowDimId>(), elemJ, fvGeom.scv(lowDimIdx));

                bulkContext_.isSet = true;
                bulkContext_.lowDimFvGeometries.emplace_back( std::move(fvGeom) );
                bulkContext_.lowDimVolVars.emplace_back( std::move(volVars) );
            }
        }
    }

    /*!
     * \brief For the assembly of the element residual of a bulk domain element
     *        we need to prepare the local views of one of the neighboring bulk
     *        domain elements. These are used later to compute the fluxes across
     *        the faces over which the coupling occurs
     * \note Since the low-dim coupling residua are fluxes stemming from
     *       the bulk domain, we have to prepare the bulk coupling context
     *       for the neighboring element (where fluxes are calculated) as well.
     */
    template< class Assembler >
    void bindCouplingContext(LowDimIdType, const Element<lowDimId>& element, const Assembler& assembler)
    {
        // reset contexts
        bulkContext_.reset();
        lowDimContext_.reset();

        // set index in context in any case
        const auto lowDimElemIdx = problem<lowDimId>().fvGridGeometry().elementMapper().index(element);
        lowDimContext_.elementIdx = lowDimElemIdx;

        const auto& map = couplingMapperPtr_->couplingMap(lowDimId, bulkId);
        auto it = map.find(lowDimElemIdx);

        // if element is coupled, actually set the context
        if (it != map.end())
        {
            // first bind the low dim context for the first neighboring bulk element
            const auto& bulkGridGeom = problem<bulkId>().fvGridGeometry();
            const auto bulkElem = bulkGridGeom.element(it->second.couplingStencil[0]);
            bindCouplingContext(bulkId, bulkElem, assembler);

            // then simply bind the local views of that first neighbor
            auto bulkFvGeom = localView(bulkGridGeom);
            auto bulkElemVolVars = localView(assembler.gridVariables(bulkId).curGridVolVars());
            auto bulkElemFluxVarsCache = localView(assembler.gridVariables(bulkId).gridFluxVarsCache());

            bulkFvGeom.bind(bulkElem);
            bulkElemVolVars.bind(bulkElem, bulkFvGeom, this->curSol()[bulkId]);
            bulkElemFluxVarsCache.bind(bulkElem, bulkFvGeom, bulkElemVolVars);

            lowDimContext_.isSet = true;
            lowDimContext_.bulkFvGeometry = std::make_unique< FVElementGeometry<bulkId> >(bulkFvGeom);
            lowDimContext_.bulkElemVolVars = std::make_unique< ElementVolumeVariables<bulkId> >(bulkElemVolVars);
            lowDimContext_.bulkElemFluxVarsCache = std::make_unique< ElementFluxVariablesCache<bulkId> >(bulkElemFluxVarsCache);
            lowDimContext_.bulkLocalResidual = std::make_unique< LocalResidual<bulkId> >(assembler.localResidual(bulkId));
        }
    }

    /*!
     * \brief After deflecting the solution of the lower-dimensional domain,
     *        we have to update the element volume variables object if the context.
     */
    template< class BulkLocalAssembler >
    void updateCouplingContext(BulkIdType domainI,
                               const BulkLocalAssembler& bulkLocalAssembler,
                               LowDimIdType domainJ,
                               IndexType<lowDimId> dofIdxGlobalJ,
                               const PrimaryVariables<lowDimId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, bulkLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // skip the rest if context is empty
        if (bulkContext_.isSet)
        {
            // since we use cc schemes: dof index = element index
            const auto& lowDimGridGeom = problem<lowDimId>().fvGridGeometry();
            const auto elementJ = lowDimGridGeom.element(dofIdxGlobalJ);

            // update vol vars in context
            const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            const auto& couplingStencil = map.find(bulkContext_.elementIdx)->second.couplingStencil;
            auto it = std::find(couplingStencil.begin(), couplingStencil.end(), dofIdxGlobalJ);

            assert(it != couplingStencil.end()); if (it == couplingStencil.end()) DUNE_THROW(Dune::InvalidStateException, "asoid");
            const auto idxInContext = std::distance(couplingStencil.begin(), it);
            const auto& lowDimScv = bulkContext_.lowDimFvGeometries[idxInContext].scv(dofIdxGlobalJ);
            const auto elemSol = elementSolution(elementJ, this->curSol()[lowDimId], lowDimGridGeom);
            bulkContext_.lowDimVolVars[idxInContext].update(elemSol, problem<lowDimId>(), elementJ, lowDimScv);
        }
    }

    /*!
     * \brief Update the coupling context for a derivative bulk -> bulk.
     *        Here, we simply have to update the solution.
     */
    template< class BulkLocalAssembler >
    void updateCouplingContext(BulkIdType domainI,
                               const BulkLocalAssembler& bulkLocalAssembler,
                               BulkIdType domainJ,
                               IndexType<bulkId> dofIdxGlobalJ,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, bulkLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief After deflecting the solution of the bulk domain, we have to update
     *        the element volume variables and transmissibilities of the neighboring
     *        bulk element stored in the context.
     */
    template< class LowDimLocalAssembler >
    void updateCouplingContext(LowDimIdType domainI,
                               const LowDimLocalAssembler& lowDimLocalAssembler,
                               BulkIdType domainJ,
                               IndexType<bulkId> dofIdxGlobalJ,
                               const PrimaryVariables<bulkId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, lowDimLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // skip the rest if context is empty
        if (lowDimContext_.isSet)
        {
            // since we use cc schemes: dof index = element index
            const auto& bulkGridGeom = problem<bulkId>().fvGridGeometry();
            const auto elementJ = bulkGridGeom.element(dofIdxGlobalJ);

            // update corresponding vol vars in context
            const auto& scv = lowDimContext_.bulkFvGeometry->scv(dofIdxGlobalJ);
            const auto elemSol = elementSolution(elementJ, this->curSol()[bulkId], bulkGridGeom);
            (*lowDimContext_.bulkElemVolVars)[dofIdxGlobalJ].update(elemSol, problem<bulkId>(), elementJ, scv);

            // update the element flux variables cache (tij depend on low dim values)
            if (dofIdxGlobalJ == bulkContext_.elementIdx)
                lowDimContext_.bulkElemFluxVarsCache->update( elementJ, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
            else
                lowDimContext_.bulkElemFluxVarsCache->update( problem<bulkId>().fvGridGeometry().element(bulkContext_.elementIdx),
                                                              *lowDimContext_.bulkFvGeometry,
                                                              *lowDimContext_.bulkElemVolVars );
        }
    }

    /*!
     * \brief After deflecting the solution of the lower-dimensional domain has been deflected
     *        during the assembly of the element residual of a lower-dimensional element, we
     *        have to communicate this to the volume variables stored in the context as well
     *        as the transmissibilities.
     */
    template< class LowDimLocalAssembler >
    void updateCouplingContext(LowDimIdType domainI,
                               const LowDimLocalAssembler& lowDimLocalAssembler,
                               LowDimIdType domainJ,
                               IndexType<lowDimId> dofIdxGlobalJ,
                               const PrimaryVariables<lowDimId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate deflected solution
        ParentType::updateCouplingContext(domainI, lowDimLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // skip the rest if context is empty
        if (lowDimContext_.isSet)
        {
            const auto& lowDimGridGeom = problem<lowDimId>().fvGridGeometry();

            assert(bulkContext_.isSet);
            assert(lowDimContext_.elementIdx == lowDimGridGeom.elementMapper().index(lowDimLocalAssembler.element()));

            // update the corresponding vol vars in the bulk context
            const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
            const auto& couplingStencil = bulkMap.find(bulkContext_.elementIdx)->second.couplingStencil;
            auto it = std::find(couplingStencil.begin(), couplingStencil.end(), lowDimContext_.elementIdx);

            assert(it != couplingStencil.end()); if (it == couplingStencil.end()) DUNE_THROW(Dune::InvalidStateException, "q90u");
            const auto idxInContext = std::distance(couplingStencil.begin(), it);
            const auto& lowDimScv = bulkContext_.lowDimFvGeometries[idxInContext].scv(lowDimContext_.elementIdx);
            const auto lowDimElement = lowDimGridGeom.element(lowDimContext_.elementIdx);
            const auto elemSol = elementSolution(lowDimElement, this->curSol()[lowDimId], lowDimGridGeom);
            bulkContext_.lowDimVolVars[idxInContext].update(elemSol, problem<lowDimId>(), lowDimElement, lowDimScv);

            // update the element flux variables cache (tij depend on low dim values)
            const auto contextElem = problem<bulkId>().fvGridGeometry().element(bulkContext_.elementIdx);
            lowDimContext_.bulkElemFluxVarsCache->update(contextElem, *lowDimContext_.bulkFvGeometry, *lowDimContext_.bulkElemVolVars);
        }
    }

    /*!
     * \brief Update the transmissibilities in the bulk domain after the coupling context changed
     * \note Specialization of the function for deactivated grid-wide volume variables caching
     */
    template< class BulkLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(BulkIdType domainI,
                                const BulkLocalAssembler& bulkLocalAssembler,
                                ElementVolumeVariables<bulkId>& elemVolVars,
                                UpdatableFluxVarCache& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed
        fluxVarsCache.update(bulkLocalAssembler.element(),
                             bulkLocalAssembler.fvGeometry(),
                             bulkLocalAssembler.curElemVolVars());
    }

    /*!
     * \brief Update the transmissibilities in the bulk domain after the coupling context changed
     * \note Specialization of the function for enabled grid-wide volume variables caching
     */
    template< class BulkLocalAssembler, class UpdatableFluxVarCache >
    void updateCoupledVariables(BulkIdType domainI,
                                const BulkLocalAssembler& bulkLocalAssembler,
                                GridVolumeVariables<bulkId>& gridVolVars,
                                UpdatableFluxVarCache& fluxVarsCache)
    {
        // update transmissibilities after low dim context has changed
        auto elemVolVars = localView(gridVolVars);
        elemVolVars.bind(bulkLocalAssembler.element(), bulkLocalAssembler.fvGeometry(), this->curSol()[bulkId]);
        fluxVarsCache.update(bulkLocalAssembler.element(), bulkLocalAssembler.fvGeometry(), elemVolVars);
    }

    /*!
     * \brief In the low dim domain the local views do not depend on coupling data.
     *        However, the overload has to be provided to avoid ambiguity.
     */
    template< class LowDimLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache >
    void updateCoupledVariables(LowDimIdType domainI,
                                const LowDimLocalAssembler& lowDimLocalAssembler,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& fluxVarsCache)
    { /*do nothing here*/ }

    //! Return a const reference to one of the problems
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const Problem<id>& problem() const { return *std::get<id-idOffset>(problemTuple_); }

    //! Return a reference to one of the problems
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    Problem<id>& problem() { return *std::get<id-idOffset>(problemTuple_); }

    //! Empty stencil to be returned for elements that aren't coupled
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getEmptyStencil(Dune::index_constant<id>) const
    { return std::get<id-idOffset>(emptyStencilTuple_); }

private:
    //! evaluates the bulk-facet exchange fluxes for a given facet element
    NumEqVector<bulkId> evalFluxToFacetElement_(const Element<bulkId>& elementI,
                                                const FVElementGeometry<bulkId>& fvGeometry,
                                                const ElementVolumeVariables<bulkId>& elemVolVars,
                                                const ElementFluxVariablesCache<bulkId>& elemFluxVarsCache,
                                                const LocalResidual<bulkId>& localResidual,
                                                IndexType<lowDimId> globalJ) const
    {
        const auto bulkElemIdx = problem<bulkId>().fvGridGeometry().elementMapper().index(elementI);

        assert(bulkContext_.isSet); if (!bulkContext_.isSet) DUNE_THROW(Dune::InvalidStateException, "aslkidj");
        assert(bulkElemIsCoupled_[bulkElemIdx]); if (!bulkElemIsCoupled_[bulkElemIdx]) DUNE_THROW(Dune::InvalidStateException, "aslkidj");

        NumEqVector<bulkId> coupledFluxes(0.0);
        const auto& map = couplingMapperPtr_->couplingMap(bulkId, lowDimId);
        const auto& couplingScvfs = map.find(bulkElemIdx)->second.getCoupledScvfs(globalJ);

        for (const auto& scvfIdx : couplingScvfs)
            coupledFluxes += localResidual.evalFlux(problem<bulkId>(),
                                                    elementI,
                                                    fvGeometry,
                                                    elemVolVars,
                                                    elemFluxVarsCache,
                                                    fvGeometry.scvf(scvfIdx));

        return coupledFluxes;
    }

    using BulkProblemPtr = std::shared_ptr< Problem<bulkId> >;
    using LowDimProblemPtr = std::shared_ptr< Problem<lowDimId> >;
    std::tuple<BulkProblemPtr, LowDimProblemPtr> problemTuple_;
    std::shared_ptr<CouplingMapper> couplingMapperPtr_;

    //! store bools for all bulk elements/scvfs that indicate if they
    //! are coupled, so that we don't have to search in the map every time
    std::vector<bool> bulkElemIsCoupled_;
    std::vector<bool> bulkScvfIsCoupled_;

    //! empty stencil to return for non-coupled elements
    using BulkStencil = typename CouplingMapper::template Stencil<bulkId>;
    using LowDimStencil = typename CouplingMapper::template Stencil<lowDimId>;
    std::tuple<BulkStencil, LowDimStencil> emptyStencilTuple_;

    //! The coupling contexts of the two domains
    BulkCouplingContext bulkContext_;
    LowDimCouplingContext lowDimContext_;
};

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation is
 *        to be used in conjunction with models using the cell-centered tpfa scheme and in problems
 *        where grids and coupling along the hierarchy from 3d cells to 1d edges is considered.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 */
template<class MDTraits, class CouplingMapper>
class CCTpfaFacetCouplingThreeDomainManager
      : public CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 0>,
        public CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 1>
{
    using BulkFacetManager = CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 0>;
    using FacetEdgeManager = CCTpfaFacetCouplingManager<MDTraits, CouplingMapper, 1>;

    // convenience aliases and instances of the domain ids
    using BulkIdType = typename MDTraits::template DomainIdx<0>;
    using FacetIdType = typename MDTraits::template DomainIdx<1>;
    using EdgeIdType = typename MDTraits::template DomainIdx<2>;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto facetId = FacetIdType();
    static constexpr auto edgeId = EdgeIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache::LocalView;

public:

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = typename std::conditional< (j > 1),
                                                            typename FacetEdgeManager::template CouplingStencilType<i, j>,
                                                            typename BulkFacetManager::template CouplingStencilType<i, j> >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the (3d) bulk domain
     * \param facetProblem The problem to be solved on the (2d) facet domain
     * \param edgeProblem The problem to be solved on the (1d) edge domain
     * \param couplingMapper The mapper object containing the connectivity between the domains
     * \tparam curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<facetId> > facetProblem,
              std::shared_ptr< Problem<edgeId> > edgeProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        BulkFacetManager::init(bulkProblem, facetProblem, couplingMapper, curSol);
        FacetEdgeManager::init(facetProblem, edgeProblem, couplingMapper, curSol);
    }

    //! Pull up functionalities from the parent classes
    using BulkFacetManager::couplingStencil;
    using FacetEdgeManager::couplingStencil;

    using BulkFacetManager::isCoupled;
    using FacetEdgeManager::isCoupled;

    using BulkFacetManager::getLowDimVolVars;
    using FacetEdgeManager::getLowDimVolVars;

    using BulkFacetManager::evalSourcesFromBulk;
    using FacetEdgeManager::evalSourcesFromBulk;

    using BulkFacetManager::evalCouplingResidual;
    using FacetEdgeManager::evalCouplingResidual;

    using BulkFacetManager::bindCouplingContext;
    using FacetEdgeManager::bindCouplingContext;

    using BulkFacetManager::updateCouplingContext;
    using FacetEdgeManager::updateCouplingContext;

    using BulkFacetManager::updateCoupledVariables;
    using FacetEdgeManager::updateCoupledVariables;

    /*!
     * \brief The coupling stencil of the bulk with the edge domain (empty stencil).
     */
    const CouplingStencilType<bulkId, edgeId>& couplingStencil(BulkIdType domainI,
                                                               const Element<bulkId>& element,
                                                               EdgeIdType domainJ) const
    { return FacetEdgeManager::getEmptyStencil(edgeId); }

    /*!
     * \brief The coupling stencil of the edge with the bulk domain (empty stencil).
     */
    const CouplingStencilType<edgeId, bulkId>& couplingStencil(EdgeIdType domainI,
                                                               const Element<edgeId>& element,
                                                               BulkIdType domainJ) const
    { return BulkFacetManager::getEmptyStencil(bulkId); }

    /*!
     * \brief updates the current solution. We have to overload this here
     *        to avoid ambiguity and update the solution in both managers
     */
    void updateSolution(const SolutionVector& sol)
    {
        BulkFacetManager::updateSolution(sol);
        FacetEdgeManager::updateSolution(sol);
    }

    /*!
     * \brief Interface for evaluating the coupling residual between the bulk and the edge domain.
     *        This is always zero as coupling only occurs between grids of codimension one. These
     *        overloads are provided by the two parent classes. However, we need this overload in
     *        order for the overload resolution not to fail.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || ((i==edgeId && j==bulkId))), int> = 0>
    typename LocalResidual<i>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<i> domainI,
                         const LocalAssembler& localAssembler,
                         Dune::index_constant<j> domainJ,
                         IndexType<j> dofIdxGlobalJ)
    {
        typename LocalResidual<i>::ElementResidualVector res(1);
        res = 0.0;
        return res;
    }

    /*!
     * \brief Interface for binding the coupling context for the facet domain. In this case
     *        we have to bind both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class Assembler >
    void bindCouplingContext(FacetIdType, const Element<facetId>& element, const Assembler& assembler)
    {
        BulkFacetManager::bindCouplingContext(facetId, element, assembler);
        FacetEdgeManager::bindCouplingContext(facetId, element, assembler);
    }

    /*!
     * \brief Interface for updating the coupling context of the facet domain. In this case
     *        we have to update both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class FacetLocalAssembler >
    void updateCouplingContext(FacetIdType domainI,
                               const FacetLocalAssembler& facetLocalAssembler,
                               FacetIdType domainJ,
                               IndexType<facetId> dofIdxGlobalJ,
                               const PrimaryVariables<facetId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        BulkFacetManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        FacetEdgeManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief Interface for updating the coupling context between the bulk and the edge domain.
     *        We do nothing here because coupling only occurs between grids of codimension one.
     *        We have to provide this overload as the overload resolution fails otherwise.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || (i==edgeId && j==bulkId)), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssembler& localAssembler,
                               Dune::index_constant<j> domainJ,
                               IndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    { /*do nothing here*/ }

    /*!
     * \brief Interface for updating the local views of the facet domain after updateCouplingContext
     *        the coupling context. In this case we have to forward the both managers as the facet
     *        domain is a part in both.
     */
    template< class FacetLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetIdType domainI,
                                const FacetLocalAssembler& facetLocalAssembler,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        BulkFacetManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
        FacetEdgeManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
    }

    //! Return a const reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    const Problem<id>& problem() const { return BulkFacetManager::template problem<id>(); }

    //! Return a reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    Problem<id>& problem() { return BulkFacetManager::template problem<id>(); }

    //! Return a const reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    const Problem<id>& problem() const { return FacetEdgeManager::template problem<id>(); }

    //! Return a reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    Problem<id>& problem() { return FacetEdgeManager::template problem<id>(); }

    /*!
     * \brief We have to provide the overload of this function in order to avoid
     *        ambiguity due to the inheritance of two classes. Models using facet
     *        coupling together with the cell-centered tpfa scheme have no extended
     *        jacobian pattern.
     */
    template<std::size_t id, class JacobianPattern>
    void extendJacobianPattern(Dune::index_constant<id> domainI, JacobianPattern& pattern) const
    {}

    /*!
     * \brief We have to provide the overload of this function in order to avoid
     *        ambiguity due to the inheritance of two classes. Models using facet
     *        coupling together with the cell-centered tpfa scheme have no extended
     *        jacobian pattern and thus no additional derivatives.
     */
    template<std::size_t i, class LocalAssemblerI, class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDomainDerivatives(Dune::index_constant<i> domainI,
                                         const LocalAssemblerI& localAssemblerI,
                                         const typename LocalAssemblerI::LocalResidual::ElementResidualVector& origResiduals,
                                         JacobianMatrixDiagBlock& A,
                                         GridVariables& gridVariables)
    {}
};

} // end namespace Dumux

#endif