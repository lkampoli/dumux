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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingManager
 */
#ifndef DUMUX_BOX_FACETCOUPLING_MANAGER_HH
#define DUMUX_BOX_FACETCOUPLING_MANAGER_HH

#include <algorithm>
#include <cassert>

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        where the coupling occurs across the facets of the bulk grid. This implementation
 *        is to be used in conjunction with models using the box scheme in the bulk domain.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 */
template<class MDTraits, class CouplingMapper, std::size_t bulkDomainId, std::size_t lowDimDomainId>
class FacetCouplingManager<MDTraits, CouplingMapper, bulkDomainId, lowDimDomainId, DiscretizationMethod::box>
: public virtual CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // convenience aliases and instances of the two domain ids
    using BulkIdType = typename MDTraits::template DomainIdx<bulkDomainId>;
    using LowDimIdType = typename MDTraits::template DomainIdx<lowDimDomainId>;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto lowDimId = LowDimIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, NumEqVector);
    template<std::size_t id> using ElementBoundaryTypes = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementBoundaryTypes);
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

    // extract corresponding grid ids from the mapper
    static constexpr int bulkDim = GridView<bulkDomainId>::dimension;
    static constexpr int lowDimDim = GridView<lowDimDomainId>::dimension;
    static constexpr auto bulkGridId = CouplingMapper::template gridId<bulkDim>();
    static constexpr auto lowDimGridId = CouplingMapper::template gridId<lowDimDim>();

    /*!
     * \brief The coupling context of the bulk domain. Contains all data of the lower-
     *        dimensional domain which is required for the computation of a bulk element
     *        residual. This boils down to the geometries of all lower-dimensional elements
     *        connected to a given bulk element, as well as their element volume variables.
     */
    struct BulkCouplingContext
    {
        bool isSet;
        IndexType< bulkId > elementIdx;
        std::vector< FVElementGeometry<lowDimId> > lowDimFvGeometries;
        std::vector< ElementVolumeVariables<lowDimId> > lowDimElemVolVars;

        void reset()
        {
            lowDimFvGeometries.clear();
            lowDimElemVolVars.clear();
            isSet = false;
        }
    };

    /*!
     * \brief The coupling context of the lower-dimensional (codim 1) domain. Contains
     *        all data of the bulk domain which is required for computation of element
     *        residuals in the lower-dimensional domain. This is essentially everything
     *        that is necessary to compute the fluxes in bulk domain entering a given
     *        lower-dimensional element. Thus, we store and bind the local views of the
     *        neighboring elements.
     *
     * \note We need unique ptrs here because the local views have no default constructor.
     */
    struct LowDimCouplingContext
    {
    private:
        // For dim == dimworld in bulk domain, we know there is always going
        // to be two bulk neighbors for a lower-dimensional element. Choose
        // the container type depending on this
        static constexpr int bulkDimWorld = GridView<bulkId>::dimensionworld;
        static constexpr bool bulkIsSurfaceGrid = bulkDim != bulkDimWorld;

        template< class T >
        using ContainerType = typename std::conditional< bulkIsSurfaceGrid,
                                                         std::vector< T >,
                                                         std::array< T, 2 > >::type;

    public:
        bool isSet;
        IndexType< lowDimId > elementIdx;
        ContainerType< std::unique_ptr<FVElementGeometry<bulkId>> > bulkFvGeometries;
        ContainerType< std::unique_ptr<ElementVolumeVariables<bulkId>> > bulkElemVolVars;
        ContainerType< std::unique_ptr<ElementFluxVariablesCache<bulkId>> > bulkElemFluxVarsCache;
        std::unique_ptr< LocalResidual<bulkId> > bulkLocalResidual;
        ContainerType< ElementBoundaryTypes<bulkId> > bulkElemBcTypes;

        void reset()
        {
            isSet = false;
            auto resetPtr = [&] (auto&& uniquePtr) { uniquePtr.reset(); };
            std::for_each(bulkFvGeometries.begin(), bulkFvGeometries.end(), resetPtr);
            std::for_each(bulkElemVolVars.begin(), bulkElemVolVars.end(), resetPtr);
            std::for_each(bulkElemFluxVarsCache.begin(), bulkElemFluxVarsCache.end(), resetPtr);
            bulkLocalResidual.reset();
        }

        // resize containers for surface grids
        template<bool s = bulkIsSurfaceGrid, std::enable_if_t<s, int> = 0>
        void resize(std::size_t numEmbedments)
        {
            bulkFvGeometries.resize(numEmbedments, nullptr);
            bulkElemVolVars.resize(numEmbedments, nullptr);
            bulkElemFluxVarsCache.resize(numEmbedments, nullptr);
            bulkElemBcTypes.resize(numEmbedments);
        }

        // no memory reservation necessary for non-surface grids
        template<bool s = bulkIsSurfaceGrid, std::enable_if_t<!s, int> = 0>
        void resize(std::size_t numEmbedments)
        {}
    };

public:

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j = (i == bulkId) ? lowDimId : bulkId>
    using CouplingStencilType = typename CouplingMapper::template Stencil< CouplingMapper::template gridId<GridView<j>::dimension>() >;

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

        // set the sub problems
        this->setSubProblem(bulkProblem, bulkId);
        this->setSubProblem(lowDimProblem, lowDimId);

        // copy the solution vector
        ParentType::updateSolution(curSol);

        // determine all bulk elements that couple to low dim elements
        const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        bulkElemIsCoupled_.resize(bulkProblem->fvGridGeometry().gridView().size(0), false);
        std::for_each( bulkMap.begin(),
                       bulkMap.end(),
                       [&] (const auto& entry) { bulkElemIsCoupled_[entry.first] = true; });
    }

    /*!
     * \brief The coupling stencil of a given bulk domain element.
     */
    const CouplingStencilType<bulkId>& couplingStencil(BulkIdType domainI,
                                                       const Element<bulkId>& element,
                                                       LowDimIdType domainJ) const
    {
        const auto eIdx = this->problem(domainI).fvGridGeometry().elementMapper().index(element);

        if (bulkElemIsCoupled_[eIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            auto it = map.find(eIdx);
            assert(it != map.end());
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
        const auto eIdx = this->problem(lowDimId).fvGridGeometry().elementMapper().index(element);

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(eIdx);
        if (it != map.end()) return it->second.couplingStencil;
        else return getEmptyStencil(bulkId);
    }

    /*!
     * \brief returns true if a bulk scvf flux depends on data in the facet domain.
     * \note for box this is the case if an scvf lies on an interior boundary
     */
    bool isCoupled(const Element<bulkId>& element,
                   const SubControlVolumeFace<bulkId>& scvf) const
    { return scvf.isOnInteriorBoundary(); }

    /*!
     * \brief returns true if a bulk scvf coincides with a facet element.
     * \note for box, this is always true for coupled scvfs
     */
    bool isOnInteriorBoundary(const Element<bulkId>& element,
                              const SubControlVolumeFace<bulkId>& scvf) const
    { return isCoupled(element, scvf); }

    /*!
     * \brief returns the vol vars of a lower-dimensional element coinciding with a bulk scvf.
     */
    const VolumeVariables<lowDimId>& getLowDimVolVars(const Element<bulkId>& element,
                                                      const SubControlVolumeFace<bulkId>& scvf) const
    {
        const auto eIdx = this->problem(bulkId).fvGridGeometry().elementMapper().index(element);
        assert(bulkContext_.isSet);
        assert(bulkElemIsCoupled_[eIdx]);

        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        const auto& couplingData = map.find(eIdx)->second;

        // search the low dim element idx this scvf is embedded in
        // and find out the local index of the scv it couples to
        unsigned int coupledScvIdx;
        auto it = std::find_if( couplingData.elementToScvfMap.begin(),
                                couplingData.elementToScvfMap.end(),
                                [&] (auto& dataPair)
                                {
                                    const auto& scvfList = dataPair.second;
                                    auto it = std::find(scvfList.begin(), scvfList.end(), scvf.index());
                                    coupledScvIdx = std::distance(scvfList.begin(), it);
                                    return it != scvfList.end();
                                } );

        assert(it != couplingData.elementToScvfMap.end());
        const auto lowDimElemIdx = it->first;
        const auto& s = map.find(bulkContext_.elementIdx)->second.couplingElementStencil;
        const auto& idxInContext = std::distance( s.begin(), std::find(s.begin(), s.end(), lowDimElemIdx) );
        assert(std::find(s.begin(), s.end(), lowDimElemIdx) != s.end());

        if (FVGridGeometry<lowDimId>::discMethod == DiscretizationMethod::box)
            return bulkContext_.lowDimElemVolVars[idxInContext][coupledScvIdx];
        else
            return bulkContext_.lowDimElemVolVars[idxInContext][lowDimElemIdx];
    }

    /*!
     * \brief returns the lower-dimensional element coinciding with a bulk scvf.
     */
    const Element<lowDimId> getLowDimElement(const Element<bulkId>& element,
                                             const SubControlVolumeFace<bulkId>& scvf) const
    {
        const auto eIdx = this->problem(bulkId).fvGridGeometry().elementMapper().index(element);
        assert(bulkContext_.isSet);
        assert(bulkElemIsCoupled_[eIdx]);

        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        const auto& couplingData = map.find(eIdx)->second;

        // search the low dim element idx this scvf is embedded in
        auto it = std::find_if( couplingData.elementToScvfMap.begin(),
                                couplingData.elementToScvfMap.end(),
                                [&] (auto& dataPair)
                                {
                                    const auto& scvfList = dataPair.second;
                                    auto it = std::find(scvfList.begin(), scvfList.end(), scvf.index());
                                    return it != scvfList.end();
                                } );

        assert(it != couplingData.elementToScvfMap.end());
        const auto lowDimElemIdx = it->first;
        return this->problem(lowDimId).fvGridGeometry().element(lowDimElemIdx);
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
        const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);

        assert(bulkContext_.isSet);
        assert(bulkElemIsCoupled_[bulkContext_.elementIdx]);
        assert(map.find(bulkContext_.elementIdx) != map.end());
        assert(bulkContext_.elementIdx == this->problem(bulkId).fvGridGeometry().elementMapper().index(bulkLocalAssembler.element()));

        const auto& fvGeometry = bulkLocalAssembler.fvGeometry();
        typename LocalResidual<bulkId>::ElementResidualVector res(fvGeometry.numScv());
        res = 0.0;

        // compute fluxes across the coupling scvfs
        const auto& couplingScvfs = map.find(bulkContext_.elementIdx)->second.dofToCouplingScvfMap.at(dofIdxGlobalJ);
        for (auto scvfIdx : couplingScvfs)
        {
            const auto& scvf = fvGeometry.scvf(scvfIdx);
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            res[insideScv.localDofIndex()] += evalBulkFluxes_( bulkLocalAssembler.element(),
                                                               bulkLocalAssembler.fvGeometry(),
                                                               bulkLocalAssembler.curElemVolVars(),
                                                               bulkLocalAssembler.elemBcTypes(),
                                                               bulkLocalAssembler.elemFluxVarsCache(),
                                                               bulkLocalAssembler.localResidual(),
                                                               std::array<IndexType<bulkId>, 1>({scvfIdx}) );
        }

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
        assert(this->problem(lowDimId).fvGridGeometry().elementMapper().index(lowDimLocalAssembler.element()) == lowDimContext_.elementIdx);

        // evaluate sources for all scvs in lower-dimensional element
        typename LocalResidual<lowDimId>::ElementResidualVector res(lowDimLocalAssembler.fvGeometry().numScv());
        res = 0.0;
        for (const auto& scv : scvs(lowDimLocalAssembler.fvGeometry()))
            res[scv.localDofIndex()] -= evalSourcesFromBulk(lowDimLocalAssembler.element(),
                                                            lowDimLocalAssembler.fvGeometry(),
                                                            lowDimLocalAssembler.curElemVolVars(),
                                                            scv);
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
        // make sure the this is called for the element of the context
        assert(this->problem(lowDimId).fvGridGeometry().elementMapper().index(element) == lowDimContext_.elementIdx);

        NumEqVector<lowDimId> sources(0.0);
        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(lowDimContext_.elementIdx);
        if (it == map.end())
            return sources;

        assert(lowDimContext_.isSet);
        const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
        for (unsigned int i = 0; i < it->second.embedments.size(); ++i)
        {
            const auto& embedment = it->second.embedments[i];

            // list of scvfs in the bulk domain whose fluxes enter this scv
            // if low dim domain uses tpfa, this is all scvfs lying on this element
            // if it uses box, it is the one scvf coinciding with the given scv
            static constexpr bool lowDimUsesBox = FVGridGeometry<lowDimId>::discMethod == DiscretizationMethod::box;
            const auto& coincidingScvfs = bulkMap.find(embedment.first)->second.elementToScvfMap.at(lowDimContext_.elementIdx);
            const auto& scvfList = lowDimUsesBox ? std::vector<IndexType<lowDimId>>{ coincidingScvfs[scv.localDofIndex()] }
                                                 : coincidingScvfs;

            sources += evalBulkFluxes_(this->problem(bulkId).fvGridGeometry().element(embedment.first),
                                       *(lowDimContext_.bulkFvGeometries[i]),
                                       *(lowDimContext_.bulkElemVolVars[i]),
                                       lowDimContext_.bulkElemBcTypes[i],
                                       *(lowDimContext_.bulkElemFluxVarsCache[i]),
                                       *lowDimContext_.bulkLocalResidual,
                                       scvfList);
        }

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
        const auto bulkElemIdx = this->problem(bulkId).fvGridGeometry().elementMapper().index(element);
        bulkContext_.elementIdx = bulkElemIdx;

        // if element is coupled, actually set the context
        if (bulkElemIsCoupled_[bulkElemIdx])
        {
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);

            auto it = map.find(bulkElemIdx); assert(it != map.end());
            const auto& elementStencil = it->second.couplingElementStencil;
            bulkContext_.lowDimFvGeometries.reserve(elementStencil.size());
            bulkContext_.lowDimElemVolVars.reserve(elementStencil.size());

            for (const auto lowDimElemIdx : elementStencil)
            {
                const auto& ldGridGeometry = this->problem(lowDimId).fvGridGeometry();

                const auto elemJ = ldGridGeometry.element(lowDimElemIdx);
                auto fvGeom = localView(ldGridGeometry);
                auto elemVolVars = localView(assembler.gridVariables(lowDimId).curGridVolVars());
                fvGeom.bindElement(elemJ);
                elemVolVars.bindElement(elemJ, fvGeom, this->curSol()[lowDimId]);

                // TODO interpolated volvars
                bulkContext_.isSet = true;
                bulkContext_.lowDimFvGeometries.emplace_back( std::move(fvGeom) );
                bulkContext_.lowDimElemVolVars.emplace_back( std::move(elemVolVars) );
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
        const auto lowDimElemIdx = this->problem(lowDimId).fvGridGeometry().elementMapper().index(element);
        lowDimContext_.elementIdx = lowDimElemIdx;

        const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
        auto it = map.find(lowDimElemIdx);

        // if element is coupled, actually set the context
        if (it != map.end())
        {
            // first bind the low dim context for the first neighboring bulk element
            // this will set up the lower-dimensional data on this current lowdim element
            // which will be enough to compute the fluxes in all neighboring elements
            const auto& bulkGridGeom = this->problem(bulkId).fvGridGeometry();
            const auto bulkElem = bulkGridGeom.element(it->second.embedments[0].first);
            bindCouplingContext(bulkId, bulkElem, assembler);

            // reserve space in the context and bind local views of neighbors
            const auto& embedments = it->second.embedments;
            const auto numEmbedments = embedments.size();
            lowDimContext_.resize(numEmbedments);
            for (unsigned int i = 0; i < numEmbedments; ++i)
            {
                auto bulkFvGeom = localView(bulkGridGeom);
                auto bulkElemVolVars = localView(assembler.gridVariables(bulkId).curGridVolVars());
                auto bulkElemFluxVarsCache = localView(assembler.gridVariables(bulkId).gridFluxVarsCache());

                const auto curBulkElem = bulkGridGeom.element(embedments[i].first);
                bulkFvGeom.bind(curBulkElem);
                bulkElemVolVars.bind(curBulkElem, bulkFvGeom, this->curSol()[bulkId]);
                bulkElemFluxVarsCache.bind(curBulkElem, bulkFvGeom, bulkElemVolVars);

                lowDimContext_.isSet = true;
                lowDimContext_.bulkElemBcTypes[i].update(this->problem(bulkId), curBulkElem, bulkFvGeom);
                lowDimContext_.bulkFvGeometries[i] = std::make_unique< FVElementGeometry<bulkId> >( std::move(bulkFvGeom) );
                lowDimContext_.bulkElemVolVars[i] = std::make_unique< ElementVolumeVariables<bulkId> >( std::move(bulkElemVolVars) );
                lowDimContext_.bulkElemFluxVarsCache[i] = std::make_unique< ElementFluxVariablesCache<bulkId> >( std::move(bulkElemFluxVarsCache) );
            }

            // finally, set the local residual
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
            const auto& map = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            const auto& couplingElemStencil = map.find(bulkContext_.elementIdx)->second.couplingElementStencil;
            const auto& ldGridGeometry = this->problem(lowDimId).fvGridGeometry();

            // find the low-dim elements in coupling stencil, where this dof is contained in
            const auto couplingElements = [&] ()
            {
                if (FVGridGeometry<lowDimId>::discMethod == DiscretizationMethod::box)
                {
                    std::vector< Element<lowDimId> > lowDimElems;
                    std::for_each( couplingElemStencil.begin(), couplingElemStencil.end(),
                                   [&] (auto lowDimElemIdx)
                                   {
                                       auto element = ldGridGeometry.element(lowDimElemIdx);
                                       for (unsigned int i = 0; i < element.geometry().corners(); ++i)
                                       {
                                           const auto dofIdx = ldGridGeometry.vertexMapper().subIndex(element, i, lowDimDim);
                                           if (dofIdxGlobalJ == dofIdx) { lowDimElems.emplace_back( std::move(element) ); break; }
                                       }
                                   } );
                    return lowDimElems;
                }
                // dof index = element index for cc schemes
                else
                    return std::vector<Element<lowDimId>>( {ldGridGeometry.element(dofIdxGlobalJ)} );
            } ();

            // update all necessary vol vars in context
            for (const auto& element : couplingElements)
            {
                // find index in coupling context
                const auto eIdxGlobal = ldGridGeometry.elementMapper().index(element);
                auto it = std::find(couplingElemStencil.begin(), couplingElemStencil.end(), eIdxGlobal);
                const auto idxInContext = std::distance(couplingElemStencil.begin(), it);
                assert(it != couplingElemStencil.end());

                auto& elemVolVars = bulkContext_.lowDimElemVolVars[idxInContext];
                const auto& fvGeom = bulkContext_.lowDimFvGeometries[idxInContext];

                // TODO interpolated vol vars
                elemVolVars.bindElement(element, fvGeom, this->curSol()[lowDimId]);
            }
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
     *        the element volume variables of the neighboring bulk elements stored
     *        in the context.
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
            assert(lowDimContext_.elementIdx == this->problem(lowDimId).fvGridGeometry().elementMapper().index(lowDimLocalAssembler.element()));

            const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
            auto it = map.find(lowDimContext_.elementIdx);

            assert(it != map.end());
            const auto& embedments = it->second.embedments;
            const auto& bulkGridGeometry = this->problem(bulkId).fvGridGeometry();

            // TODO more efficient (only one dof update per embedment)
            // update the elem volvars in context of those elements that contain globalJ
            unsigned int embedmentIdx = 0;
            for (const auto& embedment : embedments)
            {
                const auto elementJ = bulkGridGeometry.element(embedment.first);

                for (unsigned int i = 0; i < elementJ.subEntities(bulkDim); ++i)
                {
                    const auto dofIdx = bulkGridGeometry.vertexMapper().subIndex(elementJ, i, bulkDim);
                    if (dofIdx == dofIdxGlobalJ)
                    {
                        // element contains the deflected dof
                        const auto& fvGeom = *lowDimContext_.bulkFvGeometries[embedmentIdx];
                        (*lowDimContext_.bulkElemVolVars[embedmentIdx]).bindElement(elementJ, fvGeom, this->curSol()[bulkId]);
                    }
                }

                // keep track of embedments
                embedmentIdx++;
            }
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
            assert(bulkContext_.isSet);
            assert(lowDimContext_.elementIdx == this->problem(lowDimId).fvGridGeometry().elementMapper().index(lowDimLocalAssembler.element()));

            // update the corresponding vol vars in the bulk context
            const auto& bulkMap = couplingMapperPtr_->couplingMap(bulkGridId, lowDimGridId);
            const auto& couplingElementStencil = bulkMap.find(bulkContext_.elementIdx)->second.couplingElementStencil;
            auto it = std::find(couplingElementStencil.begin(), couplingElementStencil.end(), lowDimContext_.elementIdx);
            assert(it != couplingElementStencil.end());
            const auto idxInContext = std::distance(couplingElementStencil.begin(), it);

            auto& elemVolVars = bulkContext_.lowDimElemVolVars[idxInContext];
            const auto& fvGeom = bulkContext_.lowDimFvGeometries[idxInContext];
            const auto& element = lowDimLocalAssembler.element();

            // TODO interpolated update
            elemVolVars.bindElement(element, fvGeom, this->curSol()[lowDimId]);
        }
    }

    //! Empty stencil to be returned for elements that aren't coupled
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    const typename CouplingMapper::template Stencil<id>&
    getEmptyStencil(Dune::index_constant<id>) const
    { return std::get<(id == bulkId ? 0 : 1)>(emptyStencilTuple_); }

private:
    //! evaluates the bulk-facet exchange fluxes for a given facet element
    template<class BulkScvfIndices>
    NumEqVector<bulkId> evalBulkFluxes_(const Element<bulkId>& elementI,
                                        const FVElementGeometry<bulkId>& fvGeometry,
                                        const ElementVolumeVariables<bulkId>& elemVolVars,
                                        const ElementBoundaryTypes<bulkId>& elemBcTypes,
                                        const ElementFluxVariablesCache<bulkId>& elemFluxVarsCache,
                                        const LocalResidual<bulkId>& localResidual,
                                        const BulkScvfIndices& scvfIndices) const
    {

        NumEqVector<bulkId> coupledFluxes(0.0);
        for (const auto& scvfIdx : scvfIndices)
            coupledFluxes += localResidual.evalFlux(this->problem(bulkId),
                                                    elementI,
                                                    fvGeometry,
                                                    elemVolVars,
                                                    elemBcTypes,
                                                    elemFluxVarsCache,
                                                    fvGeometry.scvf(scvfIdx));
        return coupledFluxes;
    }

    std::shared_ptr<CouplingMapper> couplingMapperPtr_;

    //! store bools for all bulk elements that indicate if they
    //! are coupled, so that we don't have to search in the map every time
    std::vector<bool> bulkElemIsCoupled_;

    //! empty stencil to return for non-coupled elements
    using BulkStencil = typename CouplingMapper::template Stencil<bulkId>;
    using LowDimStencil = typename CouplingMapper::template Stencil<lowDimId>;
    std::tuple<BulkStencil, LowDimStencil> emptyStencilTuple_;

    //! The coupling contexts of the two domains
    BulkCouplingContext bulkContext_;
    LowDimCouplingContext lowDimContext_;
};

} // end namespace Dumux

#endif