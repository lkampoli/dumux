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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_VOLUMEVARIABLES_HH

#include <dune/common/exceptions.hh>
#include <dune/common/rangeutilities.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/staggered/freeflow/elementvolumevariables.hh>

namespace Dumux {

template<class P, class VV>
struct StaggeredGridDefaultGridVolumeVariablesTraits
{
    using Problem = P;
    using VolumeVariables = VV;
    using PrimaryVariables = typename VV::PrimaryVariables;

    template<class GridVolumeVariables, bool cachingEnabled>
    using LocalView = StaggeredElementVolumeVariables<GridVolumeVariables, cachingEnabled>;

    //! Use the PassKey pattern to restrict access to a member function only to a specific class
    template<class T>
    class Key { friend T; Key() {} Key(Key const&) {} };

    //! Returns the primary variables used for the boundary volVars and checks for admissible
    //! combinations for boundary conditions.
    template<class Problem, class SolutionVector, class Element, class SubControlVolumeFace>
    static PrimaryVariables getBoundaryPriVars(const Problem& problem,
                                               const SolutionVector& sol,
                                               const Element& element,
                                               const SubControlVolumeFace& scvf)
    {
        using CellCenterPrimaryVariables = typename SolutionVector::value_type;
        using Indices = typename VolumeVariables::Indices;
        static constexpr auto dim = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
        static constexpr auto offset = dim;

        const auto bcTypes = problem.boundaryTypes(element, scvf);
        PrimaryVariables boundaryPriVars(0.0);

        // make sure to not use outflow BC for momentum balance
        for(int i = 0; i < dim; ++i)
        {
            if(bcTypes.isOutflow(Indices::velocity(i)))
                DUNE_THROW(Dune::InvalidStateException, "Outflow condition cannot be used for velocity. Set only a Dirichlet value for pressure instead.");
        }

        if(bcTypes.isOutflow(Indices::pressureIdx))
            DUNE_THROW(Dune::InvalidStateException, "Outflow condition cannot be used for pressure. Set only a Dirichlet value for velocity instead.");

        // Determine the pressure value at a boundary with a Dirichlet condition for velocity.
        // This just takes the value of the adjacent inner cell.
        if(bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
        {
            if(bcTypes.isDirichlet(Indices::pressureIdx))
                DUNE_THROW(Dune::InvalidStateException, "A Dirichlet condition for velocity must not be combined with a Dirichlet condition for pressure");
            else
                boundaryPriVars[Indices::pressureIdx] = sol[scvf.insideScvIdx()][Indices::pressureIdx - offset];
                // TODO: pressure could be extrapolated to the boundary
        }

        // Determine the pressure value for a boundary with a Dirichlet condition for pressure.
        // Takes a value specified in the problem.
        if(bcTypes.isDirichlet(Indices::pressureIdx))
        {
            if(bcTypes.isDirichlet(Indices::velocity(scvf.directionIndex())))
                DUNE_THROW(Dune::InvalidStateException, "A Dirichlet condition for velocity must not be combined with a Dirichlet condition for pressure");
            else
                boundaryPriVars[Indices::pressureIdx] = problem.dirichlet(element, scvf)[Indices::pressureIdx];
        }

        // Return for isothermal single-phase systems ...
        if(CellCenterPrimaryVariables::dimension == 1)
            return boundaryPriVars;

        // ... or handle values for components, temperature, etc.
        for(int eqIdx = offset; eqIdx < PrimaryVariables::dimension; ++eqIdx)
        {
            if(eqIdx == Indices::pressureIdx)
                continue;

            if(bcTypes.isDirichlet(eqIdx))
                boundaryPriVars[eqIdx] = problem.dirichlet(element, scvf)[eqIdx];
            else if(bcTypes.isOutflow(eqIdx) || bcTypes.isSymmetry() || bcTypes.isNeumann(eqIdx))
                boundaryPriVars[eqIdx] = sol[scvf.insideScvIdx()][eqIdx - offset];
        }

        // make sure that a potential outflow condition is set for all components
        std::array<bool, VolumeVariables::numFluidComponents() - 1> isComponentOutflow;
        for(int compIdx = 1; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
        {
            const auto eqIdx = VolumeVariables::Indices::conti0EqIdx + compIdx;
            isComponentOutflow[compIdx -1] = bcTypes.isOutflow(eqIdx);
        }

        if(Dune::any_true(isComponentOutflow) && !Dune::all_true(isComponentOutflow))
            DUNE_THROW(Dune::InvalidStateException, "Outflow condition must be set for all components!");

        return boundaryPriVars;
    }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models
 */
template<class Traits, bool cachingEnabled>
class StaggeredGridVolumeVariables;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models.
          Specialization in case of storing the volume variables
 */
template<class Traits>
class StaggeredGridVolumeVariables<Traits, /*cachingEnabled*/true>
{
    using ThisType = StaggeredGridVolumeVariables<Traits, true>;
    using Problem = typename Traits::Problem;
    using PrimaryVariables = typename Traits::VolumeVariables::PrimaryVariables;

public:
    //! export the type of the indices
    using Indices = typename Traits::VolumeVariables::Indices;

    //! export the type of the VolumeVariables
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    StaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    //! Update all volume variables
    template<class FVGridGeometry, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol)
    {
        auto numScv = fvGridGeometry.numScv();
        auto numBoundaryScvf = fvGridGeometry.numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = sol[scv.dofIndex()];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

                auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));
                volumeVariables_[scv.dofIndex()].update(elemSol, problem(), element, scv);
            }
        }
    }

    //! Update the boundary volume variables.
    //! Use the PassKey pattern to grant access only to the LocalView (ElementVolumeVariables) class.
    template<class Element, class FVGeometry, class SolutionVector>
    void updateBoundary_(const Element& element, const FVGeometry& fvGeometry, const SolutionVector& sol,
                         typename Traits::template Key<LocalView>) const
    {
        // handle the boundary volume variables
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // if we are not on a boundary, skip the rest
            if (!scvf.boundary())
                continue;

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            auto boundaryPriVars = Traits::getBoundaryPriVars(problem(), sol, element, scvf);
            const auto elemSol = elementSolution<FVGeometry>(std::move(boundaryPriVars));
            volumeVariables_[scvf.outsideScvIdx()].update(elemSol, problem(), element, insideScv);
        }
    }

    const VolumeVariables& volVars(const std::size_t scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const std::size_t scvIdx)
    { return volumeVariables_[scvIdx]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volumeVariables_[scv.dofIndex()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& volVars(const SubControlVolume& scv)
    { return volumeVariables_[scv.dofIndex()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    mutable std::vector<VolumeVariables> volumeVariables_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Grid volume variables class for staggered models.
          Specialization in case of not storing the volume variables
 */
template<class Traits>
class StaggeredGridVolumeVariables<Traits, /*cachingEnabled*/false>
{
    using ThisType = StaggeredGridVolumeVariables<Traits, false>;
    using Problem = typename Traits::Problem;
    using PrimaryVariables = typename Traits::VolumeVariables::PrimaryVariables;

public:
    //! export the type of the VolumeVariables
    using VolumeVariables = typename Traits::VolumeVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    StaggeredGridVolumeVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class FVGridGeometry, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

    //! Returns the primary variables used for the boundary volVars and checks for admissible
    //! combinations for boundary conditions.
    template<class... Args>
    PrimaryVariables getBoundaryPriVars(Args&&... args) const
    {
        return Traits::getBoundaryPriVars(std::forward<Args>(args)...);
    }

private:

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
