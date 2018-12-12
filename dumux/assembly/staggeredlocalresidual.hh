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
 * \ingroup Assembly
 * \brief Calculates the element-wise residual for the staggered FV scheme
 */
#ifndef DUMUX_STAGGERED_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_LOCAL_RESIDUAL_HH

#include <dumux/common/timeloop.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief Calculates the element-wise residual for the staggered FV scheme
 */
template<class TypeTag>
class StaggeredLocalResidual
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    using CellCenterResidual = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidual = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;

    using TimeLoop = TimeLoopBase<Scalar>;

public:
    using CellCenterResidualValue = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceResidualValue = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using ElementResidualVector = CellCenterResidualValue;

    //! the constructor
    StaggeredLocalResidual(const Problem* problem,
                           const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
    {}

    //! Convenience function to evaluate the flux and source terms for the cell center residual
    CellCenterResidualValue evalFluxAndSourceForCellCenter(const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const ElementFaceVariables& elemFaceVars,
                                                           const ElementBoundaryTypes& bcTypes,
                                                           const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        CellCenterResidualValue residual(0.0);

        // evaluate the source term
        for (auto&& scv : scvs(fvGeometry))
            asImp().evalSourceForCellCenter(residual, this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, scv);

        // evaluate the flux term
        for (auto&& scvf : scvfs(fvGeometry))
            asImp().evalFluxForCellCenter(residual, this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, scvf);

        return residual;
    }

    //! Evaluate the flux terms for a cell center residual
    void evalFluxForCellCenter(CellCenterResidualValue& residual,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const ElementFaceVariables& elemFaceVars,
                               const ElementBoundaryTypes& bcTypes,
                               const ElementFluxVariablesCache& elemFluxVarsCache,
                               const SubControlVolumeFace& scvf) const
    {
        if(!scvf.boundary())
            residual += asImp_().computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);
    }

    //! Evaluate the source terms for a cell center residual
    void evalSourceForCellCenter(CellCenterResidualValue& residual,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFaceVariables& curElemFaceVars,
                                 const SubControlVolume& scv) const
    {
            const auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();

            // subtract the source term from the local rate
            auto source = asImp_().computeSourceForCellCenter(problem, element, fvGeometry, curElemVolVars, curElemFaceVars, scv);
            source *= scv.volume()*curExtrusionFactor;
            residual -= source;
    }

    //! Evaluate the storage terms for a cell center residual
    CellCenterResidualValue evalStorageForCellCenter(const Element &element,
                                                     const FVElementGeometry& fvGeometry,
                                                     const ElementVolumeVariables& prevPrevElemVolVars,
                                                     const ElementVolumeVariables& prevElemVolVars,
                                                     const ElementVolumeVariables& curElemVolVars) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        CellCenterResidualValue storage(0.0);

        for (auto&& scv : scvs(fvGeometry))
            asImp().evalStorageForCellCenter(storage, problem(), element, fvGeometry, prevPrevElemVolVars, prevElemVolVars, curElemVolVars, scv);

        return storage;
    }

    //! Evaluate the storage terms for a cell center residual
    void evalStorageForCellCenter(CellCenterResidualValue& residual,
                                  const Problem& problem,
                                  const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& prevPrevElemVolVars,
                                  const ElementVolumeVariables& prevElemVolVars,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const SubControlVolume& scv) const
    {
        CellCenterResidualValue storage(0.0);
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // mass balance within the element. this is the
        // \f$\frac{m}{\partial t}\f$ term if using implicit
        // euler as time discretization.

        // We might need a more explicit way for
        // doing the time discretization...
        auto prevCCStorage = asImp_().computeStorageForCellCenter(problem, scv, prevVolVars);
        auto curCCStorage = asImp_().computeStorageForCellCenter(problem, scv, curVolVars);

        prevCCStorage *= prevVolVars.extrusionFactor();
        curCCStorage *= curVolVars.extrusionFactor();

        if (timeLoop_->usingBdf2())
        {
            // BDF2 method
            const auto& prevPrevVolVars = prevPrevElemVolVars[scv];
            auto prevPrevCCStorage = asImp_().computeStorageForCellCenter(problem, scv, prevPrevVolVars);
            prevPrevCCStorage *= prevPrevVolVars.extrusionFactor();

            curCCStorage *= curFactor_();
            prevCCStorage *= prevFactor_();
            prevPrevCCStorage *= prevPrevFactor_();
            storage = std::move(curCCStorage);
            storage -= std::move(prevCCStorage);
            storage += std::move(prevPrevCCStorage);
        }
        else
        {
            // Implicit Euler method
            storage = std::move(curCCStorage);
            storage -= std::move(prevCCStorage);
            storage /= timeLoop_->timeStepSize();
        }

        storage *= scv.volume();
        residual += storage;
    }

    //! Evaluate the boundary conditions for a cell center residual
    void evalBoundaryForCellCenter(CellCenterResidualValue& residual,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const ElementFaceVariables& elemFaceVars,
                                   const ElementBoundaryTypes& bcTypes,
                                   const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        asImp_().evalBoundaryForCellCenter_(residual, problem, element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache);
    }

    //! for compatibility with FVLocalAssemblerBase
    template<class... Args>
    CellCenterResidualValue evalFluxAndSource(Args&&... args) const
    {
        return CellCenterResidualValue(0.0);
    }

    //! for compatibility with FVLocalAssemblerBase
    template<class... Args>
    CellCenterResidualValue evalStorage(Args&&... args) const
    {
        return CellCenterResidualValue(0.0);
    }

    /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting residual information.
     */
    // \{

    //! Convenience function to evaluate the flux and source terms for the face residual
    FaceResidualValue evalFluxAndSourceForFace(const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const ElementBoundaryTypes& bcTypes,
                                               const ElementFluxVariablesCache& elemFluxVarsCache,
                                               const SubControlVolumeFace& scvf) const
    {
        FaceResidualValue residual(0.0);
        asImp().evalSourceForFace(residual, this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        asImp().evalFluxForFace(residual, this->problem(), element, fvGeometry, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache, scvf);

        return residual;
    }

    //! Evaluate the flux terms for a face residual
    void evalFluxForFace(FaceResidualValue& residual,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFaceVariables& elemFaceVars,
                         const ElementBoundaryTypes& bcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        if(!scvf.boundary())
            residual += asImp_().computeFluxForFace(problem, element, scvf, fvGeometry, elemVolVars, elemFaceVars, elemFluxVarsCache);
    }

    //! Evaluate the source terms for a face residual
    void evalSourceForFace(FaceResidualValue& residual,
                           const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFaceVariables& elemFaceVars,
                           const SubControlVolumeFace& scvf) const
    {
        // the source term:
        auto source = asImp_().computeSourceForFace(problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto extrusionFactor = elemVolVars[scv].extrusionFactor();

        // multiply by 0.5 because we only consider half of a staggered control volume here
        source *= 0.5*scv.volume()*extrusionFactor;
        residual -= source;
    }

    //! Evaluate the storage terms for a face residual
    FaceResidualValue evalStorageForFace(const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& prevPrevElemVolVars,
                                         const ElementVolumeVariables& prevElemVolVars,
                                         const ElementVolumeVariables& curElemVolVars,
                                         const ElementFaceVariables& prevPrevElemFaceVars,
                                         const ElementFaceVariables& prevElemFaceVars,
                                         const ElementFaceVariables& curElemFaceVars,
                                         const SubControlVolumeFace& scvf) const
    {
        assert(timeLoop_ && "no time loop set for storage term evaluation");
        FaceResidualValue storage(0.0);
        asImp().evalStorageForFace(storage, problem(), element, fvGeometry, prevPrevElemVolVars, prevElemVolVars, curElemVolVars, prevPrevElemFaceVars, prevElemFaceVars, curElemFaceVars, scvf);
        return storage;
    }

    //! Evaluate the storage terms for a face residual
    void evalStorageForFace(FaceResidualValue& residual,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& prevPrevElemVolVars,
                            const ElementVolumeVariables& prevElemVolVars,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFaceVariables& prevPrevElemFaceVars,
                            const ElementFaceVariables& prevElemFaceVars,
                            const ElementFaceVariables& curElemFaceVars,
                            const SubControlVolumeFace& scvf) const
    {
        FaceResidualValue storage(0.0);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        auto prevFaceStorage = asImp_().computeStorageForFace(problem, scvf, prevElemVolVars[scv], prevElemFaceVars);
        auto curFaceStorage = asImp_().computeStorageForFace(problem, scvf, curElemVolVars[scv], curElemFaceVars);

        if (timeLoop_->usingBdf2())
        {
            // BDF2 method
            auto prevPrevFaceStorage = asImp_().computeStorageForFace(problem, scvf, prevPrevElemVolVars[scv], prevPrevElemFaceVars);

            curFaceStorage *= curFactor_();
            prevFaceStorage *= prevFactor_();
            prevPrevFaceStorage *= prevPrevFactor_();
            storage = std::move(curFaceStorage);
            storage -= std::move(prevFaceStorage);
            storage += std::move(prevPrevFaceStorage);
        }
        else
        {
            // Implicit Euler method
            storage = std::move(curFaceStorage);
            storage -= std::move(prevFaceStorage);
            storage /= timeLoop_->timeStepSize();
        }

        const auto extrusionFactor = curElemVolVars[scv].extrusionFactor();

        // multiply by 0.5 because we only consider half of a staggered control volume here
        storage *= 0.5*scv.volume()*extrusionFactor;
        residual += storage;
    }

    //! Evaluate the boundary conditions for a face residual
    void evalBoundaryForFace(FaceResidualValue& residual,
                             const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const ElementFaceVariables& elemFaceVars,
                             const ElementBoundaryTypes& bcTypes,
                             const ElementFluxVariablesCache& elemFluxVarsCache,
                             const SubControlVolumeFace& scvf) const
    {
        asImp_().evalBoundaryForFace_(residual, problem, element, fvGeometry, scvf, elemVolVars, elemFaceVars, bcTypes, elemFluxVarsCache);
    }

    //! If no solution has been set, we treat the problem as stationary.
    bool isStationary() const
    { return !timeLoop_; }

    //! the problem
    const Problem& problem() const
    { return *problem_; }

protected:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }


    TimeLoop& timeLoop()
    { return *timeLoop_; }

    const TimeLoop& timeLoop() const
    { return *timeLoop_; }

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_;

    Scalar curFactor_() const
    {
        const Scalar& dt = timeLoop_->timeStepSize();
        const Scalar& dtOld = timeLoop_->previousTimeStepSize();
        return 1.0 / dt + 1.0 / (dt + dtOld);
    }

    Scalar prevFactor_() const
    {
        const Scalar& dt = timeLoop_->timeStepSize();
        const Scalar& dtOld = timeLoop_->previousTimeStepSize();
        return (dt + dtOld) / dt / dtOld;
    }

    Scalar prevPrevFactor_() const
    {
        const Scalar& dt = timeLoop_->timeStepSize();
        const Scalar& dtOld = timeLoop_->previousTimeStepSize();
        return dt / dtOld / (dt + dtOld);
    }

};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
