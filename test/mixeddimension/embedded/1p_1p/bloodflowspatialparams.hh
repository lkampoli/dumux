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
 * \ingroup OneTests
 * \brief The spatial parameters class blood flow problem
 */
#ifndef DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH
#define DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH

#include <dune/common/fvector.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OneTests
 * \brief Definition of the spatial parameters for the blood flow problem
 */
template<class TypeTag>
class BloodFlowSpatialParams: public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dimworld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        radius_ = getParam<Scalar>("SpatialParams.Radius");
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param ipGlobal The integration point
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& ipGlobal) const
    {
        return (1 + ipGlobal[2] + 0.5*ipGlobal[2]*ipGlobal[2])/(M_PI*radius(0)*radius(0));
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param the index of the element
     */
    Scalar radius(unsigned int eIdxGlobal) const
    {
        return radius_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 1.0;
    }

private:
    Scalar radius_;
};

} // end namespace Dumux

#endif
