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
 * \ingroup OnePNCTests
 * \brief The spatial params the one-phase two-components facet coupling test
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTests
 * \brief The spatial params the single-phase two-components facet coupling test
 */
template< class FVGridGeometry, class Scalar >
class OnePSpatialParams
: public FVSpatialParamsOneP< FVGridGeometry, Scalar, OnePSpatialParams<FVGridGeometry, Scalar> >
{
    using ThisType = OnePSpatialParams< FVGridGeometry, Scalar >;
    using ParentType = FVSpatialParamsOneP< FVGridGeometry, Scalar, ThisType >;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! export the type used for permeabilities
    using PermeabilityType = Scalar;

    //! the constructor
    OnePSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    //! Return the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

private:
    PermeabilityType permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
