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
 *
 * \brief The spatial parameters class for the matrix problem
 */
#ifndef DUMUX_1P_HETEROGENEOUS_MATRIX_SPATIALPARAMS_HH
#define DUMUX_1P_HETEROGENEOUS_MATRIX_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the matrix problem
 */
template<class TypeTag>
class OnePMatrixSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

public:
    using PermeabilityType = Scalar;

    OnePMatrixSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {}

    /*!
     * \brief Return the intrinsic permeability for the given position in [m^2].
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar k = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPermeability);
        return k;
    }

    /*!
     * \brief Define the porosity in [-].
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }
};
} //end namespace

#endif