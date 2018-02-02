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
 * \ingroup NIModel
 * \brief The implicit non-isothermal model.
 *
 * This model implements a generic energy balance for single and multi-phase
 * transport problems. Currently the non-isothermal model can be used on top of
 * the 1p2c, 2p, 2p2c and 3p3c models. Comparison to simple analytical solutions
 * for pure convective and conductive problems are found in the 1p2c test. Also refer
 * to this test for details on how to activate the non-isothermal model.
 *
 * For the energy balance, local thermal equilibrium is assumed. This
 * results in one energy conservation equation for the porous solid
 * matrix and the fluids:
 \f{align*}{
 \phi \frac{\partial \sum_\alpha \varrho_\alpha u_\alpha S_\alpha}{\partial t}
 & +
 \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 -
 \sum_\alpha \text{div}
 \left\{
 \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left( \textbf{grad}\,p_\alpha - \varrho_\alpha \mbox{\bf g} \right)
 \right\} \\
    & - \text{div} \left(\lambda_{pm} \textbf{grad} \, T \right)
    - q^h = 0.
 \f}
 * where \f$h_\alpha\f$ is the specific enthalpy of a fluid phase
 * \f$\alpha\f$ and \f$u_\alpha = h_\alpha -
 * p_\alpha/\varrho_\alpha\f$ is the specific internal energy of the
 * phase.
 */

#ifndef DUMUX_NONISOTHERMAL_MODEL_HH
#define DUMUX_NONISOTHERMAL_MODEL_HH

#include <dumux/common/properties.hh>
#include "indices.hh"
#include "vtkoutputfields.hh"

namespace Dumux {
namespace Properties {

NEW_TYPE_TAG(NonIsothermal);

SET_BOOL_PROP(NonIsothermal, EnableEnergyBalance, true);

//! add the energy balance
SET_INT_PROP(NonIsothermal, NumEq, GET_PROP_VALUE(TypeTag, IsothermalNumEq) + 1);

//! indices for non-isothermal models
SET_TYPE_PROP(NonIsothermal, Indices, EnergyIndices<TypeTag, 0>);

//! indices for non-isothermal models
SET_TYPE_PROP(NonIsothermal, VtkOutputFields, EnergyVtkOutputFields<TypeTag>);

} // end namespace Properties
} // end namespace Dumux

#endif