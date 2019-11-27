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

// the header guard
#ifndef DUMUX_BIOSUB_HH
#define DUMUX_BIOSUB_HH

// including the base component
#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

// In Substrate, we define the properties of the component substrate, which is the energy and carbon source for biomass growth
template <class Scalar>
class Substrate
: public Components::Base<Scalar, Substrate<Scalar> >
{
public:
    // the name
    static std::string name()
    { return "Substrate"; }

    // the molar mass, using the molar mass of glucose
    static Scalar molarMass()
    { return 0.18016; } // kg/mol
};

} // end namespace Components
} // end namespace Dumux

#endif