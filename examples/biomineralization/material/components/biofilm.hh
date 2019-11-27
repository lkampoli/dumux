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
#ifndef DUMUX_BIOFILM_HH
#define DUMUX_BIOFILM_HH

// including the base and the generic solid component
#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

// In Biofilm, we define the properties of the solid component biofilm
template <class Scalar>
class Biofilm
: public Components::Base<Scalar, Biofilm<Scalar> >
, public Components::Solid<Scalar, Biofilm<Scalar> >
{
public:
    // the name
    static std::string name()
    { return "Biofilm"; }

    // the molar mass
    static Scalar molarMass()
    {
        Scalar molarMass = getParam<Scalar>("BioCoefficients.BiofilmMolarMass", 1);
        return molarMass;
    }

    // the density
    static Scalar solidDensity(Scalar temperature)
    {
        Scalar rho = getParam<Scalar>("BioCoefficients.RhoBiofilm", 10);
        return rho;
    }
};

} // end namespace Components
} // end namespace Dumux

#endif