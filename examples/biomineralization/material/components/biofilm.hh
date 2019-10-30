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
 * \ingroup Components
 * \brief A class for the biofilm solid component properties.
 */
#ifndef DUMUX_BIOFILM_HH
#define DUMUX_BIOFILM_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the biofilm solid component properties.
 */
template <class Scalar>
class Biofilm
: public Components::Base<Scalar, Biofilm<Scalar> >
, public Components::Solid<Scalar, Biofilm<Scalar> >
{
public:
   /*!
    * \brief A human readable name for the biofilm.
    */
    static std::string name()
    { return "Biofilm"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular biofilm.
    * This can be set in the input file, as the correct molar mass of a biofilm is hard to define and
    * any choice is arbitrary. Based on a cell mass of 2.5e-16 kg, the molar mass of "cells" would
    * be 1.5e8 kg/mol, but biofilms are more than just cells ...
    */
    static Scalar molarMass()
    {
        Scalar molarMass = getParam<Scalar>("BioCoefficients.BiofilmMolarMass", 1);
        return molarMass;
    }

   /*!
    * \brief The (dry) density \f$\mathrm{[kg/m^3]}\f$ of biofilm.
    * This is set in the input file, as biofilm densities vary over orders of magnitude.
    * A value of "10" is fitted by Ebigbo 2012 WRR
    * Hommel 2015 WRR uses "6.9".
    */
    static Scalar solidDensity(Scalar temperature)
    {
        Scalar rho = getParam<Scalar>("BioCoefficients.RhoBiofilm", 10);
        return rho;
    }
};

} // end namespace Components
} // end namespace Dumux

#endif