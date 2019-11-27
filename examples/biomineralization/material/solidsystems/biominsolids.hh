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
#ifndef DUMUX_SOLIDSYSTEMS_MICP_SOLID_PHASE_HH
#define DUMUX_SOLIDSYSTEMS_MICP_SOLID_PHASE_HH

#include <string>
#include <dune/common/exceptions.hh>

// we include all necessary solid components
#include <examples/biomineralization/material/components/biofilm.hh>
#include <dumux/material/components/calcite.hh>
#include <dumux/material/components/granite.hh>

// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
namespace Dumux {
namespace SolidSystems {

// In the BioMinSolidPhase solid system, we define all functions needed to describe the solids accounted for in our simulation
template <class Scalar>
class BioMinSolidPhase
{
public:

    // We use convenient declarations that we derive from the property system.
    using Biofilm = Components::Biofilm<Scalar>;
    using Calcite = Components::Calcite<Scalar>;
    using Granite = Components::Granite<Scalar>;

    /****************************************
     * Solid phase related static parameters
     ****************************************/
    static constexpr int numComponents = 3;
    static constexpr int numInertComponents = 1;
    static constexpr int BiofilmIdx = 0;
    static constexpr int CalciteIdx = 1;
    static constexpr int GraniteIdx = 2;

    // The component names
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case BiofilmIdx: return Biofilm::name();
            case CalciteIdx: return Calcite::name();
            case GraniteIdx: return Granite::name();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    // The solid system's name
    static std::string name()
    { return "EnzymeMinSolidPhase"; }

    // We assume incompressible solids
    static constexpr bool isCompressible(int compIdx)
    { return false; }

    // we have one inert component, the others may change
    static constexpr bool isInert()
    { return (numComponents == numInertComponents); }

    // The component molar masses
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case BiofilmIdx: return Biofilm::molarMass();
            case CalciteIdx: return Calcite::molarMass();
            case GraniteIdx: return Granite::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    // The average density
    template <class SolidState>
    static Scalar density(const SolidState& solidState)
    {
        Scalar rho1 = Biofilm::solidDensity(solidState.temperature());
        Scalar rho2 = Calcite::solidDensity(solidState.temperature());
        Scalar rho3 = Granite::solidDensity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(BiofilmIdx);
        Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (rho1*volFrac1+
                rho2*volFrac2+
                rho3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }

    // The component densities
    template <class SolidState>
    static Scalar density(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case BiofilmIdx: return Biofilm::solidDensity(solidState.temperature());
            case CalciteIdx: return Calcite::solidDensity(solidState.temperature());
            case GraniteIdx: return Granite::solidDensity(solidState.temperature());
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }

    // The component molar densities
    template <class SolidState>
    static Scalar molarDensity(const SolidState& solidState, const int compIdx)
    {
        switch (compIdx)
        {
            case BiofilmIdx: return Biofilm::solidDensity(solidState.temperature())/Biofilm::molarMass();
            case CalciteIdx: return Calcite::solidDensity(solidState.temperature())/Calcite::molarMass();
            case GraniteIdx: return Granite::solidDensity(solidState.temperature())/Calcite::molarMass();
            default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }


    // The average thermal conductivity
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        Scalar lambda1 = Biofilm::solidThermalConductivity(solidState.temperature());
        Scalar lambda2 = Calcite::solidThermalConductivity(solidState.temperature());
        Scalar lambda3 = Granite::solidThermalConductivity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(BiofilmIdx);
        Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (lambda1*volFrac1+
                lambda2*volFrac2+
                lambda3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }


    // The average heat capacity
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        Scalar c1 = Biofilm::solidHeatCapacity(solidState.temperature());
        Scalar c2 = Calcite::solidHeatCapacity(solidState.temperature());
        Scalar c3 = Granite::solidHeatCapacity(solidState.temperature());
        Scalar volFrac1 = solidState.volumeFraction(BiofilmIdx);
        Scalar volFrac2 = solidState.volumeFraction(CalciteIdx);
        Scalar volFrac3 = solidState.volumeFraction(GraniteIdx);

        return (c1*volFrac1+
                c2*volFrac2+
                c3*volFrac3)
               /(volFrac1+volFrac2+volFrac3);
    }

};

} // end namespace SolidSystems
} // end namespace Dumux

#endif
