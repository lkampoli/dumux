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
#ifndef DUMUX_BIOFLUID_SIMPLE_CHEM_SYSTEM_HH
#define DUMUX_BIOFLUID_SIMPLE_CHEM_SYSTEM_HH

#include <dumux/material/idealgas.hh>

// we include the base fluid system
#include <dumux/material/fluidsystems/base.hh>
// we include the brine adapter adapting component indices to use the brine fluidsystem
#include "icpcomplexsalinitybrine.hh"

// we include all necessary fluid components
#include <dumux/material/fluidstates/adapter.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/co2tablereader.hh>
#include <examples/biomineralization/material/components/sodiumion.hh>
#include <examples/biomineralization/material/components/chlorideion.hh>
#include <dumux/material/components/calciumion.hh>
#include <examples/biomineralization/material/components/urea.hh>
#include <dumux/material/components/o2.hh>
#include <examples/biomineralization/material/components/substrate.hh>
#include <examples/biomineralization/material/components/suspendedbiomass.hh>

// we include brine-co2 and h2o-co2 binary coefficients
#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>

#include <dumux/material/fluidsystems/nullparametercache.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

#include <assert.h>

// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
namespace Dumux
{
namespace FluidSystems
{
// In the BioMinSimpleChemistryFluid fluid system, we define all functions needed to describe the solids accounted for in our simulation
template <class Scalar,
          class CO2Table,
          class H2OType = Components::TabulatedComponent<Components::H2O<Scalar>> >
class BioMinSimpleChemistryFluid
: public Base<Scalar, BioMinSimpleChemistryFluid<Scalar, CO2Table, H2OType> >
{
    using ThisType = BioMinSimpleChemistryFluid<Scalar, CO2Table, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:

    // We use convenient declarations that we derive from the property system.
    typedef Components::CO2<Scalar, CO2Table> CO2;
    using H2O = H2OType;
    // export the underlying brine fluid system for the liquid phase, as brine is used as a "pseudo component"
    using Brine = Dumux::FluidSystems::ICPComplexSalinityBrine<Scalar, H2OType>;
    using Na = Components::SodiumIon<Scalar>;
    using Cl = Components::ChlorideIon<Scalar>;
    using Ca = Components::CalciumIon<Scalar>;
    using Urea = Components::Urea<Scalar>;
    using O2 = Components::O2<Scalar>;
    using Biosub = Components::Substrate<Scalar>;
    using Biosusp = Components::SuspBiomass<Scalar>;

    using Brine_CO2 = BinaryCoeff::Brine_CO2<Scalar, CO2Table, true>;

    // the type of parameter cache objects. this fluid system does not
    // cache anything, so it uses Dumux::NullParameterCache
    typedef Dumux::NullParameterCache ParameterCache;

    // We define phase-related indices and properties
    static const int numPhases = 2; // liquid and gas phases
    static const int wPhaseIdx = 0; // index of the liquid phase
    static const int nPhaseIdx = 1; // index of the gas phase
    static const int phase0Idx = wPhaseIdx;
    static const int phase1Idx = nPhaseIdx;

    // the phase names
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            "w",
            "n"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    // We check whether a phase is liquid
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != nPhaseIdx;
    }

    // We return whether a phase is gaseous
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == nPhaseIdx;
    }

    // We assume ideal mixtures
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    // We assume compressible fluid phases
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    // We define component-related indices and properties
    static const int numComponents = 9; // H2O, TotalC, Na, Cl, Ca,...
    static const int numSecComponents = 6;

    static const int H2OIdx = 0;
    static const int BrineIdx = 0;
    static const int TCIdx = 1;
    static const int wCompIdx = BrineIdx;
    static const int nCompIdx = TCIdx;
    static const int comp0Idx = wCompIdx;
    static const int comp1Idx = nCompIdx;

    static const int NaIdx  = 2;
    static const int ClIdx  = 3;
    static const int CaIdx  = 4;
    static const int UreaIdx  = 5;
    static const int O2Idx  = 6;
    static const int BiosubIdx  = 7;
    static const int BiosuspIdx  = 8;

    // The component names
    static std::string componentName(int compIdx)
    {

        switch (compIdx) {
        case BrineIdx: return H2O::name();
        case TCIdx: return "TotalC";
        case CaIdx: return Ca::name();
        case NaIdx: return Na::name();
        case ClIdx: return Cl::name();
        case BiosuspIdx: return Biosusp::name();
        case BiosubIdx: return Biosub::name();
        case O2Idx: return O2::name();
        case UreaIdx: return Urea::name();
        break;
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx); break;
        };
    }

    // The component molar masses
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();break;
// actually, the molar mass of brine is only needed for diffusion
// but since chloride and sodium are accounted for seperately
// only the molar mass of water is returned.
        case TCIdx: return CO2::molarMass();break;
        case CaIdx: return Ca::molarMass();break;
        case NaIdx: return Na::molarMass();break;
        case ClIdx: return Cl::molarMass();break;
        case BiosuspIdx: return Biosusp::molarMass();break;
        case BiosubIdx: return Biosub::molarMass();break;
        case O2Idx: return O2::molarMass();break;
        case UreaIdx: return Urea::molarMass();break;
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);break;
        };
    }

// The brine adapter to be able to use the brine fluidsystem which only considers a single salt, while having Na, Cl, and Ca contributing to salinity in our system.
private:
    struct BrineAdapterPolicy
    {
        using FluidSystem = Brine;

        static constexpr int phaseIdx(int brinePhaseIdx) { return wPhaseIdx; }
        static constexpr int compIdx(int brineCompIdx)
        {
            switch (brineCompIdx)
            {
                assert(brineCompIdx == Brine::H2OIdx
                        || brineCompIdx == Brine::NaIdx
                        || brineCompIdx == Brine::ClIdx
                        || brineCompIdx == Brine::CaIdx
                      );
                case Brine::H2OIdx: return H2OIdx;
                case Brine::NaIdx: return NaIdx;
                case Brine::ClIdx: return ClIdx;
                case Brine::CaIdx: return CaIdx;
                default: return 0; // this will never be reached, only needed to suppress compiler warning
            }
        }
    };

    template<class FluidState>
    using BrineAdapter = FluidStateAdapter<FluidState, BrineAdapterPolicy>;

public:


    // Initialize the fluid system
    static void init()
    {
        init(/*startTemp=*/275.15, /*endTemp=*/455.15, /*tempSteps=*/30,
             /*startPressure=*/1e4, /*endPressure=*/99e6, /*pressureSteps=*/500);
    }
    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
                     Scalar startPressure, Scalar endPressure, int pressureSteps)
    {
        std::cout << "Initializing tables for the pure-water properties.\n";
        H2O::init(startTemp, endTemp, tempSteps,
                            startPressure, endPressure, pressureSteps);
     }

    // The phase densities
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            Scalar pl = fluidState.pressure(phaseIdx);
            Scalar xlCO2 = fluidState.moleFraction(wPhaseIdx, TCIdx);
            Scalar xlH2O = fluidState.moleFraction(wPhaseIdx, H2OIdx);
            if(T < 273.15) {
                DUNE_THROW(NumericalProblem,
                        "Liquid density for Brine and CO2 is only "
                        "defined above 273.15K (is" << T << ")");
            }
            if(pl >= 2.5e8) {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined below 250MPa (is" << pl << ")");
            }

            // calling the brine adapter for brine density
            Scalar rho_brine = Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx)
                       *(H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                         + Na::molarMass()*fluidState.moleFraction(wPhaseIdx, NaIdx)
                         + Cl::molarMass()*fluidState.moleFraction(wPhaseIdx, ClIdx)
                         + Ca::molarMass()*fluidState.moleFraction(wPhaseIdx, CaIdx)
                        );
            Scalar rho_pure = H2O::liquidDensity(T, pl);
            Scalar rho_lCO2 = liquidDensityWaterCO2_(T, pl, xlH2O, xlCO2);
            Scalar contribCO2 = rho_lCO2 - rho_pure;
            return rho_brine + contribCO2;
        }

        else if (phaseIdx == nPhaseIdx)
            return H2O::gasDensity(T, fluidState.partialPressure(nPhaseIdx, H2OIdx))
                   + CO2::gasDensity(T, fluidState.partialPressure(nPhaseIdx, TCIdx));
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase molar densities
    using Base::molarDensity;
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx)
            return Brine::molarDensity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == nPhaseIdx)
            return H2O::gasMolarDensity(fluidState.temperature(phaseIdx),
                                        fluidState.partialPressure(phaseIdx, H2OIdx))
                   + CO2::gasMolarDensity(fluidState.temperature(phaseIdx),
                                          fluidState.partialPressure(phaseIdx, TCIdx));
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase viscosities
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wPhaseIdx)
            return Brine::viscosity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx == nPhaseIdx)
        {
            return CO2::gasViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);

    }

    // The component fugacity coefficients
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == nPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (pressure<0)
        {
            typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
            ComponentVector moleFractionsw;
            ComponentVector massFractionsw;

            for (int compIdx = 0; compIdx<numComponents;++compIdx)
            {
                moleFractionsw[compIdx] = fluidState.moleFraction(wPhaseIdx,compIdx);
                massFractionsw[compIdx] = fluidState.massFraction(wPhaseIdx,compIdx);
            }
        }

        assert(temperature > 0);
        assert(pressure > 0);

        // calulate the equilibrium composition for the given
        // temperature and pressure. TODO: calculateMoleFractions()
        // could use some cleanup.
        Scalar xgH2O, xlH2O;
        Scalar xlCO2, xgCO2;
        Scalar Xl_Sal =    fluidState.massFraction(wPhaseIdx, NaIdx)                    //Salinity= XNa+XCl+XCa
                        + fluidState.massFraction(wPhaseIdx, ClIdx)
                        + fluidState.massFraction(wPhaseIdx, CaIdx);
        Brine_CO2::calculateMoleFractions(temperature,
                                          pressure,
                                          Xl_Sal,
                                          /*knownPhaseIdx=*/ -1,
                                          xlCO2,
                                          xgH2O);

        xlCO2 = std::max(0.0, std::min(1.0, xlCO2));
        xgH2O = std::max(0.0, std::min(1.0, xgH2O));
        xlH2O = 1.0 - xlCO2;
        xgCO2 = 1.0 - xgH2O;

        // Only H2O, CO2, and O2 in the gas phase
        if (compIdx == BrineIdx) {
            Scalar phigH2O = 1.0;
            return phigH2O * xgH2O / xlH2O;
        }
        else if (compIdx == TCIdx)
        {
            Scalar phigCO2 = 1.0;
            return phigCO2 * xgCO2 / xlCO2;
        }
        else if (compIdx == O2Idx)
        {
            return Dumux::BinaryCoeff::H2O_O2::henry(temperature)/pressure;
        }
        // All other components are assumed not be present in the gas phase
        else
        {
            return 1e-20;
        }
    }

    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    };

    // The binary diffusion coefficients of the components and the phases main component
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);

       Scalar temperature = fluidState.temperature(phaseIdx);
       Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx) {
            assert(compIIdx == H2OIdx);
            Scalar result = 0.0;
            if(compJIdx == TCIdx)
                 result = Brine_CO2::liquidDiffCoeff(temperature, pressure);
            else if(compJIdx == O2Idx)     //TODO Test O2 diffusion coeff from binarycoeff.
                result = Dumux::BinaryCoeff::H2O_O2::liquidDiffCoeff(temperature, pressure);
            else //all other components
                result = 1.587e-9;  //[m²/s]    //J. Phys. D: Appl. Phys. 40 (2007) 2769-2776 //old Value from Anozie 1e-9

            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            assert(phaseIdx == nPhaseIdx);
            assert(compIIdx == TCIdx);
            Scalar result = 0.0;
            if(compJIdx == H2OIdx)     //TODO Test H2O diffusion coeff from binarycoeff.
             result = Brine_CO2::gasDiffCoeff(temperature, pressure);
//                 result = 1e-5;
            else if(compJIdx == O2Idx)     //TODO Test O2 diffusion coeff from binarycoeff.
                result = Dumux::BinaryCoeff::H2O_O2::gasDiffCoeff(temperature, pressure);
            else //all other components
                result = 0.0;

            Valgrind::CheckDefined(result);
            return result;
        }
    };

    // The phase enthalpies
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            Scalar XlCO2 = fluidState.massFraction(phaseIdx, TCIdx);
            Scalar XlSal = fluidState.massFraction(wPhaseIdx, NaIdx)                    //Salinity= XNa+XCl+XCa
                         + fluidState.massFraction(wPhaseIdx, ClIdx)
                         + fluidState.massFraction(wPhaseIdx, CaIdx);
            Scalar result = liquidEnthalpyBrineCO2_(temperature,
                                                    pressure,
                                                    XlSal,
                                                    XlCO2);
            Valgrind::CheckDefined(result);
            return result;
        }
        else {
            Scalar XCO2 = fluidState.massFraction(nPhaseIdx, TCIdx);
            Scalar XBrine = fluidState.massFraction(nPhaseIdx, H2OIdx);

            Scalar result = 0;
            result += XBrine * Brine::gasEnthalpy(temperature, pressure);
            result += XCO2 * CO2::gasEnthalpy(temperature, pressure);
            Valgrind::CheckDefined(result);
            return result;
        }
    };

    // The phase thermal conductivities
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx)
            return Brine::thermalConductivity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);
        else if (phaseIdx ==nPhaseIdx)
            return CO2::gasThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    // The phase heat capacities
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
            return Brine::heatCapacity(BrineAdapter<FluidState>(fluidState), Brine::liquidPhaseIdx);

        // We assume NaCl not to be present in the gas phase here
        else if (phaseIdx ==nPhaseIdx)
            return CO2::gasHeatCapacity(T, p)*fluidState.moleFraction(nPhaseIdx, TCIdx)
                   + H2O::gasHeatCapacity(T, p)*fluidState.moleFraction(nPhaseIdx, H2OIdx);

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

private:

    //
    static Scalar liquidDensityWaterCO2_(Scalar temperature,
                                         Scalar pl,
                                         Scalar xlH2O,
                                         Scalar xlCO2)
    {
        const Scalar M_CO2 = CO2::molarMass();
        const Scalar M_H2O = H2O::molarMass();

        const Scalar tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const Scalar rho_pure = H2O::liquidDensity(temperature, pl);
        xlH2O = 1.0 - xlCO2; // xlH2O is available, but in case of a pure gas phase
                             // the value of M_T for the virtual liquid phase can become very large
        const Scalar M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const Scalar V_phi =
            (37.51 +
             tempC*(-9.585e-2 +
                    tempC*(8.74e-4 -
                           tempC*5.044e-7))) / 1.0e6;
        return 1/ (xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }

    static Scalar liquidEnthalpyBrineCO2_(Scalar T,
                                          Scalar p,
                                          Scalar S,
                                          Scalar X_CO2_w)
    {
        /* X_CO2_w : mass fraction of CO2 in brine */

        /* same function as enthalpy_brine, only extended by CO2 content */

        /*Numerical coefficents from PALLISER*/
        static const Scalar f[] = {
            2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10
        };

        /*Numerical coefficents from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = {
            { 9633.6, -4080.0, +286.49 },
            { +166.58, +68.577, -4.6856 },
            { -0.90963, -0.36524, +0.249667E-1 },
            { +0.17965E-2, +0.71924E-3, -0.4900E-4 }
        };

        Scalar theta, h_NaCl;
        Scalar m, h_ls, h_ls1, d_h;
        Scalar S_lSAT, delta_h;
        int i, j;
        Scalar delta_hCO2, hg, hw;

        theta = T - 273.15;

        S_lSAT = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;
        /*Regularization*/
        if (S>S_lSAT) {
            S = S_lSAT;
        }

        hw = H2O::liquidEnthalpy(T, p) /1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        /*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                        +((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* kJ/kg */

        m = (1E3/58.44)*(S/(1-S));
        i = 0;
        j = 0;
        d_h = 0;

        for (i = 0; i<=3; i++) {
            for (j=0; j<=2; j++) {
                d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
            }
        }
        /* heat of dissolution for halite according to Michaelides 1971 */
        delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without CO2 */
        h_ls1 =(1-S)*hw + S*h_NaCl + S*delta_h; /* kJ/kg */

        /* heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
           In the relevant temperature ranges CO2 dissolution is
           exothermal */
        delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

        /* enthalpy contribution of CO2 (kJ/kg) */
        hg = CO2::liquidEnthalpy(T, p)/1E3 + delta_hCO2;

        /* Enthalpy of brine with dissolved CO2 */
        h_ls = (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1E3; /*J/kg*/

        return (h_ls);
    };
};

} // end namespace
} // end namespace

#endif