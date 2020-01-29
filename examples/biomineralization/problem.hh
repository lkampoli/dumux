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

// ## Header guard
// The header guard (or include guard) prevents compilation errors due to duplicate definitions. Here, a unique name needs to be defined for the header file:
#ifndef DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH
#define DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH

// ## Include files
// We use the dune yasp grid.
#include <dune/grid/yaspgrid.hh>

// We include the box discretization scheme.
#include <dumux/discretization/method.hh>
#include <dumux/discretization/box.hh>
// We include eval gradients to evaluate the pressure gradient for the detachment of biomass
#include <dumux/discretization/evalgradients.hh>

// We include the header which are needed for the biomineralization problem and model.
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>

// We include the necessary material files
#include <examples/biomineralization/material/fluidsystems/biominsimplechemistry.hh>
#include <examples/biomineralization/material/solidsystems/biominsolids.hh>
#include <dumux/material/components/ammonia.hh>
#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <examples/biomineralization/material/co2tableslaboratory.hh>

// We include the header that specifies all spatially variable parameters.
#include "spatialparams.hh"

#include "dumux/linear/seqsolverbackend.hh"

// ## Define basic properties for our simulation
// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don't clash with symbols from other libraries you may want to use in conjunction with Dumux. One could use these functions and classes by prefixing every use of these names by ::, but that would quickly become cumbersome and annoying. Rather, we simply import the entire Dumux namespace for general use.
namespace Dumux
{
// The problem class is forward declared.
template <class TypeTag>
class MICPColumnProblemSimpleChemistry;

// We enter the namespace Properties, which is a sub-namespace of the namespace Dumux.
namespace Properties
{
// We define the CO2 tables used as a new property to be used for the Fluidsystem
NEW_PROP_TAG(CO2Tables);

// We create new type tag for our simulation which inherits from the 2pncmin model and the box discretization
namespace TTag {
struct MICPColumnSimpleChemistryTypeTag { using InheritsFrom = std::tuple<TwoPNCMin>; };
struct MICPColumnSimpleChemistryBoxTypeTag { using InheritsFrom = std::tuple<MICPColumnSimpleChemistryTypeTag, BoxModel>; };
} // end namespace TTag

// We set the grid to a 1D Yasp Grid
template<class TypeTag>
struct Grid<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = Dune::YaspGrid<1>; };

// We set the problem  used for our simulation, defining boundary and initial conditions (see below)
template<class TypeTag>
struct Problem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = MICPColumnProblemSimpleChemistry<TypeTag>; };

// We set the CO2 tables used for our simulation
SET_TYPE_PROP(MICPColumnSimpleChemistryTypeTag, CO2Tables, Dumux::ICP::CO2Tables);

// We set the fluidSystem  used for our simulation
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::BioMinSimpleChemistryFluid<Scalar, CO2Tables, H2OTabulated>;
};

// We set the solidSystem  used for our simulation
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::BioMinSolidPhase<Scalar>;
};

// We define the spatial parameters for our simulation. The values are specified in the corresponding spatialparameters header file, which is included above.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = ICPSpatialParams<TypeTag>; };

// We set the two-phase primary variable formulation used for our simulation
template<class TypeTag>
struct Formulation<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{ static constexpr auto value = TwoPFormulation::p0s1; };

}// We leave the namespace Properties.

// ## The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a biomineralization problem in porous media, we inherit from the basic porous-media-flow problem
template <class TypeTag>
class MICPColumnProblemSimpleChemistry : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NH3 = Components::Ammonia<Scalar>; //no ammonia in Fluidsystem, but we need the molar mass for boundary conditions!
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    // We define some indices for convenience to be used later when defining the initial and boundary conditions
    enum {
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
        xwNaIdx = FluidSystem::NaIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwUreaIdx = FluidSystem::UreaIdx,
        xwO2Idx = FluidSystem::O2Idx,
        xwBiosubIdx = FluidSystem::GlucoseIdx,
        xwBiosuspIdx = FluidSystem::BiosuspIdx,
        phiBiofilmIdx = numComponents,
        phiCalciteIdx = numComponents + 1,

        //Indices of the components
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
        NaIdx = FluidSystem::NaIdx,
        ClIdx = FluidSystem::ClIdx,
        CaIdx = FluidSystem::CaIdx,
        UreaIdx = FluidSystem::UreaIdx,
        O2Idx = FluidSystem::O2Idx,
        BiosubIdx = FluidSystem::GlucoseIdx,
        BiosuspIdx = FluidSystem::BiosuspIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase State
        wPhaseOnly = Indices::firstPhaseOnly,
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
//     using Chemistry = GetPropType<TypeTag, Properties::Chemistry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

public:
    // This is the constructor of our problem class.
    MICPColumnProblemSimpleChemistry(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // We read the parameters from the params.input file.
        name_  = getParam<std::string>("Problem.Name");

        // biomass parameters
        ca1_ = getParam<Scalar>("BioCoefficients.ca1");
        ca2_ = getParam<Scalar>("BioCoefficients.ca2");
        cd1_ = getParam<Scalar>("BioCoefficients.cd1");
        dc0_ = getParam<Scalar>("BioCoefficients.dc0");
        kmue_  = getParam<Scalar>("BioCoefficients.kmue");
        F_ = getParam<Scalar>("BioCoefficients.F");
        Ke_ = getParam<Scalar>("BioCoefficients.Ke");
        Ks_ = getParam<Scalar>("BioCoefficients.Ks");
        Yield_ = getParam<Scalar>("BioCoefficients.Yield");

        //ureolysis kinetic parameters
        kub_ = getParam<Scalar>("UreolysisCoefficients.kub");
        kurease_ = getParam<Scalar>("UreolysisCoefficients.kurease");
        Ku_ = getParam<Scalar>("UreolysisCoefficients.Ku");


        //initial values
        densityW_ = getParam<Scalar>("Initial.initDensityW");
        initPressure_ = getParam<Scalar>("Initial.initPressure");

        initxwTC_ = getParam<Scalar>("Initial.initxwTC");
        initxwNa_ = getParam<Scalar>("Initial.initxwNa");
        initxwCl_ = getParam<Scalar>("Initial.initxwCl");
        initxwCa_ = getParam<Scalar>("Initial.initxwCa");
        initxwUrea_ = getParam<Scalar>("Initial.initxwUrea");
        initxwTNH_ = getParam<Scalar>("Initial.initxwTNH");
        initxwO2_ = getParam<Scalar>("Initial.initxwO2");
        initxwBiosub_ = getParam<Scalar>("Initial.initxwSubstrate");
        initxwBiosusp_ = getParam<Scalar>("Initial.initxwSuspendedBiomass");
        initCalcite_ = getParam<Scalar>("Initial.initCalcite");
        initBiofilm_ = getParam<Scalar>("Initial.initBiofilm");
        initTemperature_ = getParam<Scalar>("Initial.initTemperature");

        xwNaCorr_ = getParam<Scalar>("Initial.xwNaCorr");
        xwClCorr_ = getParam<Scalar>("Initial.xwClCorr");

        //injection values
        injQ_ = getParam<Scalar>("Injection.injVolumeflux");

        injTC_ = getParam<Scalar>("Injection.injTC");
        injNa_ = getParam<Scalar>("Injection.injNa");
        injCa_ = getParam<Scalar>("Injection.injCa");
        injUrea_ = getParam<Scalar>("Injection.injUrea");
        injTNH_ = getParam<Scalar>("Injection.injTNH");
        injO2_ = getParam<Scalar>("Injection.injO2");
        injSub_ = getParam<Scalar>("Injection.injSubstrate");
        injBiosusp_= getParam<Scalar>("Injection.injSuspendedBiomass");

        injNaCorr_ = getParam<Scalar>("Injection.injNaCorr");
        injTemperature_ = getParam<Scalar>("Injection.injTemperature");
        injPressure_ = getParam<Scalar>("Injection.injPressure");

        // We get the number of injections and the injection data file name from params.input
        numInjections_ = getParam<int>("Injection.numInjections");
        injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");

        // We resize the permeability vector contaning the permeabilities for the additional output
        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(gridGeometry->gridView().size(codim));

        // We read from the injection data file which injection type we have in each episode. We will use this in the Neumann boundary condition to set time dependend, changing boundary conditions. We do this similarly to the episode ends in the main file.
        std::ifstream injectionData;
        std::string row;
        injectionData.open( injectionParameters_); // open the Injection data file
        if (not injectionData.is_open())
        {
            std::cerr << "\n\t -> Could not open file '"
                    << injectionParameters_
                    << "'. <- \n\n\n\n";
            exit(1) ;
        }
        int tempType = 0;

        while(!injectionData.eof())
        {
            getline(injectionData, row);

            if(row == "InjectionTypes")
            {
                getline(injectionData, row);
                while(row != "#")
                {
                    if (row != "#")
                        {
                        std::istringstream ist(row);
                        ist >> tempType;
                        injType_.push_back(tempType);
                        }
                    getline(injectionData, row);
                }
            }
        }

        injectionData.close();

//   We   check the injection data against the number of injections specified in the parameter file and print an error message if the test fails
        if (injType_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection types specified in the injection data file do not match!"
                    <<"\n numInjections from parameter file = "<<numInjections_
                    <<"\n numInjTypes from injection data file = "<<injType_.size()
                    <<"\n Abort!\n";
            exit(1) ;
        }

        // We initialize the fluidsystem
        FluidSystem::init(/*startTemp=*/initTemperature_ -5.0, /*endTemp=*/initTemperature_ +5.0, /*tempSteps=*/5,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
    }

    // A function to set the time in the problem, this is done during the time loop in main.cc
    void setTime( Scalar time )
    {
        time_ = time;
    }

    // A function to set the time step size in the problem, this is done during the time loop in main.cc.
    // We need the time step size to regularize the reactive source terms.
    void setTimeStepSize( Scalar timeStepSize )
    {
        timeStepSize_ = timeStepSize;
    }

    // A function to set the episode index in the problem, this is done during the time loop in main.cc.
    // We need the episode index to choose the right Neumann boundary condition for each episode based on the injection data.
    void setEpisodeIdx( Scalar epiIdx )
    {
        episodeIdx_ = epiIdx;
    }

    // A function to return the injectionType
    int injectionType(int episodeIdx)
    {
        return injType_[episodeIdx];
    }

    // Get the problem name. It is used as a prefix for files generated by the simulation.
    const std::string name() const
    { return name_; }

    // Return the temperature
    Scalar temperature() const
    {
        return initTemperature_;
    };

    // We specify the boundary condition type.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        // We set all to Neumann, except for the top boundary, which is set to Dirichlet.
        Scalar zmax = this->gridGeometry().bBoxMax()[dim - 1];
        bcTypes.setAllNeumann();
        if(globalPos[dim - 1] > zmax - eps_)
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

    // We define the Dirichlet boundary conditions
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        // We recycle the initial conditions, but additionally enforce that substrate and oxygen necessary for biomass growth are zero.
        priVars = initial_(globalPos);
        priVars[xwBiosubIdx] = 0.0;
        priVars[xwO2Idx] = 0.0;
        return priVars;
    }


    // We define the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // We actually define the initial conditions in a private function...
        return initial_(globalPos);
    }


    // We define the Neumann boundary conditions, which are time or rather episode depended.
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        // We calculate the injected water velocity based on the volume flux injected into a 1 inch diameter column.
        Scalar diameter = 0.0254;
        Scalar waterFlux = injQ_/(3.14*diameter*diameter/4.); //[m/s]

        int injProcess = injType_[episodeIdx_];

        // We only have injection at the bottom, negative values for injection
        if(globalPos[dim - 1]<= eps_)
        {
            //We set the injection BC to the basic rinse injection and later only change the BC for those components that are different.
            values[conti0EqIdx + wCompIdx] = -waterFlux * 996/FluidSystem::molarMass(wCompIdx);
            values[conti0EqIdx + nCompIdx] = -waterFlux * injTC_*996 /FluidSystem::molarMass(nCompIdx);
            values[conti0EqIdx + xwCaIdx] = 0;
            values[conti0EqIdx + xwBiosuspIdx] = 0;
            values[conti0EqIdx + xwBiosubIdx] = -waterFlux * injSub_ /FluidSystem::molarMass(xwBiosubIdx);
            values[conti0EqIdx + xwO2Idx] = -waterFlux * injO2_ /FluidSystem::molarMass(O2Idx);
            values[conti0EqIdx + xwUreaIdx] = 0;
            values[conti0EqIdx + phiCalciteIdx] = 0;
            values[conti0EqIdx + phiBiofilmIdx] = 0;
            values[conti0EqIdx + xwNaIdx] = -waterFlux * (injNa_ + injNaCorr_) /FluidSystem::molarMass(NaIdx);
            values[conti0EqIdx + xwClIdx] = -waterFlux *injTNH_ /NH3::molarMass()     //NH4Cl --->  mol Cl = mol NH4
                                            -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);      //NaCl ---> mol Cl = mol Na

            //injProcess == -1 codes for the basic rinse injection to which we set all the BC above, so we do not change anything here.
            // Rinse does neither contain suspended biomass (Biosusp), nor calcium chloride, nor urea
            if (injProcess == -1)
            {
                //          do not change anything.
            }

            // We also have no-injection, no-flow periods, which are coded for by -99 or 9. Thus, for no flow, we set the Neumann BC to zero for all components.
            else if (injProcess == -99 || injProcess == 9)
            {
                values = 0.0; //mol/m²/s
            }

            // injProcess == 1 or injProcess == 11 codes for an injection of mineralization medium containing urea and calcium chloride.
            // Thus, we add BC terms for those components.
            // Additionally, we need to adjust the amount of water injected due to the high concentration of other components injected.
            // Finally, we need to adapt the injected NaCorr concentration to account fo the lower pH.
            else if (injProcess == 1 || injProcess == 11)
            {
                values[conti0EqIdx + wCompIdx] = - waterFlux * 0.8716 * densityW_ /FluidSystem::molarMass(wCompIdx);    //TODO 0.8716 check factor!!!
                values[conti0EqIdx + nCompIdx] = - waterFlux * injTC_ * densityW_ /FluidSystem::molarMass(nCompIdx);
                values[conti0EqIdx + xwCaIdx] = - waterFlux * injCa_/FluidSystem::molarMass(CaIdx);
                values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
                values[conti0EqIdx + xwNaIdx] = - waterFlux * injNa_ /FluidSystem::molarMass(NaIdx)
                                                - waterFlux * injNaCorr_ /FluidSystem::molarMass(NaIdx)* 0.032;
                values[conti0EqIdx + xwClIdx] = - waterFlux * injTNH_ /NH3::molarMass()               //NH4Cl --->  mol Cl = mol NH4
                                                - waterFlux * 2 * injCa_/FluidSystem::molarMass(CaIdx)             //+CaCl2 --->  mol Cl = mol Ca*2
                                                -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);      //NaCl ---> mol Cl = mol Na
            }

            // injProcess == 0 or injProcess == 3 codes for a resuscitation injection to regrow biomass.
            // It is similar to a rinse injection, but with added urea, which is what we add to the basic rinse
            else if (injProcess == 0 || injProcess == 3 )
            {
                values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
            }

            // injProcess == 2 codes for a inoculation or biomass injection.
            // It is similar to a rinse injection, but with added suspended biomass, which is what we add to the basic rinse
            else if(injProcess == 2)
            {
                values[conti0EqIdx + xwBiosuspIdx] = -waterFlux * injBiosusp_ /FluidSystem::molarMass(xwBiosuspIdx);
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << injProcess);
            }
        }
        // No flow conditions anywhere else than the bottom.
        else
        {
               values = 0.0; //mol/m²/s
        }
        return values;
    }

    // We calculate the reactive source and sink terms.
    // For the details, see the "simplified chemistry case" in the dissertation of Hommel available at https://elib.uni-stuttgart.de/handle/11682/8787
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        auto elemSol = elementSolution<Element, ElementVolumeVariables, FVElementGeometry>(element, elemVolVars, fvGeometry);
        const auto gradPw = evalGradients(element,
                                        element.geometry(),
                                        this->gridGeometry(),
                                        elemSol,
                                        scv.center())[pressureIdx];
        Scalar scvPotGradNorm = gradPw.two_norm();

        //define and compute some parameters for siplicity:
        Scalar porosity = elemVolVars[scv].porosity();
        Scalar initialPorosity = 1.0;
        for (int i=SolidSystem::numComponents-SolidSystem::numInertComponents; i<SolidSystem::numComponents ; ++i)
        {
            initialPorosity   -= elemVolVars[scv].solidVolumeFraction(i);
        }
        Scalar Sw  =  elemVolVars[scv].saturation(wPhaseIdx);
        Scalar xlSalinity = elemVolVars[scv].moleFraction(wPhaseIdx,NaIdx)
                            + elemVolVars[scv].moleFraction(wPhaseIdx,CaIdx)
                            + elemVolVars[scv].moleFraction(wPhaseIdx,ClIdx);
        Scalar densityBiofilm = elemVolVars[scv].solidComponentDensity(SolidSystem::BiofilmIdx);
        Scalar densityCalcite = elemVolVars[scv].solidComponentDensity(SolidSystem::CalciteIdx);
        Scalar cBio = elemVolVars[scv].moleFraction(wPhaseIdx, BiosuspIdx) * elemVolVars[scv].molarDensity(wPhaseIdx) * FluidSystem::molarMass(BiosuspIdx);      //[kg_suspended_Biomass/m³_waterphase]
        if(cBio < 0)
            cBio = 0;
        Scalar volFracCalcite = elemVolVars[scv].solidVolumeFraction(SolidSystem::CalciteIdx);
        if (volFracCalcite < 0)
            volFracCalcite = 0;
        Scalar volFracBiofilm = elemVolVars[scv].solidVolumeFraction(SolidSystem::BiofilmIdx);
        if (volFracBiofilm < 0)
            volFracBiofilm = 0;
        Scalar massBiofilm = densityBiofilm * volFracBiofilm;
        Scalar cSubstrate = elemVolVars[scv].moleFraction(wPhaseIdx, BiosubIdx) * elemVolVars[scv].molarDensity(wPhaseIdx) * FluidSystem::molarMass(BiosubIdx);  //[kg_substrate/m³_waterphase]
        if(cSubstrate < 0)
            cSubstrate = 0;
        Scalar cO2 = elemVolVars[scv].moleFraction(wPhaseIdx, O2Idx) * elemVolVars[scv].molarDensity(wPhaseIdx) * FluidSystem::molarMass(O2Idx);                 //[kg_oxygen/m³_waterphase]
        if(cO2 < 0)//1e-10)
            cO2 = 0;

        Scalar mUrea = moleFracToMolality(elemVolVars[scv].moleFraction(wPhaseIdx,UreaIdx), xlSalinity, elemVolVars[scv].moleFraction(wPhaseIdx,nCompIdx));  //[mol_urea/kg_H2O]
        if (mUrea < 0)
            mUrea = 0;

        // compute rate of ureolysis:
        Scalar vmax = kurease_;
        Scalar Zub = kub_ *  massBiofilm;   // [kgurease/m³]
        Scalar rurea = vmax * Zub * mUrea / (Ku_ + mUrea); //[mol/m³s]

        // compute precipitation rate of calcite, no dissolution! Simplification: rprec = rurea
        // additionally regularize the precipitation rate that we do not precipitate more calcium than available.
        Scalar rprec = rurea;
        if(rprec >
            elemVolVars[scv].moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * elemVolVars[scv].molarDensity(wPhaseIdx) / timeStepSize_)
        {
            rprec =  elemVolVars[scv].moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * elemVolVars[scv].molarDensity(wPhaseIdx) / timeStepSize_;
        }

        //compute biomass growth coefficient and rate
        Scalar mue = kmue_ * cSubstrate / (Ks_ + cSubstrate) * cO2 / (Ke_ + cO2);// [1/s]
        Scalar rgf = mue * massBiofilm;                //[kg/m³s]
        Scalar rgb = mue * porosity * Sw * cBio;   //[kg/m³s]

        // compute biomass decay coefficient and rate:
        Scalar dcf = dc0_;
        dcf += rprec * SolidSystem::molarMass(SolidSystem::CalciteIdx) /
                (densityCalcite * (initialPorosity - volFracCalcite));
        Scalar dcb = dc0_;       //[1/s]
        Scalar rdcf = dcf * massBiofilm; //[kg/m³s]
        Scalar rdcb = dcb * porosity * Sw * cBio;      //[kg/m³s]

        // compute attachment coefficient and rate:
        Scalar ka = ca1_ * volFracBiofilm + ca2_;          //[1/s]
        Scalar ra = ka * porosity * Sw * cBio;             //[kg/m³s]

        // compute detachment coefficient and rate:
        Scalar cd2 = volFracBiofilm / (initialPorosity - volFracCalcite);      //[-]
        Scalar kd = cd1_ * pow((porosity * Sw * scvPotGradNorm),0.58) + cd2 * mue;  //[1/s]
        Scalar rd = kd * massBiofilm;                      //[kg/m³s]

        // rprec[mol/m³s]
        // rurea[mol/m³s]
        // rgb + rdcb + ra + rd [kg/m³s]
        // source[kg/m³s]
        source[wCompIdx] += 0;
        source[nCompIdx] += rurea - rprec;
        source[NaIdx] += 0;
        source[ClIdx] += 0;
        source[CaIdx] += - rprec;
        source[UreaIdx] += - rurea;
        source[O2Idx] += -(rgf + rgb) *F_/Yield_ / FluidSystem::molarMass(O2Idx);
        source[BiosubIdx] += -(rgf + rgb) / Yield_ / FluidSystem::molarMass(BiosubIdx);
        source[BiosuspIdx] += (rgb - rdcb - ra + rd) / FluidSystem::molarMass(BiosuspIdx);
        source[phiBiofilmIdx] += (rgf - rdcf + ra - rd) / SolidSystem::molarMass(SolidSystem::BiofilmIdx);
        source[phiCalciteIdx] += + rprec;

        return source;
    }

    // Function to return the permeability for additional vtk output
    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }

    // Function to update the permeability for additional vtk output
    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = volVars.permeability();
            }
        }
    }

    void setGridVariables(std::shared_ptr<GridVariables> gv)
    { gridVariables_ = gv; }

private:
    // Internal method for the initial condition reused for the dirichlet conditions.
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        priVars[pressureIdx] = initPressure_ ; //70e5; // - (maxHeight - globalPos[1])*densityW_*9.81; //p_atm + rho*g*h
        priVars[switchIdx] = initxwTC_;
        priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
        priVars[xwCaIdx] = initxwCa_;
        priVars[xwUreaIdx] = initxwUrea_;
        priVars[xwO2Idx] = initxwO2_;
        priVars[xwBiosubIdx] = initxwBiosub_;
        priVars[xwBiosuspIdx] = initxwBiosusp_;
        priVars[phiBiofilmIdx] = initBiofilm_; // [m^3/m^3]
        priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]

        return priVars;
    }

    // Internal method to calculate the molality of a component based on its mole fraction.
    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(FluidSystem::H2OIdx);
        return molalityX;
    }


    // eps is used as a small value for the definition of the boundry conditions
    static constexpr Scalar eps_ = 1e-6;

    // initial condition parameters
    Scalar initPressure_;
    Scalar densityW_;//1087; // rhow=1087;
    Scalar initxwTC_;//2.3864e-7;       // [mol/mol]
    Scalar initxwNa_;//0;
    Scalar initxwCl_;//0;
    Scalar initxwCa_;//0;
    Scalar initxwUrea_;//0;
    Scalar initxwTNH_;//3.341641e-3;
    Scalar initxwO2_;//4.4686e-6;
    Scalar initxwBiosub_;//2.97638e-4;
    Scalar initxwBiosusp_;//0;
    Scalar xwNaCorr_;//2.9466e-6;
    Scalar xwClCorr_;//0;
    Scalar initBiofilm_;
    Scalar initCalcite_;
    Scalar initTemperature_;

    // biomass parameters for source/sink calculations
    Scalar ca1_;
    Scalar ca2_;
    Scalar cd1_;
    Scalar dc0_;
    Scalar kmue_ ;
    Scalar F_;
    Scalar Ke_;
    Scalar Ks_;
    Scalar Yield_;
    // urease parameters for source/sink calculations
    Scalar kub_;
    Scalar kurease_;
    Scalar Ku_;

    // injection parameters
    Scalar injQ_;
    Scalar injTC_;   // [kg/kg]
    Scalar injNa_;   // [kg/m³]
    Scalar injCa_;   // [kg/m³]      //computed from CaCl2
    Scalar injUrea_; // [kg/m³]
    Scalar injTNH_;  // [kg/m³]      //computed from NH4Cl
    Scalar injO2_;   // [kg/m³]
    Scalar injSub_;  // [kg/m³]
    Scalar injBiosusp_;  // [kg/m³]
    Scalar injNaCorr_;   // [kg/m³]
    Scalar injTemperature_;
    Scalar injPressure_;
    int numInjections_;
    std::string injectionParameters_;
    std::vector<int> injType_;

    // the problem name
    std::string name_;
    // the permeability for output
    std::vector<Scalar> permeability_;

    // timing parameters
    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
    std::shared_ptr<GridVariables> gridVariables_;
};
} //end namespace

#endif
