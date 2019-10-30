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

#ifndef DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH
#define DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/grid/yaspgrid.hh>
// #include <dune/grid/uggrid.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <examples/biomineralization/material/fluidsystems/biominsimplechemistry.hh>
#include <examples/biomineralization/material/solidsystems/biominsolids.hh>
#include <examples/biomineralization/material/components/ammonia.hh>

#include <dumux/porousmediumflow/problem.hh>
// #include <dumux/porousmediumflow/2picp/model.hh>
    //No fancy secondary components in the simple chemistry version --> 2pncmin is sufficient!
#include <dumux/porousmediumflow/2pncmin/model.hh>

#include <dumux/material/binarycoefficients/brine_co2.hh>

#include "spatialparams.hh"
#include "material/co2tableslaboratory.hh"

#include "dumux/linear/seqsolverbackend.hh"

namespace Dumux
{
template <class TypeTag>
class MICPColumnProblemSimpleChemistry;

namespace Properties
{
NEW_PROP_TAG(CO2Tables); //!< The CO2 Tables that are used

// Create new type tags
namespace TTag {
struct MICPColumnSimpleChemistryTypeTag { using InheritsFrom = std::tuple<TwoPNCMin>; };
struct MICPColumnSimpleChemistryBoxTypeTag { using InheritsFrom = std::tuple<MICPColumnSimpleChemistryTypeTag, BoxModel>; };
} // end namespace TTag

// Set the grid type
// template<class TypeTag>
// struct Grid<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = Dune::UGGrid<2>; };
template<class TypeTag>
// struct Grid<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = Dune::YaspGrid<2>; };
struct Grid<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = Dune::YaspGrid<1>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = MICPColumnProblemSimpleChemistry<TypeTag>; };

//Set the CO2 tables used.
SET_TYPE_PROP(MICPColumnSimpleChemistryTypeTag, CO2Tables, Dumux::ICP::CO2Tables);

// set the fluidSystem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::BioMinSimpleChemistryFluid<Scalar, CO2Tables, H2OTabulated>;
};

// set the solidSystem
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::BioMinSolidPhase<Scalar>;
};


// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag> { using type = ICPSpatialParams<TypeTag>; };

template<class TypeTag>
struct Formulation<TypeTag, TTag::MICPColumnSimpleChemistryTypeTag>
{ static constexpr auto value = TwoPFormulation::p0s1; };

}

/*!
 * \ingroup TwoPNCSecCompMinModel
 * \ingroup ImplicitTestProblems
 * \brief Problem for enzyme-induced calcium carbonate precipitation
 *  */
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
    enum {
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
        xwNaIdx = FluidSystem::NaIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwUreaIdx = FluidSystem::UreaIdx,
        xwO2Idx = FluidSystem::O2Idx,
        xwBiosubIdx = FluidSystem::BiosubIdx,
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
        BiosubIdx = FluidSystem::BiosubIdx,
        BiosuspIdx = FluidSystem::BiosuspIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,

        //Index of the primary component of G and L phase
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
    MICPColumnProblemSimpleChemistry(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        Dune::FMatrixPrecision<>::set_singular_limit(1e-35);

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


        numInjections_ = getParam<int>("Injection.numInjections");
        injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(gridGeometry->gridView().size(codim));
//         scvPotGradNorm_.resize(gridGeometry->gridView().size(codim));

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

        // print file to make sure it is the right file
        std::cout << "Read file: " << injectionParameters_ << " ..." << std::endl << std::endl;
        while(!injectionData.eof())
        {
            getline(injectionData, row);
            std::cout << row << std::endl;
        }
        injectionData.close();

        //read data from file
        injectionData.open(injectionParameters_);

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

//      check the injection data against the number of injections specified in the parameter file
        if (injType_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection types specified in the injection data file do not match!"
                    <<"\n numInjections from parameter file = "<<numInjections_
                    <<"\n numInjTypes from injection data file = "<<injType_.size()
                    <<"\n Abort!\n";
            exit(1) ;
        }

        FluidSystem::init(/*startTemp=*/initTemperature_ -5.0, /*endTemp=*/initTemperature_ +5.0, /*tempSteps=*/5,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
    }

    void setTime( Scalar time )
    {
        time_ = time;
    }

    void setTimeStepSize( Scalar timeStepSize )
    {
        timeStepSize_ = timeStepSize;
    }

    void setEpisodeIdx( Scalar epiIdx )
    {
        episodeIdx_ = epiIdx;
    }

    int injectionType(int episodeIdx)
    {
        return injType_[episodeIdx];
    }

   /*!
    * \name Problem parameters
    */


   /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
//    const char *name() const
    const std::string name() const
    { return name_; }

    /* \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 25 degrees Celsius.
    */
    Scalar temperature() const
    {
        return initTemperature_; //
    };

    // \}

   /*!
    * \name Boundary conditions
    */
    // \{

    /*!
    * \brief Specifies which kind of boundary condition should be
    *        used for which equation on a given boundary segment.
    */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        Scalar zmax = this->gridGeometry().bBoxMax()[dim - 1];
        bcTypes.setAllNeumann();
        if(globalPos[dim - 1] > zmax - eps_)
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

   /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        boundary segment.
    */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);

        priVars = initial_(globalPos);
        priVars[xwBiosubIdx] = 0.0;
        priVars[xwO2Idx] = 0.0;
        return priVars;
//         return initial_(globalPos);
    }

   /*!
    * \brief Evaluate the initial value for a control volume.
    *
    * \param globalPos The global position
    *
    * For this method, the \a values parameter stores primary
    * variables.
    */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return initial_(globalPos);
    }

   /*!
    * \brief Evaluate the boundary conditions for a Neumann
    *        boundary segment.
    *
    * For this method, the \a values parameter stores the mass flux
    * in normal direction of each component. Negative values mean
    * influx.
    *
    * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
    */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
//         Scalar diameter = this->gridGeometry().bBoxMax()[0];
        Scalar diameter = 0.0254;
        Scalar waterFlux = injQ_/(3.14*diameter*diameter/4.); //[m/s]

        int injProcess = injType_[episodeIdx_];
//         std::cout<<"injProcess = " <<injProcess <<", injType_[" <<episodeIdx_ <<"]=" << injType_[episodeIdx_]<<std::endl;

//         negative values for injection
        if(globalPos[dim - 1]<= eps_)
        {
            //basic rinse injection (injProcess == -1 )
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

            if (injProcess == -1)   // rinse, used as standard injection fluid
            {
                //          do not change anything.
            }

            else if (injProcess == -99 || injProcess == 9) // no injection
            {
                values = 0.0; //mol/m²/s
            }

            else if (injProcess == 1 || injProcess == 11)  //ca-rich injection: ca and urea injected additionally to rinse-fluid, Na (pH) and Cl are also different(CaCl2)
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

            else if (injProcess == 0 || injProcess == 3 )   //urea-injections: urea is injected additionally to rinse-fluid
            {
                values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
            }

            else if(injProcess == 2)        //inoculation: same as rinse, but with bacteria
            {
                values[conti0EqIdx + xwBiosuspIdx] = -waterFlux * injBiosusp_ /FluidSystem::molarMass(xwBiosuspIdx);
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << injProcess);
            }
        }
        else
        {
               values = 0.0; //mol/m²/s
        }
        return values;
    }


   /*!
    * \name Volume terms
    */
    // \{

   /*!
    * \brief Evaluate the source term for all phases within a given
    *        sub-control-volume.
    *
    * This is the method for the case where the source term is
    * potentially solution dependent and requires some quantities that
    * are specific to the fully-implicit method.
    *
    * \param values The source and sink values for the conservation equations in units of
    *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
    * \param element The finite element
    * \param fvGeometry The finite-volume geometry
    * \param elemVolVars All volume variables for the element
    * \param scv The subcontrolvolume
    *
    * For this method, the \a values parameter stores the conserved quantity rate
    * generated or annihilate per volume unit. Positive values mean
    * that the conserved quantity is created, negative ones mean that it vanishes.
    * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
    */
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
        Scalar vmax = kurease_; //According to new experimental results without pH dependence
        Scalar Zub = kub_ *  massBiofilm; //new nub_=1 ! pow(massBiofilm,nub_);        // [kgurease/m³]
        Scalar rurea = vmax * Zub * mUrea / (Ku_ + mUrea); //[mol/m³s] //no NH4 in simple chemistry!
//         Scalar rurea = vmax * Zub * mUrea / ((Ku_ + mUrea) * (1 + mNH4 / KNH4_)); //[mol/m³s]

        // compute precipitation rate of calcite, no dissolution! Simplification: rprec = rurea
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
        dcf += rprec * SolidSystem::molarMass(SolidSystem::CalciteIdx) /               //[1/s]
                (densityCalcite * (initialPorosity - volFracCalcite));
        Scalar dcb = dc0_; //no pH for simple chemistry!!       //[1/s]
//         Scalar dcb = dc0_ * (1 + (mH * mH)/KpHa_ + KpHb_/mH); //KpHb_ = 0!!!!!!        //[1/s]
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

   /*!
    * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
    */

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }

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
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
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
    /*!
        * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
        *
        * \param XwNaCl the XwNaCl [kg NaCl / kg solution]
        */
    static Scalar massTomoleFrac_(Scalar XwNaCl)
    {
        const Scalar Mw = FluidSystem::molarMass(wCompIdx);  // 18.015e-3; /* molecular weight of water [kg/mol] */
        const Scalar Ms = FluidSystem::molarMass(NaIdx) + FluidSystem::molarMass(ClIdx); // 58.44e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = XwNaCl;
        /* XwNaCl: conversion from mass fraction to mol fraction */
        const Scalar xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return xwNaCl;
    }

    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(FluidSystem::H2OIdx);
        return molalityX;
    }


    static constexpr Scalar eps_ = 1e-6;

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

    // biomass parameters
    Scalar ca1_;
    Scalar ca2_;
    Scalar cd1_;
    Scalar dc0_;
    Scalar kmue_ ;
    Scalar F_;
    Scalar Ke_;
    Scalar Ks_;
    Scalar Yield_;
    // urease parameters
    Scalar kub_;
    Scalar kurease_;
    Scalar Ku_;

    Scalar injQ_;
    Scalar angle_;

    Scalar injTC_;//5.8e-7;             // [kg/kg]
    Scalar injNa_;//0.00379;                // [kg/m³]
    Scalar injCa_;//13.3593333333;          // [kg/m³]      //computed from 0.333mol/l CaCl2
    Scalar injUrea_;//20;                   // [kg/m³]
    Scalar injTNH_;//3.183840574;//3.184;               // [kg/m³]      //computed from 10 g/l NH4Cl
    Scalar injO2_;//0.008;                  // [kg/m³]
    Scalar injSub_;//3;                 // [kg/m³]
    Scalar injBiosusp_;//0.0675;            // [kg/m³]      //2.7e8 cfu/ml (40e8cfu/ml~1g/l)
    Scalar injNaCorr_;
    Scalar injTemperature_;
    Scalar injPressure_;

    int numInjections_;
    std::string injectionParameters_;
    std::vector<int> injType_;
    std::string name_;
    std::vector<Scalar> permeability_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
    std::shared_ptr<GridVariables> gridVariables_;
};
} //end namespace

#endif
