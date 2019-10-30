/*
 * bioCarbonicAcid.hh
 *
 *  Created on: 9.8.2011
 *      Author: hommel
 */

#ifndef BIO_CARBONIC_ACID_HH_
#define BIO_CARBONIC_ACID_HH_

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>
//#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/fluidsystems/biofluidsystem.hh>
#include <dumux/implicit/2pbiomin/properties.hh>
#include <dumux/material/components/h2o.hh>

#include <cmath>
#include <iostream>
#include <dumux/common/math.hh>

namespace Dumux
{
/*!
 * \brief The equilibrium chemistry is calculated in this class. The function calculateEquilbriumChemistry is used to
 * control the Newton Solver "newton1D". The chemical functions and derivations are implemented in the private part of
 * class.
 */
template <class TypeTag, class CO2Tables>
class BioCarbonicAcid
{
    typedef GetPropType<TypeTag, Properties::Scalar> Scalar;
    typedef GetPropType<TypeTag, Properties::FluidSystem> FluidSystem;
    typedef GetPropType<TypeTag, Properties::VolumeVariables> VolumeVariables;
    typedef GetPropType<TypeTag, Properties::PrimaryVariables> PrimaryVariables;
    typedef GetPropType<TypeTag, Properties::GridView> GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef GetPropType<TypeTag, Properties::FVElementGeometry> FVElementGeometry;

    typedef BioCarbonicAcid<TypeTag, CO2Tables> ThisType;
    typedef Dumux::H2O<Scalar> H2O;


public:

    BioCarbonicAcid()
{
        try
    {

        // biomass parameters
        ca1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, ca1);
        ca2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, ca2);
        cd1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, cd1);
        dc0_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, dc0);
        kmue_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, kmue);
        F_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, F);
        Ke_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, Ke);
        KpHa_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, KpHa);
        Ks_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, Ks);
        Yield_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, BioCoefficients, Yield);

        //ureolysis kinetic parameters
        kub_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, UreolysisCoefficients, kub);
        kurease_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, UreolysisCoefficients, kurease);
        KNH4_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, UreolysisCoefficients, KNH4);
        Ku_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, UreolysisCoefficients, Ku);

    }
    catch (Dumux::ParameterException &e) {
        std::cerr << e << ". Abort!\n";
        exit(1) ;
    }
}

    static const int wPhaseIdx = FluidSystem::wPhaseIdx;
    static const int nPhaseIdx = FluidSystem::nPhaseIdx;

    static const int wCompIdx = FluidSystem::wCompIdx;
    static const int nCompIdx = FluidSystem::nCompIdx;

    static const int H2OIdx = FluidSystem::H2OIdx;
    static const int CTotIdx = FluidSystem::TCIdx;
    static const int CaIdx = FluidSystem::CaIdx;
    static const int NaIdx = FluidSystem::NaIdx;
    static const int ClIdx = FluidSystem::ClIdx;
    static const int CalciteIdx = FluidSystem::CalciteIdx;

    static const int BiofilmIdx = FluidSystem::BiofilmIdx;
    static const int BiosuspIdx = FluidSystem::BiosuspIdx;
    static const int BiosubIdx = FluidSystem::BiosubIdx;
    static const int UreaIdx = FluidSystem::UreaIdx;
    static const int O2Idx = FluidSystem::O2Idx;

    static const int TNHIdx = FluidSystem::TNHIdx;

    static const int numComponents = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numMajorComponents;
    static const int numSecComponents = FluidSystem::numSecComponents;
    static const int numTotComponents = numComponents + numSecComponents;
    static const int numPhases = FluidSystem::numPhases;
    static const int numSPhases = FluidSystem::numSPhases;

    static const int cPhaseIdx = FluidSystem::cPhaseIdx;
    static const int bPhaseIdx = FluidSystem::bPhaseIdx;

    static const int phiBiofilmIdx = numComponents;
    static const int phiCalciteIdx = numComponents + 1;

    typedef Dune::FieldVector<Scalar, 4> Vector;   // Ionic Strength with NH4/totalnh
//    typedef Dune::FieldVector<Scalar, 3> Vector;     // Ionic Strength without NH4/totalnh
    typedef Dune::FieldVector<Scalar, 2> SolVector;
    typedef Dune::FieldVector<Scalar, numTotComponents> CompVector;

    typedef CompositionalSecCompFluidState<Scalar, FluidSystem> FluidState;


    //Return equlibrium constant for chemical equation:
    //H2CO3 <--> H + HCO3
    static Scalar const1(const Scalar pw, const Scalar T)
    {
        return(pow(10,-6.3));
    }

    //Return equlibrium constant for chemical equation:
    //HCO3 <--> H + CO3
    static Scalar const2(const Scalar pw, const Scalar T)
    {
        return(pow(10,-10.3));
    }

//    Return equlibrium constant for dissolution reaction:
//    CaCO3(s) <--> Ca + CO3
    static Scalar solubilityProductCaCO(const Scalar pw, const Scalar T)
    {
//        Scalar k2_fw = 6.5e-7; //molCaCO/m3/s
//        Scalar k2_bw = 1.9e2; //molCaCO/m3/s
//        Scalar Kgg = k2_fw/k2_bw;
//        return Kgg;
        return(4.8e-9);
    }

    //Return equlibrium constant for chemical equation:
    // H2O <--> H + OH
    static Scalar constW(const Scalar pw, const Scalar T)
    {
        return(1e-14);
    }
    //Return equlibrium constant for chemical equation:
    // NH4 <--> H + NH3
    /*static*/ Scalar consta(const Scalar pw, const Scalar T)
    {
        return(pow(10,-9.29)); //pow(10,-9.25)
    }


    static Scalar massFracToMolality(const Scalar massFracX, const Scalar molarMassX, const Scalar massFracSalinity,
            const Scalar massFracC)
    {
        Scalar molalityX = massFracX/molarMassX/(1- massFracSalinity - massFracC);
        return molalityX;
    }


    /*!
     * \brief Returns the mass fraction of a component x (kg x / kg solution) for a given
     * molality fraction (mol x / mol solution)
     * The salinity and the mole Fraction of CO2 are considered
     *
     */

    static Scalar molalityToMassFrac(Scalar molalityX, Scalar molarMassX, Scalar massFracSalinity, Scalar massFracCTot)
    {
        Scalar massFracX = molalityX * molarMassX * (1 - massFracSalinity - massFracCTot);
        return massFracX;
    }


    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

    static Scalar molalityToMoleFrac(Scalar molalityX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar moleFracX = molalityX * (1 - moleFracSalinity - moleFracCTot) * FluidSystem::molarMass(H2OIdx);
        return moleFracX;
    }


   void reactionSource(PrimaryVariables &q,
           const VolumeVariables &volVars,
           const Scalar absgradpw,
           const Scalar dt)
           {
//      //define and compute some parameters for siplicity:
     Scalar porosity = volVars.porosity();
     Scalar initialPorosity = volVars.initialPorosity();
     Scalar Sw  =  volVars.saturation(wPhaseIdx);
     Scalar xlSalinity = volVars.moleFracSalinity();
     Scalar densityBiofilm = volVars.density(bPhaseIdx);
     Scalar densityCalcite = volVars.density(cPhaseIdx);
     Scalar cBio = volVars.moleFraction(wPhaseIdx, BiosuspIdx) * volVars.molarDensity(wPhaseIdx) * FluidSystem::molarMass(BiosuspIdx);      //[kg_suspended_Biomass/m³_waterphase]
     if(cBio < 0)
         cBio = 0;
     Scalar volFracCalcite = volVars.precipitateVolumeFraction(cPhaseIdx);
     if (volFracCalcite < 0)
         volFracCalcite = 0;
     Scalar volFracBiofilm = volVars.precipitateVolumeFraction(bPhaseIdx);
     if (volFracBiofilm < 0)
         volFracBiofilm = 0;
     Scalar massBiofilm = densityBiofilm * volFracBiofilm;
     Scalar cSubstrate = volVars.moleFraction(wPhaseIdx, BiosubIdx) * volVars.molarDensity(wPhaseIdx) * FluidSystem::molarMass(BiosubIdx);  //[kg_substrate/m³_waterphase]
     if(cSubstrate < 0)
         cSubstrate = 0;
     Scalar cO2 = volVars.moleFraction(wPhaseIdx, O2Idx) * volVars.molarDensity(wPhaseIdx) * FluidSystem::molarMass(O2Idx);                 //[kg_oxygen/m³_waterphase]
     if(cO2 < 0)//1e-10)
         cO2 = 0;

     Scalar mUrea = moleFracToMolality(volVars.moleFraction(wPhaseIdx,UreaIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_urea/kg_H2O]
     if (mUrea < 0)
         mUrea = 0;

     // compute rate of urealysis:
     Scalar vmax = kurease_; //According to new experimental results without pH dependence
     Scalar Zub = kub_ *  massBiofilm; //new nub_=1 ! pow(massBiofilm,nub_);        // [kgurease/m³]
     Scalar rurea = vmax * Zub * mUrea / (Ku_ + mUrea);// * (1 + mNH4 / KNH4_)); //[mol/m³s]


     // compute dissolution and precipitation rate of calcite
     Scalar rdiss = 0;
     Scalar rprec = 0;

     if(volVars.moleFraction(wPhaseIdx,CaIdx)>1e-8)
         rprec = rurea;

     //compute biomass growth coefficient and rate
     Scalar mue = kmue_ * cSubstrate / (Ks_ + cSubstrate) * cO2 / (Ke_ + cO2);// [1/s]
     Scalar rgf = mue * massBiofilm;                //[kg/m³s]
     Scalar rgb = mue * porosity * Sw * cBio;   //[kg/m³s]

     if(cO2-(rgf+rgb)*dt<0)
     {
         mue =  mue * 0.99*cO2/((rgf+rgb)*dt);
         rgf = mue * massBiofilm;
         rgb = mue * porosity * Sw * cBio;
     }

     // compute biomass decay coefficient and rate:
     Scalar dcf = dc0_;
      dcf += rprec * FluidSystem::molarMass(CalciteIdx) /               //[1/s]
             (volVars.density(cPhaseIdx) * (initialPorosity - volFracCalcite));
//   Scalar dcb = dc0_ * (1 + (mH * mH)/KpHa_ + KpHb_/mH); //KpHb_ = 0!!!!!!        //[1/s]
     Scalar dcb = dc0_;     //[1/s]
     Scalar rdcf = dcf * massBiofilm; //[kg/m³s]
     Scalar rdcb = dcb * porosity * Sw * cBio;      //[kg/m³s]

     // compute attachment coefficient and rate:
     Scalar ka = ca1_ * volFracBiofilm + ca2_;          //[1/s]
     Scalar ra = ka * porosity * Sw * cBio;             //[kg/m³s]

     // compute detachment coefficient and rate:
     Scalar cd2 = volFracBiofilm / (initialPorosity - volFracCalcite);      //[-]
     Scalar kd = cd1_ * pow((porosity * Sw * absgradpw),0.58) + cd2 * mue;  //[1/s]
     Scalar rd = kd * massBiofilm;                      //[kg/m³s]

     // rdiss+rprec[mol/m³s]
     // rurea[mol/m³s]
     // rgb + rdcb + ra + rd [kg/m³s]
     // q[kg/m³s]
     q[wCompIdx] += 0;
     q[nCompIdx] += rurea - rprec + rdiss;
     q[NaIdx] += 0;
     q[ClIdx] += 0;
     q[CaIdx] += - rprec + rdiss;
     q[UreaIdx] += - rurea;
     q[TNHIdx] += 2 * rurea;
     q[O2Idx] += -(rgf + rgb) *F_/Yield_ / FluidSystem::molarMass(O2Idx);
     q[BiosubIdx] += -(rgf + rgb) / Yield_ / FluidSystem::molarMass(BiosubIdx);
     q[BiosuspIdx] += (rgb - rdcb - ra + rd) / FluidSystem::molarMass(BiosuspIdx);
     q[phiBiofilmIdx] += (rgf - rdcf + ra - rd) / FluidSystem::molarMass(BiofilmIdx);
     q[phiCalciteIdx] += + rprec - rdiss;

}

private:
    //Value of numerical derivative at xVar
    /*static*/ Scalar equationNumDeri(Scalar xVar)
    {
        Scalar eps = 1e-8;
        Scalar xRight = xVar + eps*xVar; // x + dx
        Scalar xLeft = xVar - eps*xVar; // x - dx
        Scalar fRight = equationValue(xRight); // f(x+dx)
        Scalar fLeft = equationValue(xLeft); // f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/xVar; // {f(x+dx) - f(x-dx)}/2dx
        return df;
     }



    Scalar absolute(Scalar x)
    {
        if(x<0.0)
        {
            return x*(-1);
        }
        else return x;
    }

    Scalar sign(Scalar x)
    {
        if(x > 0.0)
        {
           return 1;
        }
        else if (x < 0.0)
        {
           return -1;
        }
        else
        {
            return 0.0;
        }
    }



    int iter_; //Number of iterations the Newton solver needs until convergence
    Scalar pressure_;
    Scalar temperature_;
    Scalar salinity_;
    Scalar h2o_;
    Scalar co2_;
    Scalar hco3_;
    Scalar co3_;
    Scalar oh_;
    Scalar h_;
    Scalar ca_;
    Scalar na_;
    Scalar cl_;
    Scalar totalnh_;
    Scalar nh4_;
    Scalar initH_;
    Scalar ionicStrength_;
    Scalar cTot_;
    Scalar gammaH_;
    Scalar gammaCO2_;
    Scalar gammaCa_;
    Scalar gammaOH_;
    Scalar gammaHCO3_;
    Scalar gammaCO3_;
    Scalar gammaNH3_;
    Scalar gammaNH4_;
    SolVector fdf_; //Solution vector for the newtons solver every equation f solved by the newton solver for an unknown x
    // has to store f(x) in fdf_[0] and df/dx in fdf[1]
    Vector molality_;
    Vector charge_;
    Scalar x_;
    Scalar y_;
    Scalar k1_;
    Scalar k2_;
    Scalar kw_;
    Scalar ka_;
    Scalar apparentk1_;
    Scalar apparentk2_;
    Scalar apparentka_;
    bool newtonOrBisection_;

    static constexpr Scalar KpHb_ = 0;//9.14e-8;//[mol/kgH2O] Kim et al. 2000 //Not implemented by Anozie!!

    // biomass parameters
        Scalar ca1_;
        Scalar ca2_;
        Scalar cd1_;
        Scalar dc0_;
        Scalar kmue_ ;
        Scalar F_;
        Scalar Ke_;
        Scalar KpHa_;
        Scalar Ks_;
        Scalar Yield_;

    // urease parameters
        Scalar kub_;
        Scalar kurease_;

        Scalar KNH4_;
        Scalar Ku_;

public:

    // biomass parameters
        Scalar ca1()    {       return ca1_; }
        Scalar ca2()    {       return ca2_; }
        Scalar cd1()    {       return cd1_; }
        Scalar dc0()    {       return dc0_; }
        Scalar kmue()    {      return kmue_; }
        Scalar F()    {         return F_; }
        Scalar Ke()    {        return Ke_; }
        Scalar KpHa()    {      return KpHa_; }
        Scalar Ks()    {        return Ks_; }
        Scalar Yield()    {     return Yield_; }

    // calcite parameters
        Scalar ac()    {        return ac_; }
        Scalar kdiss1()    {    return kdiss1_; }
        Scalar kdiss2()    {    return kdiss2_; }
        Scalar kprec()    {     return kprec_; }
        Scalar ndiss()    {     return ndiss_; }
        Scalar nprec()    {     return nprec_; }
        Scalar Asw0()    {      return Asw0_; }

    // urease parameters
        Scalar kub()    {       return kub_; }
        Scalar kurease()    {   return kurease_; }
        Scalar KNH4()    {      return KNH4_; }
        Scalar Ku()    {        return Ku_; }

public:
    Scalar kprec() const
    {   return kprec_;}
    Scalar kub() const
    {   return kub_;}
    Scalar kurease() const
    {   return kurease_;}
    Scalar nprec() const
    {   return nprec_;}
    Scalar Asw0() const
    {   return Asw0_;}
//     Scalar Keu1() const
//     {   return Keu1_;}
//     Scalar Keu2() const
//     {   return Keu2_;}
    Scalar KNH4() const
    {   return KNH4_;}
    Scalar Ku() const
    {   return Ku_;}

    /*!
     * \brief Returns the mole fraction of NaCl \f$\mathrm{[mol \ NaCl / mol \ solution]}\f$  for a given mole fraction
     *
     * \param salinity the salinity \f$\mathrm{[kg \ NaCl / kg \ solution]}\f$
     */
    static Scalar salinityToMolFrac_(Scalar salinity) {

        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = salinity;
        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return x_NaCl;
    }
};

} // end namespace

#endif