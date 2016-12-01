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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the three-phase, two-component model.
 */
#ifndef DUMUX_3P2CNI_VOLUME_VARIABLES_HH
#define DUMUX_3P2CNI_VOLUME_VARIABLES_HH

#include <dumux/implicit/model.hh>
#include <dumux/common/math.hh>

#include <dune/common/parallel/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class ThreePWaterOilVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // constraint solvers
    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        pressureIdx = Indices::pressureIdx
    };

    // present phases
    enum {
        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly
    };

    typedef typename GridView::template Codim<0>::Entity Element;

    static const Scalar R; // universial gas constant

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;


    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        bool useSimpleModel = GET_PROP_VALUE(TypeTag, UseSimpleModel);

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        int globalIdx = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);

        int phasePresence = problem.model().phasePresence(globalIdx, isOldSol);

        if(!useSimpleModel)
        {
            /* first the saturations */
            if (phasePresence == threePhases)
            {
                sw_ = priVars[switch1Idx];
                sn_ = priVars[switch2Idx];
                sg_ = 1. - sw_ - sn_;
            }
            else if (phasePresence == wPhaseOnly)
            {
                sw_ = 1.;
                sn_ = 0.;
                sg_ = 0.;
            }
            else if (phasePresence == gnPhaseOnly)
            {
                sw_ = 0.;
                sn_ = priVars[switch1Idx];
                sg_ = 1. - sn_;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                sn_ = priVars[switch2Idx];
                sw_ = 1. - sn_;
                sg_ = 0.;
            }
            else if (phasePresence == gPhaseOnly)
            {
                sw_ = 0.;
                sn_ = 0.;
                sg_ = 1.;
            }
            else if (phasePresence == wgPhaseOnly)
            {
                sw_ = priVars[switch1Idx];
                sn_ = 0.;
                sg_ = 1. - sw_;
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(sg_);

            fluidState_.setSaturation(wPhaseIdx, sw_);
            fluidState_.setSaturation(gPhaseIdx, sg_);
            fluidState_.setSaturation(nPhaseIdx, sn_);

            /* now the pressures */
            if (phasePresence == threePhases || phasePresence == gnPhaseOnly || phasePresence == gPhaseOnly || phasePresence == wgPhaseOnly)
            {
                 pg_ = priVars[pressureIdx];

                 // calculate capillary pressures
                 Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
                 Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
                 Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

                 Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
                 Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

                 pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
                 pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;
            }
            else if (phasePresence == wPhaseOnly || phasePresence == wnPhaseOnly)
            {
                 pw_ = priVars[pressureIdx];

                 // calculate capillary pressures
                 Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
                 Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
                 Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

                 Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
                 Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

                 pn_ = pw_ + pcAlpha * pcnw + (1.-pcAlpha)*pcNW1;
                 pg_ = pn_ + pcAlpha * pcgn + (1.-pcAlpha)*(pcgw - pcNW1);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(pw_);

            fluidState_.setPressure(wPhaseIdx, pw_);
            fluidState_.setPressure(gPhaseIdx, pg_);
            fluidState_.setPressure(nPhaseIdx, pn_);

            /* now the temperature */
            if (phasePresence == wPhaseOnly || phasePresence == wnPhaseOnly || phasePresence == gPhaseOnly)
            {
                 temp_ = priVars[switch1Idx];
            }
            else if (phasePresence == threePhases)
            {
                 // temp from inverse pwsat and pnsat which have to sum up to pg
                 Scalar temp = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx); // initial guess
                 fluidState_.setTemperature(temp);
                 Scalar defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                     - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                 while(std::abs(defect) > 0.01) // simply a small number chosen ...
                 {
                     Scalar deltaT = 1.e-8 * temp;
                     fluidState_.setTemperature(temp+deltaT);
                     Scalar fUp = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                     fluidState_.setTemperature(temp-deltaT);
                     Scalar fDown = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                     temp = temp - defect * 2. * deltaT / (fUp - fDown);

                     fluidState_.setTemperature(temp);
                     defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                  - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);
                 }
                 temp_ = temp;
            }
            else if (phasePresence == wgPhaseOnly)
            {
                 // temp from inverse pwsat
                 temp_ = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
            }
            else if (phasePresence == gnPhaseOnly)
            {
                 // temp from inverse pnsat
                 temp_ = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(temp_);

            fluidState_.setTemperature(temp_);

            // now comes the tricky part: calculate phase composition
            if (phasePresence == threePhases) {

                // all phases are present, phase compositions are a
                // result of the the gas <-> liquid equilibrium.
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  pg_ - partPressH2O;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Henry
                Scalar xwn = partPressNAPL / FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                // Not yet filled with real numbers for the NAPL phase
                Scalar xnw = partPressH2O / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else if (phasePresence == wPhaseOnly) {
                // only the water phase is present, water phase composition is
                // stored explicitly.

                // extract mole fractions in the water phase
                Scalar xwn = priVars[switch2Idx];
                Scalar xww = 1 - xwn;

                // write water mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xgn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;
                Scalar xgw = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;


                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                // maybe solubility would be better than this approach via Henry
                Scalar xnn = xwn * FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx) / (xgn * pg_);
                Scalar xnw = xgw*pg_ / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);

                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else if (phasePresence == gnPhaseOnly) {

                // only gas and NAPL phases are present

                Scalar xnw = priVars[switch2Idx];
                Scalar xnn = 1.-xnw;
                Scalar xgn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;
                Scalar xgw = 1.-xgn;

                // note that the water phase is actually not present
                // the values are used as switching criteria
                Scalar xww = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);

            }
            else if (phasePresence == wnPhaseOnly) {
                // water and NAPL are present, phase compositions are a
                // mole fractions of non-existing gas phase are used as switching criteria
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Henry
                Scalar xwn = partPressNAPL / FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                Scalar xnw = partPressH2O / FluidSystem::henryCoefficient(fluidState_, nPhaseIdx,wCompIdx);
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else if (phasePresence == gPhaseOnly) {
                // only the gas phase is present, gas phase composition is
                // stored explicitly here below.

                const Scalar xgn = priVars[switch2Idx];
                Scalar xgw = 1 - xgn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);

                // note that the water and NAPL phase is actually not present
                // the values are used as switching criteria
                Scalar xww = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx) / pg_;
                Scalar xnn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else if (phasePresence == wgPhaseOnly) {
                // only water and gas phases are present
                const Scalar xgn = priVars[switch2Idx];
                Scalar xgw = 1 - xgn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);


                Scalar xwn = xgn*pg_/FluidSystem::henryCoefficient(fluidState_, wPhaseIdx,nCompIdx);
                Scalar xww = 1.-xwn;

                // write mole fractions in the fluid state
                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnn = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx) / pg_;

                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else
                assert(false); // unhandled phase state
        } // end of if(!UseSimpleModel), i.e. the more complex version with six phase states
        else // use the simpler model with only two phase states
        {
            /* first the saturations */
            if (phasePresence == threePhases)
            {
                sw_ = priVars[switch1Idx];
                sn_ = priVars[switch2Idx];
                sg_ = 1. - sw_ - sn_;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                sn_ = priVars[switch2Idx];
                sw_ = 1. - sn_;
                sg_ = 0.;
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(sg_);

            fluidState_.setSaturation(wPhaseIdx, sw_);
            fluidState_.setSaturation(gPhaseIdx, sg_);
            fluidState_.setSaturation(nPhaseIdx, sn_);

            /* now the pressures */
            if (phasePresence == threePhases)
            {
                 pg_ = priVars[pressureIdx];

                 // calculate capillary pressures
                 Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
                 Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
                 Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

                 Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
                 Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

                 pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
                 pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;
            }
            else if (phasePresence == wnPhaseOnly)
            {
                 pw_ = priVars[pressureIdx];

                 // calculate capillary pressures
                 Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
                 Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
                 Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

                 Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
                 Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

                 pn_ = pw_ + pcAlpha * pcnw + (1.-pcAlpha)*pcNW1;
                 pg_ = pn_ + pcAlpha * pcgn + (1.-pcAlpha)*(pcgw - pcNW1);
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(pw_);

            fluidState_.setPressure(wPhaseIdx, pw_);
            fluidState_.setPressure(gPhaseIdx, pg_);
            fluidState_.setPressure(nPhaseIdx, pn_);

            /* now the temperature */
            if (phasePresence == wnPhaseOnly)
            {
                 temp_ = priVars[switch1Idx];
            }
            else if (phasePresence == threePhases)
            {
                 if(sn_<=1.e-10) // this threshold values is chosen arbitrarily as a small number
                 {
                     Scalar tempOnlyWater = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
                     temp_ = tempOnlyWater;
                 }
                 if(sw_<=1.e-10) // this threshold values is chosen arbitrarily as a small number
                 {
                     Scalar tempOnlyNAPL = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
                     temp_ = tempOnlyNAPL;
                 }
                 else
                 {
                     // temp from inverse pwsat and pnsat which have to sum up to pg
                     Scalar tempOnlyNAPL = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, nCompIdx);
                     Scalar tempOnlyWater = FluidSystem::inverseVaporPressureCurve(fluidState_, gPhaseIdx, wCompIdx);
                     fluidState_.setTemperature(tempOnlyWater);
                     Scalar defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                         - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                     Scalar temp = tempOnlyWater; // initial guess
                     int counter = 0;
                     while(std::abs(defect) > 0.01) // simply a small number chosen ...
                     {
                         Scalar deltaT = 1.e-6; // fixed number, but T should always be in the order of a few hundred Kelvin
                         fluidState_.setTemperature(temp+deltaT);
                         Scalar fUp = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                          - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                         fluidState_.setTemperature(temp-deltaT);
                         Scalar fDown = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                          - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);

                         temp = temp - defect * 2. * deltaT / (fUp - fDown);

                         fluidState_.setTemperature(temp);
                         defect = pg_ - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx)
                                      - FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);
                         counter +=1;
                         if (counter>10) break;
                     }
                     if ((sw_>1.e-10)&&(sw_<0.01))
                         temp = temp + (sw_ - 1.e-10) * (temp - tempOnlyNAPL) / (0.01 - 1.e-10);
                     if ((sn_>1.e-10)&&(sn_<0.01))
                         temp = temp + (sn_ - 1.e-10) * (temp - tempOnlyWater) / (0.01 - 1.e-10);
                     temp_ = temp;
                 }
            }
            else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
            Valgrind::CheckDefined(temp_);

            fluidState_.setTemperature(temp_);

            // now comes the tricky part: calculate phase composition
            if (phasePresence == threePhases) {

                // all phases are present, phase compositions are a
                // result of the the gas <-> liquid equilibrium.
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  pg_ - partPressH2O;
                // regularized evaporation for small liquid phase saturations
                // avoids negative saturations of liquid phases
                if (sw_<0.02) partPressH2O *= sw_/0.02;
                if (partPressH2O < 0.) partPressH2O = 0;
                if (sn_<0.02) partPressNAPL *= sn_ / 0.02;
                if (partPressNAPL < 0.) partPressNAPL = 0;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // Immiscible liquid phases, mole fractions are just dummy values
                Scalar xwn = 0;
                Scalar xww = 1.-xwn;

                // Not yet filled with real numbers for the NAPL phase
                Scalar xnw = 0;
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else if (phasePresence == wnPhaseOnly) {
                // mole fractions of non-existing gas phase are used as switching criteria
                Scalar partPressH2O = FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, wCompIdx);
                Scalar partPressNAPL =  FluidSystem::partialPressureGas(fluidState_, gPhaseIdx, nCompIdx);;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;

                // immiscible liquid phases, mole fractions are just dummy values
                Scalar xwn = 0;
                Scalar xww = 1.-xwn;

                Scalar xnw = 0;
                Scalar xnn = 1.-xnw;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
            else
                assert(false); // unhandled phase state
            }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);

            Scalar kr;
            kr = MaterialLaw::kr(materialParams, phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx),
                                 fluidState_.saturation(gPhaseIdx));
            mobility_[phaseIdx] = kr / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // material dependent parameters for NAPL adsorption
        bulkDensTimesAdsorpCoeff_ =
            MaterialLaw::bulkDensTimesAdsorpCoeff(materialParams);

        /* ATTENTION: The conversion to effective diffusion parameters
         *            for the porous media happens at another place!
         */

        // diffusivity coefficents
        diffusionCoefficient_[gPhaseIdx] = FluidSystem::diffusionCoefficient(fluidState_, gPhaseIdx);

        diffusionCoefficient_[wPhaseIdx] = FluidSystem::diffusionCoefficient(fluidState_, wPhaseIdx);

        /* no diffusion in NAPL phase considered  at the moment, dummy values */
        diffusionCoefficient_[nPhaseIdx] = 1.e-10;

        Valgrind::CheckDefined(diffusionCoefficient_);

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                         fvGeometry,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        // permeability
        permeability_ = problem.spatialParams().intrinsicPermeability(element,
                                                                          fvGeometry,
                                                                          scvIdx);
        Valgrind::CheckDefined(permeability_);

        // the enthalpies (internal energies are directly calculated in the fluidstate
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fluidState_, phaseIdx);
            fluidState_.setEnthalpy(phaseIdx, h);
        }


        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

     /*!
     * \brief Returns the mass fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mole fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    Scalar permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the diffusivity coefficient matrix
     */
    Dune::FieldVector<Scalar, numPhases> diffusionCoefficient() const
    { return diffusionCoefficient_; }

    /*!
     * \brief Returns the adsorption information
     */
    Scalar bulkDensTimesAdsorpCoeff() const
    { return bulkDensTimesAdsorpCoeff_; }

     /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return fluidState_.enthalpy(phaseIdx); };

protected:

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    { }

    Scalar sw_, sg_, sn_, pg_, pw_, pn_, temp_;

    Scalar moleFrac_[numPhases][numComponents];
    Scalar massFrac_[numPhases][numComponents];

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar permeability_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    Scalar bulkDensTimesAdsorpCoeff_; //!< the basis for calculating adsorbed NAPL
    /* We need a tensor here !! */
    //!< Binary diffusion coefficients of the 3 components in the phases
    Dune::FieldVector<Scalar, numPhases> diffusionCoefficient_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

template <class TypeTag>
const typename ThreePWaterOilVolumeVariables<TypeTag>::Scalar ThreePWaterOilVolumeVariables<TypeTag>::R = Constants<typename GET_PROP_TYPE(TypeTag, Scalar)>::R;

} // end namespace

#endif
