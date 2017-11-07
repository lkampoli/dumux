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
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
#ifndef DUMUX_FUELCELL_PROBLEM_HH
#define DUMUX_FUELCELL_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/2pnc/implicit/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>

#include "fuelcellspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class FuelCellProblem;

namespace Properties
{
NEW_TYPE_TAG(FuelCellProblem, INHERITS_FROM(TwoPNC, FuelCellSpatialParams));
NEW_TYPE_TAG(FuelCellBoxProblem, INHERITS_FROM(BoxModel, FuelCellProblem));
NEW_TYPE_TAG(FuelCellCCProblem, INHERITS_FROM(CCTpfaModel, FuelCellProblem));

// Set the grid type
SET_TYPE_PROP(FuelCellProblem, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(FuelCellProblem, Problem, FuelCellProblem<TypeTag>);
// Set the primary variable combination for the 2pnc model
SET_INT_PROP(FuelCellProblem, Formulation, TwoPNCFormulation::pnsw);

// Set fluid configuration
SET_PROP(FuelCellProblem, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    static const bool useComplexRelations = true;
public:
    using type = FluidSystems::H2ON2O2<Scalar, useComplexRelations>;
};

// Set the transport equation that is replaced by the total mass balance
SET_INT_PROP(FuelCellProblem, ReplaceCompEqIdx, 3);
}


/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class FuelCellProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Sources = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);




    // Select the electrochemistry method
    using ElectroChemistry = typename Dumux::ElectroChemistry<TypeTag, ElectroChemistryModel::Ochs>;

    enum
    {
        numComponents = FluidSystem::numComponents,
        numSecComponents = FluidSystem::numSecComponents,
    };
    // phase indices
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    // component indices
    enum
    {
        wCompIdx = FluidSystem::wCompIdx, //major component of the liquid phase
        nCompIdx = FluidSystem::nCompIdx, //major component of the gas phase
    };
    // privar indices
    enum
    {
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        switchIdx = Indices::switchIdx, //liquid saturation or mole fraction
        conti0EqIdx = Indices::conti0EqIdx
    };

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { dofCodim = isBox ? dim : 0 };
public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FuelCellProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");
        temperature_        = getParam<Scalar>("Problem.InitialTemperature");

        name_               = getParam<std::string>("Problem.Name");

        pO2Inlet_           = getParam<Scalar>("ElectroChemistry.pO2Inlet");

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        currentDensity_.resize(fvGridGeometry->gridView().size(dofCodim));
        reactionSourceH2O_.resize(fvGridGeometry->gridView().size(dofCodim));
        reactionSourceO2_.resize(fvGridGeometry->gridView().size(dofCodim));
    }

    /*!
     * \name Problem parameters
     */

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    //! \copydoc Dumux::ImplicitProblem::source()
    Sources source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        Sources values(0.0);
        const auto& globalPos = scv.dofPosition();

        //reaction sources from electro chemistry
        if(inReactionLayer_(globalPos))
        {
            const auto& volVars = elemVolVars[scv];
            auto currentDensity = ElectroChemistry::calculateCurrentDensity(volVars);
            ElectroChemistry::reactionSource(values, currentDensity);
        }

        return values;
    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();

        if (onUpperBoundary_(globalPos)){
            bcTypes.setAllDirichlet();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        auto priVars = initial_(globalPos);

        if(onUpperBoundary_(globalPos))
        {
            Scalar pg = 1.0e5;
            priVars[pressureIdx] = pg;
            priVars[switchIdx] = 0.3;//Sl for bothPhases
            priVars[switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
        }

        return priVars;
    }

    /*!
     * \name Volume terms
     */


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }


    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    const std::vector<Scalar>& getCurrentDensity()
    {
        return currentDensity_;
    }

    const std::vector<Scalar>& getReactionSourceH2O()
    {
        return reactionSourceH2O_;
    }

    const std::vector<Scalar>& getReactionSourceO2()
    {
        return reactionSourceO2_;
    }


    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            ElementSolutionVector elemSol(element, curSol, this->fvGridGeometry());

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto& globalPos = scv.dofPosition();
                const auto dofIdxGlobal = scv.dofIndex();

                if(inReactionLayer_(globalPos))
                {
                    //reactionSource Output
                    PrimaryVariables source;
                    auto i = ElectroChemistry::calculateCurrentDensity(volVars);
                    ElectroChemistry::reactionSource(source, i);

                    reactionSourceH2O_[dofIdxGlobal] = source[wPhaseIdx];
                    reactionSourceO2_[dofIdxGlobal] = source[numComponents-1];

                    //Current Output in A/cm^2
                    currentDensity_[dofIdxGlobal] = i/10000;
                }
                else
                {
                    reactionSourceH2O_[dofIdxGlobal] = 0.0;
                    reactionSourceO2_[dofIdxGlobal] = 0.0;
                    currentDensity_[dofIdxGlobal] = 0.0;
                }
            }
        }
    }



private:

    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::bothPhases);

        Scalar pg = 1.0e5;
        priVars[pressureIdx] = pg;
        priVars[switchIdx] = 0.3;//Sl for bothPhases
        priVars[switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases

        return priVars;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    bool inReactionLayer_(const GlobalPosition& globalPos) const
    { return globalPos[1] < 0.1*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]) + eps_; }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    int nTemperature_;
    int nPressure_;
    std::string name_ ;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar pO2Inlet_;
    std::vector<double> currentDensity_;
    std::vector<double> reactionSourceH2O_;
    std::vector<double> reactionSourceO2_;
};

} //end namespace Dumux

#endif
