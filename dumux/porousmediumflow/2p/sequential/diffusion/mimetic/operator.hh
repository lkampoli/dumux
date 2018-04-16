// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup SequentialTwoPModel
 * \brief An assembler for the Jacobian matrix based on mimetic FD.
 */
#ifndef DUMUX_MIMETICOPERATOR2P_HH
#define DUMUX_MIMETICOPERATOR2P_HH

#include "croperator.hh"
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/mimetic/properties.hh>

namespace Dumux
{
/*!
 * \brief Levelwise assembler
 *
 * \ingroup SequentialTwoPModel
 *
 * This class serves as a base class for local assemblers. It provides
 * space and access to the local stiffness matrix. The actual assembling is done
 * in a derived class via the virtual assemble method.
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class MimeticOperatorAssemblerTwoP: public CROperatorAssemblerTwoP<TypeTag>
{
    using ParentType = CROperatorAssemblerTwoP<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
    };
    using LocalStiffness = typename GET_PROP_TYPE(TypeTag, LocalStiffness);

    using Element = typename GridView::template Codim<0>::Entity;

    using CellData = typename GET_PROP_TYPE(TypeTag, CellData);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using SolutionTypes = typename GET_PROP(TypeTag, SolutionTypes);
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation),
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        saturationIdx = Indices::saturationIdx,
        satEqIdx = Indices::satEqIdx,
        pressureEqIdx = Indices::pressureEqIdx
    };

    using FieldVector = Dune::FieldVector<Scalar, dimWorld>;

public:

    MimeticOperatorAssemblerTwoP(const GridView& gridView) :
            ParentType(gridView)
    {
    }

    // TODO doc me!
    template<class Vector>
    void calculatePressure(LocalStiffness& loc, Vector& u, Problem& problem)
    {
        Dune::FieldVector<Scalar, 2 * dim> velocityW(0);
        Dune::FieldVector<Scalar, 2 * dim> velocityNw(0);
        Dune::FieldVector<Scalar, 2 * dim> pressTrace(0);
        Dune::FieldVector<Scalar, 2 * dim> gravPotTrace(0);

        const auto firstElement = *this->gridView_.template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem.referencePressure(firstElement));
        fluidState.setPressure(nPhaseIdx, problem.referencePressure(firstElement));
        fluidState.setTemperature(problem.temperature(firstElement));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        Scalar densityDiff = FluidSystem::density(fluidState, nPhaseIdx) - FluidSystem::density(fluidState, wPhaseIdx);
        Scalar viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
        Scalar viscosityNw = FluidSystem::viscosity(fluidState, nPhaseIdx);

        //reset velocity
        for (int i = 0; i < problem.gridView().size(0); i++)
        {
            problem.variables().cellData(i).fluxData().resetVelocity();
        }

        // run over all level elements
        for (const auto& element : elements(this->gridView_))
        {
            int eIdxGlobal = problem.variables().index(element);

            CellData& cellData = problem.variables().cellData(eIdxGlobal);
            FieldVector globalPos = element.geometry().center();

            // get local to global id map and pressure traces
            for (const auto& intersection : intersections(problem.gridView(), element))
            {
                int indexInInside = intersection.indexInInside();

                int fIdxGlobal = this->faceMapper_.subIndex(element, indexInInside, 1);

                pressTrace[indexInInside] = u[fIdxGlobal];

                gravPotTrace[indexInInside] = (problem.bBoxMax() - intersection.geometry().center()) * problem.gravity() * densityDiff;
            }

            switch (pressureType)
            {
            case pw:
            {
                Scalar potW = loc.constructPressure(element, pressTrace);
                Scalar gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * densityDiff;
                Scalar potNw = potW + gravPot;

                cellData.setPotential(wPhaseIdx, potW);
                cellData.setPotential(nPhaseIdx, potNw);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, wPhaseIdx);

                cellData.setPressure(wPhaseIdx, potW - gravPot);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, nPhaseIdx);

                cellData.setPressure(nPhaseIdx, potNw - gravPot);

                break;
            }
            case pn:
            {
                Scalar potNw = loc.constructPressure(element, pressTrace);
                Scalar  gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * densityDiff;
                Scalar potW = potNw - gravPot;

                cellData.setPotential(nPhaseIdx, potNw);
                cellData.setPotential(wPhaseIdx, potW);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, wPhaseIdx);

                cellData.setPressure(wPhaseIdx, potW - gravPot);

                gravPot = (problem.bBoxMax() - globalPos) * problem.gravity() * FluidSystem::density(fluidState, nPhaseIdx);

                cellData.setPressure(nPhaseIdx, potNw - gravPot);

                break;
            }
            }

            //velocity reconstruction: !!! The velocity which is not reconstructed from the primary
            //pressure variable can be slightly wrong and not conservative!!!!
            // -> Should not be used for transport!!
            switch (pressureType)
            {
            case pw:
            {
                loc.constructVelocity(element, velocityW, pressTrace, cellData.potential(wPhaseIdx));
                pressTrace += gravPotTrace;
                loc.constructVelocity(element, velocityNw, pressTrace, cellData.potential(nPhaseIdx));

                break;
            }
            case pn:
            {
                loc.constructVelocity(element, velocityW, pressTrace, cellData.potential(nPhaseIdx));
                pressTrace -= gravPotTrace;
                loc.constructVelocity(element, velocityNw, pressTrace, cellData.potential(wPhaseIdx));

                break;
            }
            }

            for (const auto& intersection : intersections(problem.gridView(), element))
            {
                int idxInInside = intersection.indexInInside();

                cellData.fluxData().setUpwindPotential(wPhaseIdx, idxInInside, velocityW[idxInInside]);
                cellData.fluxData().setUpwindPotential(nPhaseIdx, idxInInside, velocityNw[idxInInside]);

                Scalar mobilityW = 0;
                Scalar mobilityNw = 0;

                if (intersection.neighbor())
                {
                    int neighborIdx = problem.variables().index(intersection.outside());

                    CellData& cellDataNeighbor = problem.variables().cellData(neighborIdx);

                    mobilityW =
                            (velocityW[idxInInside] >= 0.) ? cellData.mobility(wPhaseIdx) :
                                    cellDataNeighbor.mobility(wPhaseIdx);
                    mobilityNw =
                            (velocityNw[idxInInside] >= 0.) ? cellData.mobility(nPhaseIdx) :
                                    cellDataNeighbor.mobility(nPhaseIdx);

                    if (velocityW[idxInInside] >= 0.)
                    {
                        FieldVector velocity(intersection.centerUnitOuterNormal());
                        velocity *= mobilityW/(mobilityW+mobilityNw) * velocityW[idxInInside];
                        cellData.fluxData().addVelocity(wPhaseIdx, idxInInside, velocity);
                        cellDataNeighbor.fluxData().addVelocity(wPhaseIdx, intersection.indexInOutside(), velocity);
                    }
                    if (velocityNw[idxInInside] >= 0.)
                    {
                        FieldVector velocity(intersection.centerUnitOuterNormal());
                        velocity *= mobilityNw/(mobilityW+mobilityNw) * velocityNw[idxInInside];
                        cellData.fluxData().addVelocity(nPhaseIdx, idxInInside, velocity);
                        cellDataNeighbor.fluxData().addVelocity(nPhaseIdx, intersection.indexInOutside(), velocity);
                    }

                    cellData.fluxData().setVelocityMarker(idxInInside);
                }
                else
                {
                    BoundaryTypes bctype;
                    problem.boundaryTypes(bctype, intersection);
                    if (bctype.isDirichlet(satEqIdx))
                    {
                        PrimaryVariables boundValues(0.0);
                        problem.dirichlet(boundValues, intersection);

                        if (velocityW[idxInInside] >= 0.)
                        {
                            mobilityW = cellData.mobility(wPhaseIdx);
                        }
                        else
                        {
                            mobilityW = MaterialLaw::krw(problem.spatialParams().materialLawParams(element),
                                    boundValues[saturationIdx]) / viscosityW;
                        }

                        if (velocityNw[idxInInside] >= 0.)
                        {
                            mobilityNw = cellData.mobility(nPhaseIdx);
                        }
                        else
                        {
                            mobilityNw = MaterialLaw::krn(problem.spatialParams().materialLawParams(element),
                                    boundValues[saturationIdx]) / viscosityNw;
                        }
                    }
                    else
                    {
                        mobilityW = cellData.mobility(wPhaseIdx);
                        mobilityNw = cellData.mobility(nPhaseIdx);
                    }

                    FieldVector velocity(intersection.centerUnitOuterNormal());
                    velocity *= mobilityW / (mobilityW + mobilityNw) * velocityW[idxInInside];
                    cellData.fluxData().setVelocity(wPhaseIdx, idxInInside, velocity);

                    velocity = 0;
                    velocity = intersection.centerUnitOuterNormal();
                    velocity *= mobilityNw / (mobilityW + mobilityNw) * velocityNw[idxInInside];
                    cellData.fluxData().setVelocity(nPhaseIdx, idxInInside, velocity);
                    cellData.fluxData().setVelocityMarker(idxInInside);
                }
            }
        }
    }
};
}
#endif
