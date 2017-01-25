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
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1947)
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>


#include <dumux/linear/amgbackend.hh>

// solve Navier-Stokes equations
#define ENABLE_NAVIERSTOKES 1


namespace Dumux
{
template <class TypeTag>
class KovasznayTestProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<KovasznayTestProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
NEW_TYPE_TAG(KovasznayTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

SET_PROP(KovasznayTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(KovasznayTestProblem, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(KovasznayTestProblem, Problem, Dumux::KovasznayTestProblem<TypeTag> );

SET_BOOL_PROP(KovasznayTestProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(KovasznayTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(KovasznayTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(KovasznayTestProblem, ProblemEnableGravity, true);

#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(KovasznayTestProblem, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(KovasznayTestProblem, EnableInertiaTerms, false);
#endif
}

/*!
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the staggered grid (Kovasznay 1947)
 */
template <class TypeTag>
class KovasznayTestProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    KovasznayTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView), eps_(1e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        kinematicViscosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, LiquidKinematicViscosity);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);

        using CellArray = std::array<unsigned int, dimWorld>;
        const CellArray numCells = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                      CellArray,
                                                      Grid,
                                                      Cells);
        cellSizeX_ = this->bBoxMax()[0] / numCells[0];
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }


    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (isLowerLeftCell_(globalPos))
            values.setDirichlet(massBalanceIdx);
        else
            values.setNeumann(massBalanceIdx);

        return values;
    }


    using ParentType::dirichletAtPos;
    BoundaryValues dirichletAtPos(const GlobalPosition & globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        BoundaryValues values;
        values[pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        values[velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
        values[velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    CellCenterPrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        return CellCenterPrimaryVariables(0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    using ParentType::initial;
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values;
        values[pressureIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
    {
        //Here we calculate the analytical solution
        const auto numElements = this->gridView().size(0);
        auto& pExact = *(this->resultWriter().allocateManagedBuffer(numElements));
        auto& velocityExact = *(this->resultWriter()).template allocateManagedBuffer<double, dimWorld>(numElements);

        using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
        typename DofTypeIndices::FaceIdx faceIdx;

        const auto someElement = *(elements(this->gridView()).begin());

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                const auto& globalPos = scv.dofPosition();

                pExact[dofIdxGlobal] = dirichletAtPos(globalPos)[pressureIdx];

                GlobalPosition velocityVector(0.0);
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    auto dirIdx = scvf.directionIndex();
                    velocityVector[dirIdx] += 0.5*dirichletAtPos(globalPos)[faceIdx][dirIdx];
                }
                velocityExact[dofIdxGlobal] = velocityVector;
            }
        }
        this->resultWriter().attachDofData(pExact, "p_exact", false);
        this->resultWriter().attachDofData(velocityExact,  "velocity_exact", false, dim);
    }

    /*!
     * \brief Calculate the L2 error between the analytical solution and the numerical approximation.
     *
     */
    BoundaryValues calculateL2Error()
    {
        BoundaryValues sumError(0.0), sumReference(0.0), l2NormAbs(0.0), l2NormRel(0.0);
        BoundaryValues analyticalSolution(0.0);
        BoundaryValues numericalSolution(0.0); // TODO type
        int numElements(0);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                    auto dofIdxGlobal = scv.dofIndex();
                    const auto& globalPos = scv.dofPosition();

                    analyticalSolution = dirichletAtPos(globalPos);

                    numericalSolution[pressureIdx] = this->model().curSol()[cellCenterIdx][pressureIdx][dofIdxGlobal]; // Multitypeblockvector
                    numericalSolution[velocityXIdx] = this->model().curSol()[faceIdx][velocityXIdx][dofIdxGlobal];
                    numericalSolution[velocityYIdx] = this->model().curSol()[faceIdx][velocityYIdx][dofIdxGlobal];

                    sumError[pressureIdx] += (analyticalSolution[pressureIdx] - numericalSolution[pressureIdx]) * (analyticalSolution[pressureIdx] - numericalSolution[pressureIdx]);
                    sumError[velocityXIdx] += (analyticalSolution[velocityXIdx] - numericalSolution[velocityXIdx]) * (analyticalSolution[velocityXIdx] - numericalSolution[velocityXIdx]);
                    sumError[velocityYIdx] += (analyticalSolution[velocityYIdx] - numericalSolution[velocityYIdx]) * (analyticalSolution[velocityYIdx] - numericalSolution[velocityYIdx]);

                    sumReference[pressureIdx] += (analyticalSolution[pressureIdx]) * (analyticalSolution[pressureIdx]);
                    sumReference[velocityXIdx] += (analyticalSolution[velocityXIdx]) * (analyticalSolution[velocityXIdx]);
                    sumReference[velocityYIdx] += (analyticalSolution[velocityYIdx]) * (analyticalSolution[velocityYIdx]);

                    ++numElements;
            }
        }

//    	l2NormAbs = std::sqrt(sumError / numElements);
        l2NormRel[pressureIdx] = std::sqrt(sumError[pressureIdx] / numElements / sumReference[pressureIdx]);
        l2NormRel[velocityXIdx] = std::sqrt(sumError[velocityXIdx] / numElements / sumReference[velocityXIdx]);
        l2NormRel[velocityYIdx] = std::sqrt(sumError[velocityYIdx] / numElements / sumReference[velocityYIdx]);

        return l2NormRel; // TODO return both errors??
    }

    /*!
     * \brief Write the L2 error into an output file
     *
     */
    void writeOutput(const bool verbose = true)
    {
        ParentType::writeOutput(verbose);

        BoundaryValues l2error = calculateL2Error();

        std::cout.precision(8);
        std::cout << "** L2 error for "
                    << std::setw(6) << this->gridView().size(0)
            " elements: "
            std::scientific
            "L2(p) = "
            l2error[pressureIdx]
            ", L2(vx) = "
            l2error[velocityXIdx]
            << ", L2(vy) = "
            l2error[velocityYIdx]
            std::endl;
    }

private:

    bool isLowerLeftCell_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < (this->bBoxMin()[0] + 0.5*cellSizeX_ + eps_);
    }

    Scalar eps_;
    Scalar cellSizeX_;
    std::string name_;

    Scalar kinematicViscosity_;
    Scalar lambda_;
};
} //end namespace

#endif
