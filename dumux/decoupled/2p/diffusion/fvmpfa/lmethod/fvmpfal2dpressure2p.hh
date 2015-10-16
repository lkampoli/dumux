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
#ifndef DUMUX_FVMPFAL2DPRESSURE2P_HH
#define DUMUX_FVMPFAL2DPRESSURE2P_HH

// dumux environment
#include <dumux/decoupled/common/fv/fvpressure.hh>
#include <dumux/decoupled/common/fv/mpfa/mpfalinteractionvolume.hh>
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include "fvmpfal2dtransmissibilitycalculator.hh"

/**
 * @file
 * @brief Finite volume MPFA L-method discretization of a two-phase pressure equation of the sequential IMPES model.
 */

namespace Dumux
{
//! \ingroup FVPressure2p
/*! \brief Finite volume MPFA L-method discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * Finite volume MPFA L-method discretization of the equations
 * \f[ - \text{div}\, \boldsymbol v_t = - \text{div}\, (\lambda_t \boldsymbol K \textbf{grad}\,
 * \Phi_w + f_n \lambda_t \boldsymbol K \textbf{grad}\, \Phi_{cap}   ) = 0, \f]
 * or
 * \f[ - \text{div}\, \boldsymbol v_t = - \text{div}\, (\lambda_t \boldsymbol K \textbf{grad}\,
 * \Phi_n - f_w \lambda_t \boldsymbol K \textbf{grad}\, \Phi_{cap}   ) = 0. \f]
 *  At Dirichlet boundaries a two-point flux approximation is used.
 * \f[ \Phi = g \;  \text{on} \; \Gamma_1, \quad \text{and} \quad
 * -\text{div}\, \boldsymbol v_t \cdot \mathbf{n} = J \;  \text{on}  \; \Gamma_2. \f]
 *
 *  Here, \f$ \Phi_\alpha \f$ denotes the potential of phase \f$ \alpha \f$, \f$ \boldsymbol K \f$ the intrinsic permeability,
 * \f$ \lambda_t \f$ the total mobility, \f$ f_\alpha \f$ the phase fractional flow function.
 *
 * More details on the equations can be found in
 *
 * Wolff 2013: http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/
 *
 * M. Wolff, Y. Cao, B. Flemisch, R. Helmig, and B. Wohlmuth (2013a). Multi-point flux
 * approximation L-method in 3D: numerical convergence and application to two-phase
 * flow through porous media. In P. Bastian, J. Kraus, R. Scheichl, and M. Wheeler,
 * editors, Simulation of Flow in Porous Media - Applications in Energy and Environment. De Gruyter.
 *
 * M. Wolff, B. Flemisch, R. Helmig, I. Aavatsmark.
 * Treatment of tensorial relative permeabilities with multipoint flux approximation.
 * International Journal of Numerical Analysis and Modeling (9), pp. 725-744, 2012.
 *
 *  Remark1: only for 2-D quadrilateral grid
 *
 *  Remark2: implemented for UGGrid, ALUGrid, or YaspGrid
 *
 *\tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class FvMpfaL2dPressure2p: public FVPressure<TypeTag>
{
    typedef FVPressure<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;

    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef typename GET_PROP_TYPE(TypeTag, GridTypeIndices) GridTypeIndices;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    enum
    {
        globalCorner = 2,
        globalEdge = 3,
        neumannNeumann = 0,
        dirichletDirichlet = 1,
        dirichletNeumann = 2,
        neumannDirichlet = 3
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*! \brief Type of the interaction volume objects
     *
     * Type of the interaction volume objects used to store the geometric information which is needed
     * to calculated the transmissibility matrices of one MPFA interaction volume.
     *
     */
    typedef Dumux::FVMPFALInteractionVolume<TypeTag> InteractionVolume;
    typedef Dumux::FvMpfaL2dTransmissibilityCalculator<TypeTag> TransmissibilityCalculator;
private:

    typedef std::vector<InteractionVolume> GlobalInteractionVolumeVector;
    typedef std::vector<Dune::FieldVector<bool, 2 * dim> > InnerBoundaryVolumeFaces;

    // helper function that finds the correct neighboring intersections
    Intersection getNextIntersection_(const Element&, const IntersectionIterator&);

    //initializes the matrix to store the system of equations
    friend class FVPressure<TypeTag>;
    void initializeMatrix();

    void storeInteractionVolumeInfo();

    //function which assembles the system of equations to be solved
    void assemble();

public:

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    /*! \brief Updates interaction volumes
     *
     * Globally rebuilds the MPFA interaction volumes.
     *
     */
    void updateInteractionVolumeInfo()
    {
        interactionVolumes_.clear();
        innerBoundaryVolumeFaces_.clear();

        interactionVolumes_.resize(problem_.gridView().size(dim));
        innerBoundaryVolumeFaces_.resize(problem_.gridView().size(0), Dune::FieldVector<bool, 2 * dim>(false));

        storeInteractionVolumeInfo();
    }

    /*! \brief Initializes the pressure model
     *
     * \copydetails FVPressure::initialize()
     */
    void initialize()
    {
        ParentType::initialize();

        const auto element = *problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
        fluidState.setTemperature(problem_.temperature(element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

        updateMaterialLaws();

        interactionVolumes_.resize(problem_.gridView().size(dim));
        innerBoundaryVolumeFaces_.resize(problem_.gridView().size(0), Dune::FieldVector<bool, 2 * dim>(false));

        storeInteractionVolumeInfo();

        assemble();
        this->solve();

        storePressureSolution();

        return;
    }

    /*! \brief Globally stores the pressure solution
     *
     */
    void storePressureSolution()
    {
        for (const auto& element : Dune::elements(problem_.gridView()))
        {
            storePressureSolution(element);
        }
    }

    /*! \brief Stores the pressure solution of a cell
     *
     * \param element Dune grid element
     */
    void storePressureSolution(const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        switch (pressureType_)
        {
        case pw:
        {
            Scalar potW = this->pressure()[eIdxGlobal];

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potPc = cellData.capillaryPressure()
                    + gravityDiff * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potW + potPc);

            Scalar pressW = potW - gravityDiff * density_[wPhaseIdx];

            cellData.setPressure(wPhaseIdx, pressW);
            cellData.setPressure(nPhaseIdx, pressW + cellData.capillaryPressure());

            break;
        }
        case pn:
        {
            Scalar potNw = this->pressure()[eIdxGlobal];

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potPc = cellData.capillaryPressure()
                    + gravityDiff * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            cellData.setPotential(nPhaseIdx, potNw);
            cellData.setPotential(wPhaseIdx, potNw - potPc);

            Scalar pressNw = potNw - gravityDiff * density_[nPhaseIdx];

            cellData.setPressure(wPhaseIdx, pressNw - cellData.capillaryPressure());
            cellData.setPressure(nPhaseIdx, pressNw);

            break;
        }
        }
        cellData.fluxData().resetVelocity();
    }

    /*! \brief Pressure update
     *
     * \copydetails FVPressure::update()
     *
     */
    void update()
    {
        //error bounds for error term for incompressible models to correct unphysical saturation
        //over/undershoots due to saturation transport
        timeStep_ = problem_.timeManager().timeStepSize();
        maxError_ = 0.0;
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            Scalar sat = 0;
            switch (saturationType_)
            {
            case sw:
                sat = problem_.variables().cellData(i).saturation(wPhaseIdx);
                break;
            case sn:
                sat = problem_.variables().cellData(i).saturation(nPhaseIdx);
                break;
            }
            if (sat > 1.0)
            {
                maxError_ = std::max(maxError_, (sat - 1.0) / timeStep_);
            }
            if (sat < 0.0)
            {
                maxError_ = std::max(maxError_, (-sat) / timeStep_);
            }
        }

        assemble();
        this->solve();

        storePressureSolution();

        return;
    }


    /*! \brief Adds pressure output to the output file
     *
     * Adds the potential, the potential and the capillary pressure to the output.
     * If the VtkOutputLevel is equal to zero (default) only primary variables are written,
     * if it is larger than zero also secondary variables are written.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        ScalarSolutionType *potential = writer.allocateManagedBuffer(size);

        (*potential) = this->pressure();

        if (pressureType_ == pw)
        {
            writer.attachCellData(*potential, "wetting potential");
        }

        if (pressureType_ == pn)
        {
            writer.attachCellData(*potential, "nonwetting potential");
        }

        if (vtkOutputLevel_ > 0)
        {
            ScalarSolutionType *pressure = writer.allocateManagedBuffer(size);
            ScalarSolutionType *pressureSecond = writer.allocateManagedBuffer(size);
            ScalarSolutionType *potentialSecond = writer.allocateManagedBuffer(size);
            ScalarSolutionType *pc = writer.allocateManagedBuffer(size);

            for (const auto& element : Dune::elements(problem_.gridView()))
            {
                int idx = problem_.variables().index(element);
                CellData& cellData = problem_.variables().cellData(idx);

                (*pc)[idx] = cellData.capillaryPressure();

                if (pressureType_ == pw)
                {
                    (*pressure)[idx] = cellData.pressure(wPhaseIdx);
                    (*potentialSecond)[idx] = cellData.potential(nPhaseIdx);
                    (*pressureSecond)[idx] = cellData.pressure(nPhaseIdx);
                }

                if (pressureType_ == pn)
                {
                    (*pressure)[idx] = cellData.pressure(nPhaseIdx);
                    (*potentialSecond)[idx] = cellData.potential(wPhaseIdx);
                    (*pressureSecond)[idx] = cellData.pressure(wPhaseIdx);
                }
            }

            if (pressureType_ == pw)
            {
                writer.attachCellData(*pressure, "wetting pressure");
                writer.attachCellData(*pressureSecond, "nonwetting pressure");
                writer.attachCellData(*potentialSecond, "nonwetting potential");
            }

            if (pressureType_ == pn)
            {
                writer.attachCellData(*pressure, "nonwetting pressure");
                writer.attachCellData(*pressureSecond, "wetting pressure");
                writer.attachCellData(*potentialSecond, "wetting potential");
            }

            writer.attachCellData(*pc, "capillary pressure");
        }

        return;
    }

    //! Constructs a FvMpfaL2dPressure2p object
    /**
     * \param problem A problem class object
     */
    FvMpfaL2dPressure2p(Problem& problem) :
            ParentType(problem), problem_(problem), transmissibilityCalculator_(problem),
                    gravity_(problem.gravity()),
                    maxError_(0.), timeStep_(1.)
    {
        if (pressureType_ != pw && pressureType_ != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (GET_PROP_VALUE(TypeTag, EnableCompressibility))
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported!");
        }
        if (dim != 2)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }

        ErrorTermFactor_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, ErrorTermFactor);
        ErrorTermLowerBound_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, ErrorTermLowerBound);
        ErrorTermUpperBound_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, ErrorTermUpperBound);

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, Vtk, OutputLevel);
    }

private:
    Problem& problem_;
    TransmissibilityCalculator transmissibilityCalculator_;

protected:

    GlobalInteractionVolumeVector interactionVolumes_;//!< Global Vector of interaction volumes
    InnerBoundaryVolumeFaces innerBoundaryVolumeFaces_; //!< Vector marking faces which intersect the boundary

private:
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static constexpr Scalar threshold_ = 1e-15;
    //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
    //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    //! gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);

    /* Volume correction term to correct for unphysical saturation overshoots/undershoots.
     * These can occur if the estimated time step for the explicit transport was too large.
     * Correction by an artificial source term allows to correct this errors due to wrong
     * time-stepping without losing mass conservation. The error term looks as follows:
     * \f[
     *  q_{error} = \begin{cases}
     *          S < 0 &  \frac{S}{\Delta t} V \\
     *          S > 1 &  \frac{(S - 1)}{\Delta t} V \\
     *          0 \le S \le 1 & 0
     *      \end{cases}
     *  \f]
     *
     *  \param cellData The IMPES <tt>CellData</tt> object of the current cell.
     *
     *  \return The scalar value of the error term.
    */
    Scalar evaluateErrorTerm_(CellData& cellData)
    {
        //error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        // error reduction routine: volumetric error is damped and inserted to right hand side
        Scalar sat = 0;
        switch (saturationType_)
        {
        case sw:
            sat = cellData.saturation(wPhaseIdx);
            break;
        case sn:
            sat = cellData.saturation(nPhaseIdx);
            break;
        }

        Scalar error = (sat > 1.0) ? sat - 1.0 : 0.0;
        if (sat < 0.0)
        {
            error = sat;
        }
        error /= timeStep_;

        Scalar errorAbs = std::abs(error);

        if ((errorAbs * timeStep_ > 1e-6) && (errorAbs > ErrorTermLowerBound_ * maxError_)
                && (!problem_.timeManager().willBeFinished()))
        {
            return ErrorTermFactor_ * error;
        }
        return 0.0;
    }

};

template<class TypeTag>
typename FvMpfaL2dPressure2p<TypeTag>::Intersection
  FvMpfaL2dPressure2p<TypeTag>::getNextIntersection_(const Element& element,
                                                     const IntersectionIterator& isIt)
{
    auto isItBegin = problem_.gridView().ibegin(element);
    const auto isEndIt = problem_.gridView().iend(element);

    auto tempIsIt = isIt;
    auto nextIsIt = ++tempIsIt;

    // get 'nextIsIt'
    switch (GET_PROP_VALUE(TypeTag, GridImplementation))
    {
        // for YaspGrid
        case GridTypeIndices::yaspGrid:
        {
            if (nextIsIt == isEndIt)
            {
                nextIsIt = isItBegin;
            }
            else
            {
                nextIsIt = ++tempIsIt;

                if (nextIsIt == isEndIt)
                {
                    auto tempIsItBegin = isItBegin;
                    nextIsIt = ++tempIsItBegin;
                }
            }

            break;
        }
        // for ALUGrid and UGGrid
        case GridTypeIndices::aluGrid:
        case GridTypeIndices::ugGrid:
        {
            if (nextIsIt == isEndIt)
                nextIsIt = isItBegin;

            break;
        }
        default:
        {
            DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
            break;
        }
    }

    return *nextIsIt;
}

template<class TypeTag>
void FvMpfaL2dPressure2p<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        const auto isEndIt = problem_.gridView().iend(element);
        for (auto isIt = problem_.gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;
            auto nextIntersection = getNextIntersection_(element, isIt);

            if (intersection.neighbor())
                rowSize++;

            // try to find the common neighbor of intersection's and nextIntersection's outside
            if (intersection.neighbor() && nextIntersection.neighbor())
            {
                for (const auto& innerIntersection
                     : Dune::intersections(problem_.gridView(), *intersection.outside()))
                    for (const auto& innerNextIntersection
                         : Dune::intersections(problem_.gridView(), *nextIntersection.outside()))
                    {
                        if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                        {
                            if (innerIntersection.outside() == innerNextIntersection.outside()
                                && innerIntersection.outside() != intersection.inside())
                            {
                                rowSize++;
                            }
                        }
                    }
            }
        } // end of 'for' IntersectionIterator

        // set number of indices in row eIdxGlobalI to rowSize
        this->A_.setrowsize(eIdxGlobalI, rowSize);

    } // end of 'for' ElementIterator

    // indicate that size of all rows is defined
    this->A_.endrowsizes();

    // determine position of matrix entries
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // add diagonal index
        this->A_.addindex(eIdxGlobalI, eIdxGlobalI);

        // run through all intersections with neighbors
        const auto isEndIt = problem_.gridView().iend(element);
        for (auto isIt = problem_.gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;
            auto nextIntersection = getNextIntersection_(element, isIt);

            if (intersection.neighbor())
            {
                // access neighbor
                int eIdxGlobalJ = problem_.variables().index(intersection.outside());

                // add off diagonal index
                // add index (row,col) to the matrix
                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);
            }

            if (intersection.neighbor() && nextIntersection.neighbor())
            {
                for (const auto& innerIntersection
                     : Dune::intersections(problem_.gridView(), *intersection.outside()))
                    for (const auto& innerNextIntersection
                         : Dune::intersections(problem_.gridView(), *nextIntersection.outside()))
                    {
                        if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                        {
                            auto innerOutside = innerIntersection.outside();

                            if (innerOutside == innerNextIntersection.outside()
                                && innerOutside != intersection.inside())
                            {
                                int eIdxGlobalJ = problem_.variables().index(innerOutside);

                                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);
                            }
                        }
                    }
            }
        } // end of 'for' IntersectionIterator
    } // end of 'for' ElementIterator

    // indicate that all indices are defined, check consistency
    this->A_.endindices();

    return;
}
//                 Indices used in a interaction volume of the MPFA-o method
//                 ___________________________________________________
//                 |                        |                        |
//                 | nuxy: cell geometry    |       nxy: face normal |
//                 |     vectors (see MPFA) |                        |
//                 |                        |                        |
//                 |            4-----------3-----------3            |
//                 |            | --> nu43  |  nu34 <-- |            |
//                 |            | |nu41    1|--> n43   ||nu32        |
//                 |            | v   ^     |0     ^   v|            |
//                 |____________4__0__|n14__|__n23_|_1__2____________|
//                 |            |    1    0 |     0     |
//                 |            | ^         |1   nu23 ^ |            |
//                 |            | |nu14    0|--> n12  | |            |
//                 |            | -->nu12   |   nu21<-- |            |
//                 |            1-----------1-----------2            |
//                 |          elementnumber |inter-                  |
//                 |                        |face-                   |
//                 |                        |number                  |
//                 |________________________|________________________|

// only for 2-D general quadrilateral
template<class TypeTag>
void FvMpfaL2dPressure2p<TypeTag>::storeInteractionVolumeInfo()
{
    BoundaryTypes bcType;

    // run through all elements
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        // get common geometry information for the following computation

        int eIdxGlobal1 = problem_.variables().index(element);

        const auto isIt12End = problem_.gridView().iend(element);
        for (auto isIt12 = problem_.gridView().ibegin(element); isIt12 != isIt12End; ++isIt12)
        {
            const auto& intersection12 = *isIt12;
            auto intersection14 = getNextIntersection_(element, isIt12);

            int indexInInside12 = intersection12.indexInInside();
            int indexInInside14 = intersection14.indexInInside();

            // get the intersection node /bar^{x_3} between 'intersection12'
            // and 'intersection14', denoted as 'corner1234'
            const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());

            GlobalPosition corner1234(0);

            int globalVertIdx1234 = 0;

            // get the global coordinate and global vertex index of corner1234
            for (int i = 0; i < intersection12.geometry().corners(); ++i)
            {
                bool finished = false;

                int localVertIdx12corner = referenceElement.subEntity(indexInInside12, 1, i, dim);

                int globalVertIdx12corner = problem_.variables().vertexMapper().subIndex(element, localVertIdx12corner, dim);

                for (int j = 0; j < intersection14.geometry().corners(); ++j)
                {
                    int localVertIdx14corner = referenceElement.subEntity(indexInInside14, 1, j, dim);

                    int globalVertIdx14corner = problem_.variables().vertexMapper().subIndex(element, localVertIdx14corner, dim);

                    if (globalVertIdx12corner == globalVertIdx14corner)
                    {
                        corner1234 = element.geometry().corner(localVertIdx12corner);

                        globalVertIdx1234 = globalVertIdx12corner;

                        finished = true;
                        break;
                    }
                }

                if (finished)
                {
                    break;
                }
            }

            if (interactionVolumes_[globalVertIdx1234].isStored())
            {
                continue;
            }
            else
            {
                interactionVolumes_[globalVertIdx1234].setStored();
            }

            interactionVolumes_[globalVertIdx1234].setCenterPosition(corner1234);

            //store pointer 1
            interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection12.indexInInside(), 0, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInInside(), 0, 1);

            // center of face in global coordinates, i.e., the midpoint of edge 'intersection12'
            const GlobalPosition& globalPosFace12 = intersection12.geometry().center();

            // get face volume
            Scalar faceVol12 = intersection12.geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'intersection12'
            DimVector unitOuterNormal12 = intersection12.centerUnitOuterNormal();

            // center of face in global coordinates, i.e., the midpoint of edge 'intersection14'
            const GlobalPosition& globalPosFace41 = intersection14.geometry().center();

            // get face volume
            Scalar faceVol41 = intersection14.geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'intersection14': for numbering of n see Aavatsmark, Eigestad
            DimVector unitOuterNormal14 = intersection14.centerUnitOuterNormal();

            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 0, 1);
            //get the normals of from cell 2 and 4
            unitOuterNormal14 *= -1;
            unitOuterNormal12 *= -1;
            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);
            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 0, 1);

            // handle interior face
            if (intersection12.neighbor())
            {
                // access neighbor cell 2 of 'intersection12'
                auto element2 = intersection12.outside();

                int eIdxGlobal2 = problem_.variables().index(element2);

                //store pointer 2
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element2, 1);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection12.indexInOutside(), 1, 1);
                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 1, 1);
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 1, 1);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace12, 1, 1);

                // 'intersection14' is an interior face
                if (intersection14.neighbor())
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    auto element4 = intersection14.outside();

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element4, 3);
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInOutside(), 3, 0);

                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3, 0);

                    // cell 3
                    GlobalPosition globalPos3(0);

                    GlobalPosition globalPosFace23(0);
                    GlobalPosition globalPosFace34(0);

                    for (const auto& intersection23
                         : Dune::intersections(problem_.gridView(), element2))
                    {
                        bool finished = false;

                        for (const auto& intersection43
                             : Dune::intersections(problem_.gridView(), element4))
                        {
                            if (intersection23.neighbor() && intersection43.neighbor())
                            {
                                auto element32 = intersection23.outside();
                                auto element34 = intersection43.outside();

                                // find the common neighbor cell between cell 2 and cell 3, except cell 1
                                if (element32 == element34 && element32 != element)
                                {
                                    //store pointer 3
                                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element32, 2);

                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection23.indexInInside(), 1,
                                            0);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection23.indexInOutside(), 2,
                                            1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection43.indexInInside(), 3,
                                            1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection43.indexInOutside(), 2,
                                            0);

                                    // get global coordinate of neighbor cell 4 center
                                    globalPos3 = element32.geometry().center();

                                    globalPosFace23 = intersection23.geometry().center();
                                    globalPosFace34 = intersection43.geometry().center();

                                    Scalar faceVol23 = intersection23.geometry().volume() / 2.0;
                                    Scalar faceVol34 = intersection43.geometry().volume() / 2.0;

                                    // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                    DimVector unitOuterNormal23 = intersection23.centerUnitOuterNormal();

                                    DimVector unitOuterNormal43 = intersection43.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                    unitOuterNormal23 *= -1;
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                    unitOuterNormal43 *= -1;
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 3, 1);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                }

                // 'intersection14' is on the boundary
                else
                {
                    problem_.boundaryTypes(bcType, intersection14);
                    PrimaryVariables boundValues(0.0);

                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 3);
                    if (bcType.isNeumann(pressEqIdx))
                    {
                        problem_.neumann(boundValues, intersection14);
                        boundValues *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 3);
                    }
                    if (bcType.hasDirichlet())
                    {
                        problem_.dirichlet(boundValues, intersection14);
                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 3);
                    }

                    // get common geometry information for the following computation
                    // get the information of the face 'isIt23' between cell2 and cell4 (locally numbered)

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt23'
                    GlobalPosition globalPosFace23(0);

                    // get face volume
                    Scalar faceVol23 = 0;

                    // get outer normal vector scaled with half volume of face 'isIt23'
                    DimVector unitOuterNormal23(0);

                    bool finished = false;

                    for (const auto& intersection2
                         : Dune::intersections(problem_.gridView(), element2))
                    {
                        if (intersection2.boundary())
                        {
                            for (int i = 0; i < intersection2.geometry().corners(); ++i)
                            {
                                int localVertIdx2corner = referenceElement.subEntity(intersection2.indexInInside(), dim - 1, i,
                                        dim);

                                int globalVertIdx2corner = problem_.variables().index(
                                        *((element2).template subEntity < dim > (localVertIdx2corner)));

                                if (globalVertIdx2corner == globalVertIdx1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInInside(), 1,
                                            0);

                                    globalPosFace23 = intersection2.geometry().center();

                                    faceVol23 = intersection2.geometry().volume() / 2.0;

                                    unitOuterNormal23 = intersection2.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1, 0);

                                    problem_.boundaryTypes(bcType, intersection2);
                                    PrimaryVariables boundValues(0.0);

                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 1);
                                    if (bcType.isNeumann(pressEqIdx))
                                    {
                                        problem_.neumann(boundValues, intersection2);
                                        boundValues *= faceVol23;
                                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 1);
                                    }
                                    if (bcType.hasDirichlet())
                                    {
                                        problem_.dirichlet(boundValues, intersection2);
                                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 1);
                                    }

                                    interactionVolumes_[globalVertIdx1234].setOutsideFace(2);

                                    innerBoundaryVolumeFaces_[eIdxGlobal1][intersection12.indexInInside()] = true;
                                    innerBoundaryVolumeFaces_[eIdxGlobal2][intersection12.indexInOutside()] = true;

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(
                                Dune::NotImplemented,
                                "fvmpfao2pfaboundpressure2p.hh, l. 997: boundary shape not available as interaction volume shape");
                    }
                }
            }

            // handle boundary face 'intersection12'
            else
            {
                problem_.boundaryTypes(bcType, intersection12);
                PrimaryVariables boundValues(0.0);

                interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 0);
                if (bcType.isNeumann(pressEqIdx))
                {
                    problem_.neumann(boundValues, intersection12);
                    boundValues *= faceVol12;
                    interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 0);
                }
                if (bcType.hasDirichlet())
                {
                    problem_.dirichlet(boundValues, intersection12);
                    interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 0);
                }

                // 'intersection14' is on boundary
                if (intersection14.boundary())
                {
                    problem_.boundaryTypes(bcType, intersection14);
                    PrimaryVariables boundValues(0.0);

                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 3);
                    if (bcType.isNeumann(pressEqIdx))
                    {
                        problem_.neumann(boundValues, intersection14);
                        boundValues *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 3);
                    }
                    if (bcType.hasDirichlet())
                    {
                        problem_.dirichlet(boundValues, intersection14);
                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 3);
                    }

                    interactionVolumes_[globalVertIdx1234].setOutsideFace(1);
                    interactionVolumes_[globalVertIdx1234].setOutsideFace(2);
                }

                // 'intersection14' is inside
                else
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    auto element4 = intersection14.outside();
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInOutside(), 3, 0);

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element4, 3);

                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3, 0);

                    int eIdxGlobal4 = problem_.variables().index(element4);

                    bool finished = false;

                    // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                    for (const auto& intersection4
                         : Dune::intersections(problem_.gridView(), element4))
                    {
                        if (intersection4.boundary())
                        {
                            for (int i = 0; i < intersection4.geometry().corners(); ++i)
                            {
                                int localVertIdx4corner = referenceElement.subEntity(intersection4.indexInInside(), dim - 1, i,
                                        dim);

                                int globalVertIdx4corner = problem_.variables().index(
                                        *((element4).template subEntity < dim > (localVertIdx4corner)));

                                if (globalVertIdx4corner == globalVertIdx1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInInside(), 3,
                                            1);

                                    const GlobalPosition& globalPosFace34 = intersection4.geometry().center();

                                    Scalar faceVol34 = intersection4.geometry().volume() / 2.0;

                                    DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 3, 1);

                                    problem_.boundaryTypes(bcType, intersection4);
                                    PrimaryVariables boundValues(0.0);

                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 2);
                                    if (bcType.isNeumann(pressEqIdx))
                                    {
                                        problem_.neumann(boundValues, intersection4);
                                        boundValues *= faceVol34;
                                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 2);
                                    }
                                    if (bcType.hasDirichlet())
                                    {
                                        problem_.dirichlet(boundValues, intersection4);
                                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 2);
                                    }

                                    interactionVolumes_[globalVertIdx1234].setOutsideFace(1);

                                    innerBoundaryVolumeFaces_[eIdxGlobal1][intersection14.indexInInside()] = true;
                                    innerBoundaryVolumeFaces_[eIdxGlobal4][intersection14.indexInOutside()] = true;

                                    // get absolute permeability of neighbor cell 2
                                    DimMatrix K4(
                                            problem_.spatialParams().intrinsicPermeability(element4));

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(
                                Dune::NotImplemented,
                                "fvmpfao2pfaboundpressure2p.hh, l. 1164: boundary shape not available as interaction volume shape");
                    }
                }
            }

        } // end all intersections
    } // end grid traversal

    return;
}

// only for 2-D general quadrilateral
template<class TypeTag>
void FvMpfaL2dPressure2p<TypeTag>::assemble()
{
    // initialization: set global matrix this->A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    // run through all vertices
    for (const auto& vertex : Dune::vertices(problem_.gridView()))
    {
        int vIdxGlobal = problem_.variables().index(vertex);

        InteractionVolume& interactionVolume = interactionVolumes_[vIdxGlobal];

        if (interactionVolume.isInnerVolume())
        {

            auto element1 = interactionVolume.getSubVolumeElement(0);
            auto element2 = interactionVolume.getSubVolumeElement(1);
            auto element3 = interactionVolume.getSubVolumeElement(2);
            auto element4 = interactionVolume.getSubVolumeElement(3);

            // get global coordinate of cell centers
            const GlobalPosition& globalPos1 = element1.geometry().center();
            const GlobalPosition& globalPos2 = element2.geometry().center();
            const GlobalPosition& globalPos3 = element3.geometry().center();
            const GlobalPosition& globalPos4 = element4.geometry().center();

            // cell volumes
            Scalar volume1 = element1.geometry().volume();
            Scalar volume2 = element2.geometry().volume();
            Scalar volume3 = element3.geometry().volume();
            Scalar volume4 = element4.geometry().volume();

            // cell index
            int eIdxGlobal1 = problem_.variables().index(element1);
            int eIdxGlobal2 = problem_.variables().index(element2);
            int eIdxGlobal3 = problem_.variables().index(element3);
            int eIdxGlobal4 = problem_.variables().index(element4);

            //get the cell Data
            CellData& cellData1 = problem_.variables().cellData(eIdxGlobal1);
            CellData& cellData2 = problem_.variables().cellData(eIdxGlobal2);
            CellData& cellData3 = problem_.variables().cellData(eIdxGlobal3);
            CellData& cellData4 = problem_.variables().cellData(eIdxGlobal4);

            // evaluate right hand side
            PrimaryVariables source(0.0);
            problem_.source(source, element1);
            this->f_[eIdxGlobal1] += volume1 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element2);
            this->f_[eIdxGlobal2] += volume2 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element3);
            this->f_[eIdxGlobal3] += volume3 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element4);
            this->f_[eIdxGlobal4] += volume4 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

            this->f_[eIdxGlobal1] += evaluateErrorTerm_(cellData1) * volume1 / (4.0);
            this->f_[eIdxGlobal2] += evaluateErrorTerm_(cellData2) * volume2 / (4.0);
            this->f_[eIdxGlobal3] += evaluateErrorTerm_(cellData3) * volume3 / (4.0);
            this->f_[eIdxGlobal4] += evaluateErrorTerm_(cellData4) * volume4 / (4.0);

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda1(cellData1.mobility(wPhaseIdx));
            lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda2(cellData2.mobility(wPhaseIdx));
            lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda3(cellData3.mobility(wPhaseIdx));
            lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
            lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

            std::vector<DimVector > lambda(2*dim);

            lambda[0][0] = lambdaTotal1;
            lambda[0][1] = lambdaTotal1;
            lambda[1][0] = lambdaTotal2;
            lambda[1][1] = lambdaTotal2;
            lambda[2][0] = lambdaTotal3;
            lambda[2][1] = lambdaTotal3;
            lambda[3][0] = lambdaTotal4;
            lambda[3][1] = lambdaTotal4;

            //add capillary pressure and gravity terms to right-hand-side
            //calculate capillary pressure velocity
            Dune::FieldVector<Scalar, 2 * dim> pc(0);
            pc[0] = cellData1.capillaryPressure();
            pc[1] = cellData2.capillaryPressure();
            pc[2] = cellData3.capillaryPressure();
            pc[3] = cellData4.capillaryPressure();

            Dune::FieldVector<Scalar, 2 * dim> gravityDiff(0);

//            std::cout<<"maxPos = "<<problem_.bBoxMax()<<"\n";

            gravityDiff[0] = (problem_.bBoxMax() - globalPos1) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[1] = (problem_.bBoxMax() - globalPos2) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[2] = (problem_.bBoxMax() - globalPos3) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[3] = (problem_.bBoxMax() - globalPos4) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            pc += gravityDiff;

            Dune::FieldVector<Scalar, 2 * dim> pcFlux(0);

            Scalar pcPotential12 = 0;
            Scalar pcPotential14 = 0;
            Scalar pcPotential32 = 0;
            Scalar pcPotential34 = 0;

            DimVector Tu(0);
            Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);
            Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> T(0);

            int lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 2, 3);

            if (lType == TransmissibilityCalculator::rightTriangle)
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][0];
                this->A_[eIdxGlobal1][eIdxGlobal3] += T[1][1];
                this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][2];

                this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][0];
                this->A_[eIdxGlobal2][eIdxGlobal3] -= T[1][1];
                this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][2];

                u[0] = pc[1];
                u[1] = pc[2];
                u[2] = pc[0];

                T.mv(u, Tu);

                pcFlux[0] = Tu[1];
                pcPotential12 = Tu[1];
            }
            else
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][0];
                this->A_[eIdxGlobal1][eIdxGlobal4] += T[1][1];
                this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][2];

                this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][0];
                this->A_[eIdxGlobal2][eIdxGlobal4] -= T[1][1];
                this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][2];

                u[0] = pc[0];
                u[1] = pc[3];
                u[2] = pc[1];

                T.mv(u, Tu);

                pcFlux[0] = Tu[1];
                pcPotential12 = Tu[1];
            }

            lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 1, 2, 3, 0);

            if (lType == TransmissibilityCalculator::rightTriangle)
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal2][eIdxGlobal3] += T[1][0];
                this->A_[eIdxGlobal2][eIdxGlobal4] += T[1][1];
                this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][2];

                this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][0];
                this->A_[eIdxGlobal3][eIdxGlobal4] -= T[1][1];
                this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][2];

                u[0] = pc[2];
                u[1] = pc[3];
                u[2] = pc[1];

                T.mv(u, Tu);

                pcFlux[1] = Tu[1];
                pcPotential32 = -Tu[1];
            }
            else
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][0];
                this->A_[eIdxGlobal2][eIdxGlobal1] += T[1][1];
                this->A_[eIdxGlobal2][eIdxGlobal3] += T[1][2];

                this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][0];
                this->A_[eIdxGlobal3][eIdxGlobal1] -= T[1][1];
                this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][2];

                u[0] = pc[1];
                u[1] = pc[0];
                u[2] = pc[2];

                T.mv(u, Tu);

                pcFlux[1] = Tu[1];
                pcPotential32 = -Tu[1];
            }

            lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 2, 3, 0, 1);

            if (lType == TransmissibilityCalculator::rightTriangle)
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal3][eIdxGlobal4] += T[1][0];
                this->A_[eIdxGlobal3][eIdxGlobal1] += T[1][1];
                this->A_[eIdxGlobal3][eIdxGlobal3] += T[1][2];

                this->A_[eIdxGlobal4][eIdxGlobal4] -= T[1][0];
                this->A_[eIdxGlobal4][eIdxGlobal1] -= T[1][1];
                this->A_[eIdxGlobal4][eIdxGlobal3] -= T[1][2];

                u[0] = pc[3];
                u[1] = pc[0];
                u[2] = pc[2];

                T.mv(u, Tu);

                pcFlux[2] = Tu[1];
                pcPotential34 = Tu[1];
            }
            else
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal3][eIdxGlobal3] += T[1][0];
                this->A_[eIdxGlobal3][eIdxGlobal2] += T[1][1];
                this->A_[eIdxGlobal3][eIdxGlobal4] += T[1][2];

                this->A_[eIdxGlobal4][eIdxGlobal3] -= T[1][0];
                this->A_[eIdxGlobal4][eIdxGlobal2] -= T[1][1];
                this->A_[eIdxGlobal4][eIdxGlobal4] -= T[1][2];

                u[0] = pc[2];
                u[1] = pc[1];
                u[2] = pc[3];

                T.mv(u, Tu);

                pcFlux[2] = Tu[1];
                pcPotential34 = Tu[1];
            }

            lType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 3, 0, 1, 2);

            if (lType == TransmissibilityCalculator::rightTriangle)
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal4][eIdxGlobal1] += T[1][0];
                this->A_[eIdxGlobal4][eIdxGlobal2] += T[1][1];
                this->A_[eIdxGlobal4][eIdxGlobal4] += T[1][2];

                this->A_[eIdxGlobal1][eIdxGlobal1] -= T[1][0];
                this->A_[eIdxGlobal1][eIdxGlobal2] -= T[1][1];
                this->A_[eIdxGlobal1][eIdxGlobal4] -= T[1][2];

                u[0] = pc[0];
                u[1] = pc[1];
                u[2] = pc[3];

                T.mv(u, Tu);

                pcFlux[3] = Tu[1];
                pcPotential14 = -Tu[1];
            }
            else
            {
                if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                {
                    T *= 2;
                }
                this->A_[eIdxGlobal4][eIdxGlobal4] += T[1][0];
                this->A_[eIdxGlobal4][eIdxGlobal3] += T[1][1];
                this->A_[eIdxGlobal4][eIdxGlobal1] += T[1][2];

                this->A_[eIdxGlobal1][eIdxGlobal4] -= T[1][0];
                this->A_[eIdxGlobal1][eIdxGlobal3] -= T[1][1];
                this->A_[eIdxGlobal1][eIdxGlobal1] -= T[1][2];

                u[0] = pc[3];
                u[1] = pc[2];
                u[2] = pc[0];

                T.mv(u, Tu);

                pcFlux[3] = Tu[1];
                pcPotential14 = -Tu[1];
            }

            if (pc[0] == 0 && pc[1] == 0 && pc[2] == 0 && pc[3] == 0)
            {
                continue;
            }

            //compute mobilities of face 1
            Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
            lambda12Upw[wPhaseIdx] = (pcPotential12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda12Upw[nPhaseIdx] = (pcPotential12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 4
            Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
            lambda14Upw[wPhaseIdx] = (pcPotential14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda14Upw[nPhaseIdx] = (pcPotential14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

            //compute mobilities of face 2
            Dune::FieldVector<Scalar, numPhases> lambda32Upw(0.0);
            lambda32Upw[wPhaseIdx] = (pcPotential32 >= 0) ? lambda3[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda32Upw[nPhaseIdx] = (pcPotential32 >= 0) ? lambda3[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 3
            Dune::FieldVector<Scalar, numPhases> lambda34Upw(0.0);
            lambda34Upw[wPhaseIdx] = (pcPotential34 >= 0) ? lambda3[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda34Upw[nPhaseIdx] = (pcPotential34 >= 0) ? lambda3[nPhaseIdx] : lambda4[nPhaseIdx];

            for (int i = 0; i < numPhases; i++)
            {
                Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
                Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
                Scalar lambdaT32 = lambda32Upw[wPhaseIdx] + lambda32Upw[nPhaseIdx];
                Scalar lambdaT34 = lambda34Upw[wPhaseIdx] + lambda34Upw[nPhaseIdx];
                Scalar fracFlow12 = (lambdaT12 > threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
                Scalar fracFlow14 = (lambdaT14 > threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
                Scalar fracFlow32 = (lambdaT32 > threshold_) ? lambda32Upw[i] / (lambdaT32) : 0.0;
                Scalar fracFlow34 = (lambdaT34 > threshold_) ? lambda34Upw[i] / (lambdaT34) : 0.0;

                Dune::FieldVector<Scalar, 2 * dim> pcFluxReal(pcFlux);

                pcFluxReal[0] *= fracFlow12;
                pcFluxReal[1] *= fracFlow32;
                pcFluxReal[2] *= fracFlow34;
                pcFluxReal[3] *= fracFlow14;

//                if (std::isnan(pcFluxReal.two_norm()))
//                                std::cout<<"pcFlux = "<<pcFlux<<"\n";

                switch (pressureType_)
                {
                case pw:
                {
                    if (i == nPhaseIdx)
                    {
                        //add capillary pressure term to right hand side
                        this->f_[eIdxGlobal1] -= (pcFluxReal[0] - pcFluxReal[3]);
                        this->f_[eIdxGlobal2] -= (pcFluxReal[1] - pcFluxReal[0]);
                        this->f_[eIdxGlobal3] -= (pcFluxReal[2] - pcFluxReal[1]);
                        this->f_[eIdxGlobal4] -= (pcFluxReal[3] - pcFluxReal[2]);

                    }
                    break;
                }
                case pn:
                {
                    if (i == wPhaseIdx)
                    {
                        //add capillary pressure term to right hand side
                        this->f_[eIdxGlobal1] += (pcFluxReal[0] - pcFluxReal[3]);
                        this->f_[eIdxGlobal2] += (pcFluxReal[1] - pcFluxReal[0]);
                        this->f_[eIdxGlobal3] += (pcFluxReal[2] - pcFluxReal[1]);
                        this->f_[eIdxGlobal4] += (pcFluxReal[3] - pcFluxReal[2]);
                    }
                    break;
                }
                }
            }
        }

        // at least one face on boundary!
        else
        {
            for (int elemIdx = 0; elemIdx < 2 * dim; elemIdx++)
            {
                bool isOutside = false;
                for (int fIdx = 0; fIdx < dim; fIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);
                    if (interactionVolume.isOutsideFace(intVolFaceIdx))
                    {
                        isOutside = true;
                        break;
                    }
                }
                if (isOutside)
                {
                    continue;
                }

                auto element = interactionVolume.getSubVolumeElement(elemIdx);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos = element.geometry().center();

                // cell volumes
                Scalar volume = element.geometry().volume();

                // cell index
                int eIdxGlobal = problem_.variables().index(element);

                //get the cell Data
                CellData& cellData = problem_.variables().cellData(eIdxGlobal);

                //permeability vector at boundary
                DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(element));

                // evaluate right hand side
                PrimaryVariables source(0);
                problem_.source(source, element);
                this->f_[eIdxGlobal] += volume / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

                this->f_[eIdxGlobal] += evaluateErrorTerm_(cellData) * volume / (4.0);

                //get mobilities of the phases
                Dune::FieldVector<Scalar, numPhases> lambda(cellData.mobility(wPhaseIdx));
                lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

                Scalar pc = cellData.capillaryPressure();

                Scalar gravityDiff = (problem_.bBoxMax() - globalPos) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                pc += gravityDiff; //minus because of gravity definition!

                for (int fIdx = 0; fIdx < dim; fIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);

                    if (interactionVolume.isBoundaryFace(intVolFaceIdx))
                    {

                        if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressEqIdx))
                        {
                            int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                            const ReferenceElement& referenceElement = ReferenceElements::general(
                                    element.geometry().type());

                            const LocalPosition& localPos = referenceElement.position(boundaryFaceIdx, 1);

                            const GlobalPosition& globalPosFace = element.geometry().global(localPos);

                            DimVector distVec(globalPosFace - globalPos);
                            Scalar dist = distVec.two_norm();
                            DimVector unitDistVec(distVec);
                            unitDistVec /= dist;

                            Scalar faceArea = interactionVolume.getFaceArea(elemIdx, fIdx);

                            // get pc and lambda at the boundary
                            Scalar satWBound = cellData.saturation(wPhaseIdx);
                            //check boundary sat at face 1
                            if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(satEqIdx))
                            {
                                Scalar satBound = interactionVolume.getDirichletValues(intVolFaceIdx)[saturationIdx];
                                switch (saturationType_)
                                {
                                case sw:
                                {
                                    satWBound = satBound;
                                    break;
                                }
                                case sn:
                                {
                                    satWBound = 1 - satBound;
                                    break;
                                }
                                }

                            }

                            Scalar pcBound = MaterialLaw::pc(
                                    problem_.spatialParams().materialLawParams(element), satWBound);

                            Scalar gravityDiffBound = (problem_.bBoxMax() - globalPosFace) * gravity_
                                    * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                            pcBound += gravityDiffBound;

                            Dune::FieldVector<Scalar, numPhases> lambdaBound(
                                    MaterialLaw::krw(problem_.spatialParams().materialLawParams(element),
                                            satWBound));
                            lambdaBound[nPhaseIdx] = MaterialLaw::krn(
                                    problem_.spatialParams().materialLawParams(element), satWBound);
                            lambdaBound[wPhaseIdx] /= viscosity_[wPhaseIdx];
                            lambdaBound[nPhaseIdx] /= viscosity_[nPhaseIdx];

                            Scalar potentialBound = interactionVolume.getDirichletValues(intVolFaceIdx)[pressureIdx];
                            Scalar gdeltaZ = (problem_.bBoxMax()-globalPosFace) * gravity_;

                            //calculate potential gradients
                            Scalar potentialDiffW = 0;
                            Scalar potentialDiffNw = 0;
                            switch (pressureType_)
                            {
                            case pw:
                            {
                                potentialBound += density_[wPhaseIdx]*gdeltaZ;
                                potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBound) / dist;
                                potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBound - pcBound) / dist;
                                break;
                            }
                            case pn:
                            {
                                potentialBound += density_[nPhaseIdx]*gdeltaZ;
                                potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBound + pcBound) / dist;
                                potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBound) / dist;
                                break;
                            }
                            }

                            Scalar lambdaTotal = (potentialDiffW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                            lambdaTotal += (potentialDiffNw >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                            DimVector permTimesNormal(0);
                            permeability.mv(unitDistVec, permTimesNormal);

                            //calculate current matrix entry
                            Scalar entry = lambdaTotal * (unitDistVec * permTimesNormal) / dist * faceArea;

                            //calculate right hand side
                            Scalar pcFlux = 0;

                            switch (pressureType_)
                            {
                            case pw:
                            {
                                // calculate capillary pressure gradient
                                DimVector pcGradient = unitDistVec;
                                pcGradient *= (pc - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                pcFlux = 0.5 * (lambda[nPhaseIdx] + lambdaBound[nPhaseIdx])
                                        * (permTimesNormal * pcGradient) * faceArea;

                                break;
                            }
                            case pn:
                            {
                                // calculate capillary pressure gradient
                                DimVector pcGradient = unitDistVec;
                                pcGradient *= (pc - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                pcFlux = 0.5 * (lambda[wPhaseIdx] + lambdaBound[wPhaseIdx])
                                        * (permTimesNormal * pcGradient) * faceArea;

                                break;

                            }
                            }

                            // set diagonal entry and right hand side entry
                            this->A_[eIdxGlobal][eIdxGlobal] += entry;
                            this->f_[eIdxGlobal] += entry * potentialBound;

                            if (pc == 0 && pcBound == 0)
                            {
                                continue;
                            }

                            for (int i = 0; i < numPhases; i++)
                            {
                                switch (pressureType_)
                                {
                                case pw:
                                {
                                    if (i == nPhaseIdx)
                                    {
                                        //add capillary pressure term to right hand side
                                        this->f_[eIdxGlobal] -= pcFlux;
                                    }
                                    break;
                                }
                                case pn:
                                {
                                    if (i == wPhaseIdx)
                                    {
                                        //add capillary pressure term to right hand side
                                        this->f_[eIdxGlobal] += pcFlux;
                                    }
                                    break;
                                }
                                }
                            }

                        }
                        else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx))
                        {
                            Scalar J = interactionVolume.getNeumannValues(intVolFaceIdx)[wPhaseIdx] / density_[wPhaseIdx];
                            J += interactionVolume.getNeumannValues(intVolFaceIdx)[nPhaseIdx] / density_[nPhaseIdx];

                            this->f_[eIdxGlobal] -= J;
                        }
                        else
                        {
                            std::cout << "interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx)"
                                    << interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx) << "\n";
                            DUNE_THROW(Dune::NotImplemented,
                                    "No valid boundary condition type defined for pressure equation!");
                        }
                    }
                }
            }

        } // end boundaries

    } // end vertex iterator

    // only do more if we have more than one process
    if (problem_.gridView().comm().size() > 1)
    {
        // set ghost and overlap element entries
        for (const auto& element : Dune::elements(problem_.gridView()))
        {
            if (element.partitionType() == Dune::InteriorEntity)
                continue;

            // get the global index of the cell
            int eIdxGlobalI = problem_.variables().index(element);

            this->A_[eIdxGlobalI] = 0.0;
            this->A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            this->f_[eIdxGlobalI] = this->pressure()[eIdxGlobalI];
        }
    }

    return;
}

/*! \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 */
template<class TypeTag>
void FvMpfaL2dPressure2p<TypeTag>::updateMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        Scalar satW = cellData.saturation(wPhaseIdx);

        Scalar pc = MaterialLaw::pc(problem_.spatialParams().materialLawParams(element), satW);

        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), satW)
                / viscosity_[wPhaseIdx];
        Scalar mobilityNw = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), satW)
                / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNw);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));
    }
    return;
}

}
// end of Dune namespace
#endif
