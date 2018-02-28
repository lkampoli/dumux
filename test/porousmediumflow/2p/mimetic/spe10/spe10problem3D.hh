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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_MIMETIC_SPE10_HH
#define DUMUX_MIMETIC_SPE10_HH

#include "spe103DWetting.hh"
#include "spe103DNonwetting.hh"

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p/implicit/chopnewtoncontroller.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "spe10spatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class TwoPSpe10Problem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
#if PROBLEM==1
NEW_TYPE_TAG(TwoPSpe10Problem, INHERITS_FROM(CCTpfaModel, TwoP, Spe10SpatialParams));
#elif PROBLEM==2
NEW_TYPE_TAG(TwoPSpe10Problem, INHERITS_FROM(BoxModel, TwoP, Spe10SpatialParams));
#endif

SET_TYPE_PROP(TwoPSpe10Problem, Grid, Dune::YaspGrid<3>);

// Set the problem property
SET_TYPE_PROP(TwoPSpe10Problem, Problem, TwoPSpe10Problem<TypeTag>);

// Set the wetting phase
SET_PROP(TwoPSpe10Problem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SPE10Wetting<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TwoPSpe10Problem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SPE10Nonwetting<Scalar> > type;
};

SET_TYPE_PROP(TwoPSpe10Problem, SpatialParams, Spe10SpatialParams<TypeTag> );

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(TwoPSpe10Problem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

SET_TYPE_PROP(TwoPSpe10Problem, LinearSolver, AMGBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TwoPSpe10Problem, ProblemEnableGravity, false);

SET_BOOL_PROP(TwoPSpe10Problem, EnableGlobalFVGeometryCache, false);

SET_BOOL_PROP(TwoPSpe10Problem, EnableGlobalFluxVariablesCache, false);
SET_BOOL_PROP(TwoPSpe10Problem, EnableGlobalVolumeVariablesCache, false);

SET_TYPE_PROP(TwoPSpe10Problem, NewtonController, TwoPChopNewtonController<TypeTag> );
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class TwoPSpe10Problem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoPSpe10Problem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, EpisodeLength);

        this->timeManager().startNextEpisode(episodeLength_);

        pIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureIn);
        pOut_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureOut);

        IPoreVol_ = 0.0;
        P1PoreVol_ = 0.0;
        P2PoreVol_ = 0.0;
        P3PoreVol_ = 0.0;
        P4PoreVol_ = 0.0;

        xMin_[0] = 0;
        xMin_[1] = 0;
        xMin_[2] = 0;

        xMax_[0] = 365.76;
        xMax_[1] = 670.56;
        xMax_[2] = 51.816;
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        ParentType::init();

//        Scalar IX = 185.928;
//        Scalar IY = 336.804;
//
//        Scalar P1X = 0.0;
//        Scalar P1Y = 0.0;
//
//        Scalar P2X = 365.76;
//        Scalar P2Y = 0.0;
//
//        Scalar P3X = 365.76;
//        Scalar P3Y = 670.56;
//
//        Scalar P4X = 0.0;
//        Scalar P4Y = 670.56;
//
//        Scalar hX = 6.096;
//        Scalar hY = 3.048;
//        Scalar hZ = 0.6096;
//
//        for (const auto& element : elements(this->gridView()))
//        {
//            auto fvGeometry = localView(this->model().globalFvGeometry());
//            fvGeometry.bindElement(element);
//
//            for (auto&& scv : scvs(fvGeometry))
//            {
//                auto poro = this->spatialParams().porosity(element, scv,
//                        this->model().elementSolution(element, this->model().curSol()));
//
//                auto center = scv.center();
//
//                Scalar x = center[0];
//                Scalar y = center[1];
//                Scalar z = center[2];
//
//                if(std::abs(x-IX) < hX - eps_ && std::abs(y-IY) < hY - eps_)
//                {
//                    IPoreVol_ += poro*scv.volume();
//                }
//                else if(std::abs(x-P1X) < hX - eps_ && std::abs(y-P1Y) < hY - eps_)
//                {
//                    P1PoreVol_ += poro*scv.volume();
//                }
//                else if(std::abs(x-P2X) < hX - eps_ && std::abs(y-P2Y) < hY - eps_)
//                {
//                    P2PoreVol_ += poro*scv.volume();
//                }
//                else if(std::abs(x-P3X) < hX - eps_ && std::abs(y-P3Y) < hY - eps_)
//                {
//                    P3PoreVol_ += poro*scv.volume();
//                }
//                else if(std::abs(x-P4X) < hX - eps_ && std::abs(y-P4Y) < hY - eps_)
//                {
//                    P4PoreVol_ += poro*scv.volume();
//                }
//            }
//        }
//
//        std::cout << IPoreVol_ << std::endl;
//        std::cout << P1PoreVol_ << std::endl;
//        std::cout << P2PoreVol_ << std::endl;
//        std::cout << P3PoreVol_ << std::endl;
//        std::cout << P4PoreVol_ << std::endl;

        IPoreVol_ = 195.838;
        P1PoreVol_ = 147.786;
        P2PoreVol_ = 138.041;
        P3PoreVol_ = 175.53;
        P4PoreVol_ = 147.592;

        //readInitialValues();
    }

//    void readInitialValues()
//    {
//        std::ifstream file;
//        std::string iname = "spe10_tpfa_full3D";
//        iname.append("0_sn.txt");
//        file.open(iname, std::ios::in);
//        if (file.fail())
//            throw std::ios_base::failure(std::strerror(errno));
//
//        std::vector<Scalar> valSn;
//        int lineNumber = 0;
//        while (!file.eof()) {
//            Scalar val = 0;
//            file >> val;
//            valSn.push_back(val);
//            lineNumber++;
//         }
//
//        file.close();
//
//        iname = "spe10_tpfa_full3D";
//        iname.append("0_pw.txt");
//        file.open(iname, std::ios::in);
//        if (file.fail())
//            throw std::ios_base::failure(std::strerror(errno));
//
//        std::vector<Scalar> valpw;
//        lineNumber = 0;
//        while (!file.eof()) {
//            Scalar val = 0;
//            file >> val;
//            valpw.push_back(val);
//            lineNumber++;
//         }
//
//        file.close();
//        std::cout << "lineNumber: " << lineNumber << std::endl;
//        std::cout << "numEle: " << this->gridView().size(0) << std::endl;
//
//        auto& uCur = this->model().curSol();
//        for (const auto& element : elements(this->gridView()))
//        {
//            int eIdx = this->elementMapper().index(element);
//            int idxGlobal = this->spatialParams().getIndex_(element.geometry().center());
//            uCur[eIdx][snIdx] = valSn[idxGlobal];
//            uCur[eIdx][pwIdx] = valpw[idxGlobal];
//        }
//
//    }

    void readInitialValues()
    {
        std::ifstream file;
        std::string iname = "spe103D_vtu13_sw.txt";;
        file.open(iname, std::ios::in);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        std::vector<Scalar> valSw;
        int lineNumber = 0;
        while (!file.eof()) {
            Scalar val = 0;
            file >> val;
            valSw.push_back(val);
            lineNumber++;
         }

        file.close();

        iname = "spe103D_vtu13_pn.txt";
        file.open(iname, std::ios::in);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        std::vector<Scalar> valpw;
        lineNumber = 0;
        while (!file.eof()) {
            Scalar val = 0;
            file >> val;
            valpw.push_back(val);
            lineNumber++;
         }

        file.close();
        std::cout << "lineNumber: " << lineNumber << std::endl;
        std::cout << "numEle: " << this->gridView().size(0) << std::endl;

        auto& uCur = this->model().curSol();
        auto& uPrev = this->model().prevSol();
        for (const auto& element : elements(this->gridView()))
        {
            int eIdx = this->elementMapper().index(element);
            uCur[eIdx][snIdx] = 1.0 - valSw[eIdx];
            uCur[eIdx][pwIdx] = valpw[eIdx];
            uPrev[eIdx][snIdx] = 1.0 - valSw[eIdx];
            uPrev[eIdx][pwIdx] = valpw[eIdx];
        }

    }


    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }


    PrimaryVariables source(const Element &element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

//        auto globalPos = scv.center();
//        Scalar x = globalPos[0];
//        Scalar y = globalPos[1];
//        Scalar z = globalPos[2];
//
//        Scalar IX = 185.928;
//        Scalar IY = 336.804;
//
//        Scalar P1X = 0.0;
//        Scalar P1Y = 0.0;
//
//        Scalar P2X = 365.76;
//        Scalar P2Y = 0.0;
//
//        Scalar P3X = 365.76;
//        Scalar P3Y = 670.56;
//
//        Scalar P4X = 0.0;
//        Scalar P4Y = 670.56;
//
//        Scalar hX = 6.096;
//        Scalar hY = 3.048;
//        Scalar hZ = 0.6096;
//
//        Scalar lZ = this->bBoxMax()[dim-1];
//
//        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
//        fluidState.setTemperature(temperature_);
//        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
//        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);
//        Scalar densityN = FluidSystem::density(fluidState, FluidSystem::nPhaseIdx);
//        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);
//
//        auto poro = this->spatialParams().porosity(element, scv,
//                this->model().elementSolution(element, this->model().curSol()));
//
//        auto K = this->spatialParams().permeability(element, scv,
//                this->model().elementSolution(element, this->model().curSol()));
//
//        Scalar pi = 4.0*atan(1.0);
//        Scalar rw = 0.2;
//        Scalar re = 0.14*std::sqrt(hX*hX + hY*hY);
//
//        Scalar pbhI = 6.8948e7;
//        Scalar pbhP = 2.7579e7;
//
//        Scalar pn = elemVolVars[scv].pressure(nPhaseIdx);
//        Scalar pw = elemVolVars[scv].pressure(wPhaseIdx);
//        Scalar mobN = elemVolVars[scv].mobility(nPhaseIdx);
//        Scalar viscosityW = elemVolVars[scv].fluidState().viscosity(wPhaseIdx);
//
//        if(std::abs(x-IX) < hX - eps_ && std::abs(y-IY) < hY - eps_)
//        {
//            values[contiWEqIdx] = (2.0*pi*hZ*K[0][0]*densityW)/(std::log(re/rw)*viscosityW) * (pbhI - pw);
//        }
//        else if(std::abs(x-P1X) < hX - eps_ && std::abs(y-P1Y) < hY - eps_)
//        {
//            values[contiNEqIdx] = (2.0*pi*hZ*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - pn);
//        }
//        else if(std::abs(x-P2X) < hX - eps_ && std::abs(y-P2Y) < hY - eps_)
//        {
//            values[contiNEqIdx] = (2.0*pi*hZ*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - pn);
//        }
//        else if(std::abs(x-P3X) < hX - eps_ && std::abs(y-P3Y) < hY - eps_)
//        {
//            values[contiNEqIdx] = (2.0*pi*hZ*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - pn);
//        }
//        else if(std::abs(x-P4X) < hX - eps_ && std::abs(y-P4Y) < hY - eps_)
//        {
//            values[contiNEqIdx] = (2.0*pi*hZ*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - pn);
//        }
//
//        values[contiWEqIdx] /= element.geometry().volume();
//        values[contiNEqIdx] /= element.geometry().volume();

        return values;

    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        Scalar hX = 6.096;
        Scalar hY = 3.048;

        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar z = globalPos[2];

        if((std::abs(xMin_[0]-x) < hX && std::abs(xMin_[1]-y) < hY && std::abs(xMin_[2]-z) > 1.0e-3) ||
           (std::abs(xMax_[0]-x) < hX && std::abs(xMax_[1]-y) < hY && std::abs(xMax_[2]-z) > 1.0e-3))
        {
            values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0);

        Scalar pbhI = 6.8948e7;
        Scalar pbhP = 2.7579e7;

        Scalar hX = 6.096;
        Scalar hY = 3.048;

        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar z = globalPos[2];

        if(std::abs(xMin_[0]-x) < hX && std::abs(xMin_[1]-y) < hY && std::abs(xMin_[2]-z) > 1.0e-3)
        {
            values[pwIdx] = pbhI;
            values[snIdx] = 0.0;
        }
        else if(std::abs(xMax_[0]-x) < hX && std::abs(xMax_[1]-y) < hY && std::abs(xMax_[2]-z) > 1.0e-3)
        {
            values[pwIdx] = pbhP;
            values[snIdx] = 1.0;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        return values;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{


    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, 1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, 1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

//        auto fvGeometry = localView(this->model().globalFvGeometry());
//        fvGeometry.bindElement(element);
//        const auto& scv = fvGeometry.scv(this->elementMapper(element));
//        const auto globalPos = scv.center();

        Scalar depth = this->bBoxMax()[dim-1] - globalPos[dim-1];

        values[pwIdx] = 1.0e7 - densityW*this->gravity()[dim-1]*depth;
        values[snIdx] = 1.0;

        return values;
    }

    bool shouldWriteOutput() const
    {
        return this->timeManager().timeStepIndex() == 0 ||
               this->timeManager().episodeWillBeFinished() ||
               this->timeManager().willBeFinished();
    }

    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
//    void postTimeStep()
//    {
//        auto uCur = this->model().curSol();
//        int numEle = this->gridView().size(0);
//        std::vector<Scalar> valSn(numEle,0.0);
//        std::vector<Scalar> valpw(numEle,0.0);
//
//        for (const auto& element : elements(this->gridView()))
//        {
//            int eIdx = this->elementMapper().index(element);
//            int idxGlobal = this->spatialParams().getIndex_(element.geometry().center());
//            valSn[idxGlobal] = uCur[eIdx][snIdx];
//            valpw[idxGlobal] = uCur[eIdx][pwIdx];
//        }
//
//        std::ofstream file;
//        std::string outname = "spe10_tpfa_full3D";
//        outname.append(std::to_string(this->timeManager().timeStepIndex()));
//        outname.append("_sn.txt");
//        file.open(outname, std::ios::out);
//        if (file.fail())
//            throw std::ios_base::failure(std::strerror(errno));
//
//        for(int i=0; i<valSn.size(); i++)
//        {
//            file << valSn[i];
//            if(!(i == valSn.size() -1))
//                file << "\n";
//        }
//
//        file.close();
//
//
//        outname = "spe10_tpfa_full3D";
//        outname.append(std::to_string(this->timeManager().timeStepIndex()));
//        outname.append("_pw.txt");
//        file.open(outname, std::ios::out);
//        if (file.fail())
//            throw std::ios_base::failure(std::strerror(errno));
//
//        for(int i=0; i<valpw.size(); i++)
//        {
//            file << valpw[i];
//            if(!(i == valpw.size() -1))
//                file << "\n";
//        }
//
//        file.close();
//    }

        void postTimeStep()
        {
            ParentType::postTimeStep();

            auto uCur = this->model().curSol();

            std::ofstream fileInj;
            std::string outname = "spe10_tpfa_full3D";
            outname.append(std::to_string(this->timeManager().timeStepIndex()));
            outname.append("_solInj.txt");
            fileInj.open(outname, std::ios::out);
            if (fileInj.fail())
                throw std::ios_base::failure(std::strerror(errno));

            std::ofstream fileProd;
            outname = "spe10_tpfa_full3D";
            outname.append(std::to_string(this->timeManager().timeStepIndex()));
            outname.append("_solProd.txt");
            fileProd.open(outname, std::ios::out);
            if (fileProd.fail())
                throw std::ios_base::failure(std::strerror(errno));

            for (const auto& element : elements(this->gridView()))
            {
                auto globalPos = element.geometry().center();
                Scalar hX = 6.096;
                Scalar hY = 3.048;

                Scalar x = globalPos[0];
                Scalar y = globalPos[1];

                if(std::abs(xMin_[0]-x) < hX && std::abs(xMin_[1]-y) < hY)
                {
                    int eIdx = this->elementMapper().index(element);
                    Scalar Sn = 0.0;
                    Scalar pw = 0.0;
#if PROBLEM==1
                    Sn = uCur[eIdx][snIdx];
                    pw = uCur[eIdx][pwIdx];
#elif PROBLEM==2
                    const LocalFiniteElementCache feCache;
                    const auto geometryI = element.geometry();
                    Dune::GeometryType geomType = geometryI.type();

                    GlobalPosition centerI = geometryI.local(geometryI.center());
                    const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
                    std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
                    localFiniteElement.localBasis().evaluateFunction(centerI, shapeVal);

                    for (int i = 0; i < shapeVal.size(); ++i)
                    {
                        int dofIdxGlobal = this->model().dofMapper().subIndex(element, i, dim);
                        Sn += shapeVal[i]*this->model().curSol()[dofIdxGlobal][snIdx];
                        pw += shapeVal[i]*this->model().curSol()[dofIdxGlobal][pwIdx];
                    }
#endif
                    fileInj << Sn <<  "\n";
                    fileInj << pw <<  "\n";
                }
                else if(std::abs(xMax_[0]-x) < hX && std::abs(xMax_[1]-y) < hY)
                {
                    int eIdx = this->elementMapper().index(element);
                    Scalar Sn = 0.0;
                    Scalar pw = 0.0;
#if PROBLEM==1
                    Sn = uCur[eIdx][snIdx];
                    pw = uCur[eIdx][pwIdx];
#elif PROBLEM==2
                    const LocalFiniteElementCache feCache;
                    const auto geometryI = element.geometry();
                    Dune::GeometryType geomType = geometryI.type();

                    GlobalPosition centerI = geometryI.local(geometryI.center());
                    const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
                    std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
                    localFiniteElement.localBasis().evaluateFunction(centerI, shapeVal);

                    for (int i = 0; i < shapeVal.size(); ++i)
                    {
                        int dofIdxGlobal = this->model().dofMapper().subIndex(element, i, dim);
                        Sn += shapeVal[i]*this->model().curSol()[dofIdxGlobal][snIdx];
                        pw += shapeVal[i]*this->model().curSol()[dofIdxGlobal][pwIdx];
                    }
#endif
                    fileInj << Sn <<  "\n";
                    fileInj << pw <<  "\n";
                }
            }


            fileInj.close();
            fileProd.close();
        }

    // \}

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& eIdx = outputModule.createScalarField("eIdx", 0);
#if PROBLEM==2
        auto& SnCell = outputModule.createScalarField("SnCell", 0);
#endif
        for (const auto& element : elements(this->gridView()))
        {
            int idx = this->elementMapper().index(element);
            eIdx[idx] = idx;
#if PROBLEM==2
            const LocalFiniteElementCache feCache;
            const auto geometryI = element.geometry();
            Dune::GeometryType geomType = geometryI.type();

            GlobalPosition centerI = geometryI.local(geometryI.center());
            const LocalFiniteElement &localFiniteElement = feCache.get(geomType);
            std::vector<Dune::FieldVector<Scalar, 1> > shapeVal;
            localFiniteElement.localBasis().evaluateFunction(centerI, shapeVal);

            Scalar SnI = 0.0;
            for (int i = 0; i < shapeVal.size(); ++i)
            {
                int dofIdxGlobal = this->model().dofMapper().subIndex(element, i, dim);
                SnI += shapeVal[i]*this->model().curSol()[dofIdxGlobal][snIdx];
            }
            SnCell[idx] = SnI;
#endif
        }
    }

private:
    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
    std::string name_;

    Scalar pIn_;
    Scalar pOut_;
    Scalar episodeLength_;
    Scalar IPoreVol_;
    Scalar P1PoreVol_;
    Scalar P2PoreVol_;
    Scalar P3PoreVol_;
    Scalar P4PoreVol_;
    GlobalPosition xMax_;
    GlobalPosition xMin_;
};
} //end namespace

#endif