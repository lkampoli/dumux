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
#ifndef DUMUX_ICP_SPATIAL_PARAMS_HH
#define DUMUX_ICP_SPATIAL_PARAMS_HH

// We include the basic spatial parameters for finite volumes file from which we will inherit
#include <dumux/material/spatialparams/fv.hh>
// We include the files for the two-phase laws
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
// We include the laws for changing porosity and permeability
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
namespace Dumux{

//In the ICPSpatialParams class we define all functions needed to describe the spatial distributed parameters.
template<class TypeTag>
class ICPSpatialParams
: public FVSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                         GetPropType<TypeTag, Properties::Scalar>,
                         ICPSpatialParams<TypeTag>>
{
    // We introduce using declarations that are derived from the property system which we need in this class
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ICPSpatialParams<TypeTag>>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using CoordScalar = typename GridView::ctype;
    enum { dimWorld=GridView::dimensionworld };

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // In the constructor be read some values from the `params.input` and initialize the material params laws
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

    ICPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {

        try
        {
        referencePorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity", 0.4);
        referencePermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability", 2.e-10);

        }
        catch (Dumux::ParameterException &e)
        {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }

        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(0.05);

        // parameters for the Brooks-Corey law
        materialParams_.setPe(1e4);
        materialParams_.setLambda(2.0);
    }

    // We return the reference or initial porosity
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }


    // We return the volume fraction of the inert component
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 1.0-referencePorosity_; }

    // We return the updated porosity
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol, referencePorosity_); }

    // We return the updated permeability
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        const auto poro = porosity(element, scv, elemSol);
        return permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, poro);
    }

    // Return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialParams_; }

    // Define the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    MaterialLawParams materialParams_;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PorosityLaw = PorosityPrecipitation<Scalar, ModelTraits::numFluidComponents(), SolidSystem::numComponents - SolidSystem::numInertComponents>;
    PorosityLaw poroLaw_;
    PermeabilityKozenyCarman<PermeabilityType> permLaw_;
    Scalar referencePorosity_;
    PermeabilityType referencePermeability_ = 0.0;
};

}

#endif
