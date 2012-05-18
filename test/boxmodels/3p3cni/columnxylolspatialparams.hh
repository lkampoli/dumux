// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-     by Holger Class                                 *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Definition of the spatial parameters for the column problem.
 */
#ifndef DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH
#define DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH

#include <dumux/boxmodels/3p3c/3p3cindices.hh>
#include <dumux/material/spatialparams/boxspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3pparams.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ColumnSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ColumnSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ColumnSpatialParams, SpatialParams, Dumux::ColumnSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ColumnSpatialParams, MaterialLaw)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef ParkerVanGen3P<Scalar> type;
};
}

/*!
 * \ingroup ThreePThreeCNIModel
 *
 * \brief Definition of the spatial parameters for the column problem
 */
template<class TypeTag>
class ColumnSpatialParams : public BoxSpatialParams<TypeTag>
{
    typedef BoxSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> DimVector;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;


public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    ColumnSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        // intrinsic permeabilities
        fineK_ = 1.4e-11;
        coarseK_ = 1.4e-8;

        // porosities
        finePorosity_ = 0.46;
        coarsePorosity_ = 0.46;

        // specific heat capacities
        fineHeatCap_ = 850.;
        coarseHeatCap_ = 84000.;

        // residual saturations
        fineMaterialParams_.setSwr(0.12);
        fineMaterialParams_.setSwrx(0.12);
        fineMaterialParams_.setSnr(0.10);
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.12);
        coarseMaterialParams_.setSwrx(0.12);
        coarseMaterialParams_.setSnr(0.10);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.5);
        fineMaterialParams_.setVgN(4.0);
        coarseMaterialParams_.setVgN(4.0);

        coarseMaterialParams_.setkrRegardsSnr(true);
        fineMaterialParams_.setkrRegardsSnr(true);

        // parameters for adsorption
        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);
    }

    ~ColumnSpatialParams()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution The global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    };

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeomtry,
                    const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeomtry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }


    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeomtry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeomtry,
                                               const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeomtry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeomtry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    double heatCapacity(const Element &element,
                        const FVElementGeometry &fvGeomtry,
                        const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeomtry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineHeatCap_ * 2650 // density of sand [kg/m^3]
                * (1 - porosity(element, fvGeomtry, scvIdx));
        else
            return coarseHeatCap_ * 2650 // density of sand [kg/m^3]
                * (1 - porosity(element, fvGeomtry, scvIdx));
    }

    /*!
     * \brief Calculate the heat flux \f$[W/m^2]\f$ through the
     *        rock matrix based on the temperature gradient \f$[K / m]\f$
     *
     * This is only required for non-isothermal models.
     *
     * \param heatFlux The resulting heat flux vector
     * \param fluxDat The flux variables
     * \param elemVolVars The volume variables
     * \param tempGrad The temperature gradient
     * \param element The current finite element
     * \param fvGeomtry The finite volume geometry of the current element
     * \param faceIdx The local index of the sub-control volume face where
     *                    the matrix heat flux should be calculated
     */
    void matrixHeatFlux(DimVector &heatFlux,
                        const FluxVariables &fluxDat,
                        const ElementVolumeVariables &elemVolVars,
                        const DimVector &tempGrad,
                        const Element &element,
                        const FVElementGeometry &fvGeomtry,
                        int faceIdx) const
    {
        static const Scalar ldry = 0.35;
        static const Scalar lSw1 = 1.8;
        static const Scalar lSn1 = 0.65;

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fvGeomtry.subContVolFace[faceIdx].i;
        const int j = fvGeomtry.subContVolFace[faceIdx].j;
        Scalar Sw = std::max(0.0, (elemVolVars[i].saturation(wPhaseIdx) +
                                   elemVolVars[j].saturation(wPhaseIdx)) / 2);
        Scalar Sn = std::max(0.0, (elemVolVars[i].saturation(nPhaseIdx) +
                                   elemVolVars[j].saturation(nPhaseIdx)) / 2);

        // the heat conductivity of the matrix. in general this is a
        // tensorial value, but we assume isotropic heat conductivity.
        Scalar heatCond = ldry + sqrt(Sw) * (lSw1-ldry) + sqrt(Sn) * (lSn1-ldry);

        // the matrix heat flux is the negative temperature gradient
        // times the heat conductivity.
        heatFlux = tempGrad;
        heatFlux *= -heatCond;
    }

private:
    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        if (0.90 <= pos[1])
            return true;
        else return false;
    };

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar fineHeatCap_;
    Scalar coarseHeatCap_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
};

}

#endif
