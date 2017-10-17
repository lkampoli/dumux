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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief The spatial parameters for the 2pDFM Problem which uses the
 *        twophase discrete fracture model.
 */
#ifndef DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH
#define DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/2pdfm/implicit/model.hh>
#include <dumux/io/artgridcreator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/spatialparams/implicit.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPDFMSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPDFMSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPDFMSpatialParams, SpatialParams, TwoPDFMSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPDFMSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the 2PDFMProblem which uses the
 *        twophase model
 */
template<class TypeTag>
class TwoPDFMSpatialParams : public ImplicitSpatialParams<TypeTag>
{

    template<int dim>
    struct FaceLayout
    {
        bool contains (Dune::GeometryType geomType)
        {
            return geomType.dim() == dim - 1;
        }
    };

    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, FaceLayout> FaceMapper;

    enum {
        dim = GridView::dimension
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    TwoPDFMSpatialParams(const GridView& gridView)
        : ParentType(gridView), gridView_(gridView),
        faceMapper_(gridView), vertexMapper_(gridView),
        fractureMapper_(gridView)
    {
        inactivateFractures_ = false;

        Scalar mD = 1e-12 * 1e-3; //miliDarcy

        swrf_    = 0.00;
        swrm_    = 0.00;
        SnrF_    = 0.00;
        SnrM_    = 0.00;
        pdf_     = 1000; //2.5*1e4;
        pdm_     = 2000; //2.5*1e4;
        lambdaF_ = 2.0;
        lambdaM_ = 2.0;

        rockMatrixMaterialParams_.setSwr(swrm_);
        rockMatrixMaterialParams_.setSnr(SnrM_);
        fractureMaterialParams_.setSwr(swrf_);
        fractureMaterialParams_.setSnr(SnrF_);

        rockMatrixMaterialParams_.setPe(pdm_);
        rockMatrixMaterialParams_.setLambda(lambdaM_);
        fractureMaterialParams_.setPe(pdf_);
        fractureMaterialParams_.setLambda(lambdaF_);

        KMatrix_   = 1 * mD; //m^2
        KFracture_ = 1e5 * mD; //m^2

        porosityMatrix_   =    0.25;
        porosityFracture_ = 0.10;
        fractureWidth_    = 1e-2;

        fractureMapper_.map();
    }

    /*!
     * \brief Intrinsic permability
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Intrinsic permeability
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 int scvIdx) const
    {
        return KMatrix_;
    }

    /*!
     * \brief Intrinsic permeability of fractures.
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     */
    Scalar intrinsicPermeabilityFracture(const Element &element,
                                         const FVElementGeometry &fvGeometry,
                                         int scvIdx) const
    {
        return KFracture_;
    }
    /*!
     * \brief Porosity
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    { return porosityMatrix_; }

    /*!
     * \brief Porosity Fracture
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity Fracture
     */
    Scalar porosityFracture(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        return porosityFracture_;
    }
    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
    {
        return rockMatrixMaterialParams_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object of the Fracture
     */
    const MaterialLawParams& materialLawParamsFracture(const Element &element,
                                                    const FVElementGeometry &fvGeometry,
                                                    int scvIdx) const
    {
        DUNE_UNUSED int vIdxGlobal = vertexMapper_.subIndex(element, scvIdx, dim);

        // be picky if called for non-fracture vertices
        assert(isVertexFracture(vIdxGlobal));

        return fractureMaterialParams_;
    }

    /*!
     * \brief Checks whether vertex is a fracture.
     *
     * \param element The current element
     * \param localVertexIdx Vertex index to be checked
     */
    bool isVertexFracture(const Element &element, int localVertexIdx) const
    {
        if (inactivateFractures_)
        {
            return false;
        }
        int vIdxGlobal = vertexMapper_.subIndex(element, localVertexIdx, dim);
        return fractureMapper_.isDuneFractureVertex(vIdxGlobal);
    }

    /*!
     * \brief Checks whether vertex is a fracture.
     *
     * \param vIdxGlobal Vertex index to be checked
     */
    bool isVertexFracture(int vIdxGlobal) const
    {
        if (inactivateFractures_)
        {
            return false;
        }
        return fractureMapper_.isDuneFractureVertex(vIdxGlobal);
    }

    /*!
     * \brief Checks whether element edge is a fracture.
     *
     * \param element The current element
     * \param localFaceIdx Face index to be checked
     */
    bool isEdgeFracture(const Element &element, int localFaceIdx) const
    {
        int fIdxGlobal = faceMapper_.subIndex(element, localFaceIdx, 1);
        return fractureMapper_.isDuneFractureEdge(fIdxGlobal);
    }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \param globalFaceIdx Global face index of which the width is returned
     */
    Scalar fractureWidth(int globalFaceIdx) const
    {
        return fractureWidth_;
    }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \param element The current element
     * \param localFaceIdx Local face index of which the width is returned
     */
    Scalar fractureWidth(const Element &element, int localFaceIdx) const
    {
        return fractureWidth_;
    }

    Scalar swrf_;
    Scalar swrm_;
    Scalar SnrF_;
    Scalar SnrM_;
    Scalar lambdaF_;
    Scalar lambdaM_;
    Scalar pdf_;
    Scalar pdm_;

private:
    Scalar KMatrix_;
    Scalar KFracture_;
    Scalar porosityMatrix_;
    Scalar porosityFracture_;

    Scalar fractureWidth_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams rockMatrixMaterialParams_;
    bool inactivateFractures_;

    const GridView gridView_;
    const FaceMapper faceMapper_;
    const VertexMapper vertexMapper_;

    FractureMapper<TypeTag> fractureMapper_;
};

} // end namespace
#endif // DUMUX_TEST_2PDFM_SPATIAL_PARAMETERS_HH