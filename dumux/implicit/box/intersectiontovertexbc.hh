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
 * \brief Convert intersection boundary types to vertex boundary types
 */
#ifndef DUMUX_INTERSECTIONTOVERTEXBC_HH
#define DUMUX_INTERSECTIONTOVERTEXBC_HH

#include <dumux/implicit/box/properties.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief Convert intersection boundary types to vertex boundary types
 */
template<typename TypeTag>
class IntersectionToVertexBC
{
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,
    };

    typedef typename GridView::template Codim<dim>::Entity Vertex;

    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

public:
    IntersectionToVertexBC(const Problem& problem)
    : problem_(problem)
    {
        vertexBC.resize(problem.vertexMapper().size());
        for (int vIdx = 0; vIdx < vertexBC.size(); vIdx++)
            vertexBC[vIdx].setAllNeumann();

        for (const auto& element : elements(problem.gridView())) {
            Dune::GeometryType geomType = element.geometry().type();
            const ReferenceElement &refElement = ReferenceElements::general(geomType);

            for (const auto& intersection : intersections(problem.gridView(), element)) {
                if (!intersection.boundary())
                    continue;

                BoundaryTypes bcTypes;
                problem.boundaryTypes(bcTypes, intersection);

                if (!bcTypes.hasDirichlet())
                    continue;

                int fIdx = intersection.indexInInside();
                int numFaceVerts = refElement.size(fIdx, 1, dim);
                for (int faceVertexIdx = 0; faceVertexIdx < numFaceVerts; ++faceVertexIdx)
                {
                    int vIdx = refElement.subEntity(fIdx, 1, faceVertexIdx, dim);
                    int vIdxGlobal = problem.vertexMapper().subIndex(element, vIdx, dim);
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        if (bcTypes.isDirichlet(eqIdx))
                            vertexBC[vIdxGlobal].setDirichlet(eqIdx);
                }
            }
        }
    }

    void boundaryTypes(BoundaryTypes& values, const Vertex& vertex) const
    {
        values.setAllNeumann();
        int vIdxGlobal = problem_.vertexMapper().index(vertex);

        const BoundaryTypes& bcTypes = vertexBC[vIdxGlobal];

        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            if (bcTypes.isDirichlet(eqIdx))
                values.setDirichlet(eqIdx);
    }

private:
    const Problem& problem_;
    std::vector<BoundaryTypes> vertexBC;
};
}

#endif
