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
 * \brief Tests the grid creator class for models using facet coupling.
 */
#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/io/vtkfunction.hh>
#include <dumux/multidomain/facet/gridcreator.hh>
#include <dumux/multidomain/facet/gridvertexadapter.hh>
#include <dumux/multidomain/facet/box/vertexmapper.hh>

#ifndef BULKGRIDTYPE // default to alu grid if not provided by CMake
#define BULKGRIDTYPE Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>
#endif

#ifndef FACETGRIDTYPE // default to foam grid if not provided by CMake
#define FACETGRIDTYPE Dune::FoamGrid<2, 3>
#endif

// returns a ficticious displacement for a given position
template<class GlobalPosition>
GlobalPosition getDisplacement(const GlobalPosition& pos)
{
    static const auto dim = GlobalPosition::size();

    GlobalPosition d(1.0);

    if (pos[0] < 0.5 && pos[1] < 0.5)
    {
        d *= -1.0;
        if (dim == 3 && pos[2] > 0.5)
            d[2] = 1.0;
        d /= d.two_norm();
        return d;
    }
    else if (pos[0] > 0.5 && pos[1] < 0.5)
    {
        d *= -1.0;
        d[0] = 1.0;
        if (dim == 3 && pos[2] > 0.5)
        {
            d[0] = 1.0;
            d[2] = 1.0;
        }
        d /= d.two_norm();
        return d;
    }
    else if (pos[0] > 0.5 && pos[1] > 0.5)
    {
        if (dim == 3 && pos[2] < 0.5)
            d[2] = -1.0;
        d /= d.two_norm();
        return d;
    }
    else if (pos[0] < 0.5 && pos[1] > 0.5)
    {
        d *= -1.0;
        d[1] = 1.0;
        if (dim == 3 && pos[2] > 0.5)
            d[2] = 1.0;
        d /= d.two_norm();
        return d;
    }
    else
        DUNE_THROW(Dune::InvalidStateException, "Invalid position provided: " << pos);
}

//! main function
int main (int argc, char *argv[]) try
{
    // initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // we must have one single command line argument
    if (argc != 2)
        DUNE_THROW(Dune::IOError, "You must provide one command line argument, which is the grid file");

    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = FACETGRIDTYPE;
    using GridCreator = Dumux::FacetCouplingGridCreator<BulkGrid, FacetGrid>;
    GridCreator gridCreator;
    gridCreator.makeGrids(argv[1]);

    // dimension of the bulk grid
    static constexpr int bulkDim = BulkGrid::dimension;

    // check grid sizes
    const auto& bulkGridView = gridCreator.grid<0>().leafGridView();
    const auto& facetGridView = gridCreator.grid<1>().leafGridView();

    // data converter
    Dumux::FacetGridVertexAdapter<GridCreator, 0, 1> facetGridVertexAdapter(gridCreator);

    // enriched vertex dof mapper
    Dumux::EnrichedVertexDofMapper<typename BulkGrid::LeafGridView> mapper(bulkGridView);

    // print out number of dofs managed by the mapper before enrichment
    const auto numDofsBefore = mapper.size();
    std::cout << "Number of dofs before enrichment: " << numDofsBefore << std::endl;
    if (bulkDim == 3 && numDofsBefore != 153)
        DUNE_THROW(Dune::InvalidStateException, "3d test: Number of dofs before enrichment is expected to be 153 but is " << numDofsBefore);
    if (bulkDim == 2 && numDofsBefore != 169)
        DUNE_THROW(Dune::InvalidStateException, "2d test: Number of dofs before enrichment is expected to be 169 but is " << numDofsBefore);

    // enrich nodes subject to facet grid
    mapper.enrich(facetGridView, facetGridVertexAdapter, true);
    const auto numDofsAfter = mapper.size();
    std::cout << "Number of dofs after enrichment: " << numDofsAfter << std::endl;
    if (bulkDim == 3 && numDofsAfter != 164)
        DUNE_THROW(Dune::InvalidStateException, "3d test: Number of dofs after enrichment is expected to be 164 but is " << numDofsAfter);
    if (bulkDim == 2 && numDofsAfter != 180)
        DUNE_THROW(Dune::InvalidStateException, "2d test: Number of dofs after enrichment is expected to be 180 but is " << numDofsAfter);

    // vector for output of ficticious displacement
    static constexpr int bulkDimWorld = BulkGrid::dimensionworld;
    using Displacement = Dune::FieldVector<double, bulkDimWorld>;
    std::vector< Displacement > displacement(mapper.size(), Displacement(0.0));
    for (const auto& e : elements(bulkGridView))
    {
        for (int i = 0; i < e.geometry().corners(); ++i)
        {
            if ( mapper.isEnriched(e.template subEntity<bulkDim>(i)) )
                displacement[ mapper.subIndex(e, i, BulkGrid::dimension) ] = getDisplacement(e.geometry().center());
            else
                displacement[ mapper.subIndex(e, i, BulkGrid::dimension) ] = 0.0;
        }
    }

    // nonconforming vectorial vtk function
    using VTKFunction = Dumux::Vtk::VectorP1NonConformingVTKFunction< typename BulkGrid::LeafGridView,
                                                                      Dumux::EnrichedVertexDofMapper<typename BulkGrid::LeafGridView>,
                                                                      std::vector< Displacement > >;
    auto displacementFunction = std::make_shared< VTKFunction >(bulkGridView, mapper, displacement, "displacement", BulkGrid::dimensionworld);

    // write .vtk file for the bulk grid
    using BulkWriter = Dune::VTKWriter<typename BulkGrid::LeafGridView>;
    BulkWriter bulkWriter(bulkGridView, Dune::VTK::nonconforming);
    bulkWriter.addVertexData(displacementFunction);
    bulkWriter.write("displacement_" + std::to_string(BulkGrid::dimension) + "d");

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e)
{
    std::cout << e << std::endl;
    return 1;
}
catch (...)
{
    std::cout << "Unknown exception thrown" << std::endl;
    return 1;
}
