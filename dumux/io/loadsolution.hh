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
 * \ingroup InputOutput
 * \brief read from a file into a solution vector
 */
#ifndef DUMUX_IO_LOADSOLUTION_HH
#define DUMUX_IO_LOADSOLUTION_HH

#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/io/vtk/vtkreader.hh>

namespace Dumux {

namespace Detail {
//! helper struct detecting if a PrimaryVariables object has a state() function
struct hasState
{
    template<class PrimaryVariables>
    auto operator()(PrimaryVariables&& priVars)
    -> decltype(priVars.state())
    {}
};

} // end namespace Detail

/*!
 * \ingroup InputOutput
 * \brief helper function to read from a file into a solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadSolutionFromVtkFile(const std::string fileName,
                             const VTKReader::DataType& dataType,
                             PvNamesFunc&& pvNamesFunc,
                             SolutionVector& sol)
-> typename std::enable_if_t<!decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);
    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;

    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        const auto pvName = pvNamesFunc(pvIdx);
        const auto vec = vtu.readData<std::vector<Scalar>>(pvName, dataType);
        if (vec.size() != sol.size())
            DUNE_THROW(Dune::IOError, "Size mismatch between solution vector and read data (" << sol.size() << " != " << vec.size() << ")");

        for (std::size_t i = 0; i < sol.size(); ++i)
            sol[i][pvIdx] = vec[i];
    }
}

/*!
 * \ingroup InputOutput
 * \brief helper function to read from a file into a solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadSolutionFromVtkFile(const std::string fileName,
                             const VTKReader::DataType& dataType,
                             PvNamesFunc&& pvNamesFunc,
                             SolutionVector& sol)
-> typename std::enable_if_t<decltype(isValid(Detail::hasState())(sol[0]))::value, void>
{
    VTKReader vtu(fileName);
    const auto vec = vtu.readData<std::vector<int>>("phase presence", dataType);
    std::unordered_set<int> states;
    for (size_t i = 0; i < sol.size(); ++i)
    {
        const int state = vec[i];
        sol[i].setState(state);
        states.insert(state);
    }

    using PrimaryVariables = typename SolutionVector::block_type;
    using Scalar = typename PrimaryVariables::field_type;
    for (size_t pvIdx = 0; pvIdx < PrimaryVariables::dimension; ++pvIdx)
    {
        if (pvNamesFunc(pvIdx, 1) == pvNamesFunc(pvIdx, 2))
        {
            const auto vec = vtu.readData<std::vector<Scalar>>(pvNamesFunc(pvIdx, 1), dataType);
            for (size_t i = 0; i < sol.size(); ++i)
                sol[i][pvIdx] = vec[i];
        }
        else
        {
            std::unordered_map<int, std::vector<Scalar>> switchedPvsSol;
            for (const auto& state : states)
                switchedPvsSol[state] = vtu.readData<std::vector<Scalar>>(pvNamesFunc(pvIdx, state), dataType);

            for (size_t i = 0; i < sol.size(); ++i)
                sol[i][pvIdx] = switchedPvsSol[sol[i].state()][i];
        }
    }
}

/*!
 * \ingroup InputOutput
 * \brief helper function to read from a file into a solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadSolutionFromVtkFile(const std::string fileName,
                             const VTKReader::DataType& dataType,
                             PvNamesFunc&& pvNamesFunc,
                             SolutionVector& sol)
-> typename std::enable_if_t<decltype(isMultiTypeBlockVector<SolutionVector>())::value, void>
{}

/*!
 * \ingroup InputOutput
 * \brief helper function to read from two files into a staggered solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadStaggeredSolutionFromVtkFiles(const std::string baseFileName,
                                       PvNamesFunc&& pvNamesFunc,
                                       SolutionVector& sol)
-> typename std::enable_if_t<!decltype(isMultiTypeBlockVector<SolutionVector>())::value, void>
{}

/*!
 * \ingroup InputOutput
 * \brief helper function to read from two files into a staggered solution vector
 */
template <class SolutionVector, class PvNamesFunc>
auto loadStaggeredSolutionFromVtkFiles(const std::string baseFileName,
                                       PvNamesFunc&& pvNamesFunc,
                                       SolutionVector& sol)
-> typename std::enable_if_t<decltype(isMultiTypeBlockVector<SolutionVector>())::value, void>
{

    // assume that the first component contains the cell data
    auto& cellSol = sol[Dune::index_constant<0>{}];
    using CellPrimaryVariables = typename std::decay_t<decltype(cellSol)>::block_type;
    using Scalar = typename CellPrimaryVariables::field_type;
    VTKReader cellVtk(baseFileName + ".vtu");
    for (size_t pvIdx = 0; pvIdx < CellPrimaryVariables::dimension; ++pvIdx)
    {
        const auto vec = cellVtk.readData<std::vector<Scalar>>(pvNamesFunc(pvIdx),
                                                               VTKReader::DataType::cellData);
        for (size_t i = 0; i < cellSol.size(); ++i)
            cellSol[i][pvIdx] = vec[i];
    }

    // assume that the second component contains the face data
    auto& faceSol = sol[Dune::index_constant<1>{}];
    using FacePrimaryVariables = typename std::decay_t<decltype(faceSol)>::block_type;
    auto nameSize = baseFileName.size();
    // assume that baseFileName contains numbers like '-000123' in the end
    VTKReader faceVtk(baseFileName.substr(0, nameSize - 6) + "-face" + baseFileName.substr(nameSize - 6) + ".vtp");
    const auto vec = faceVtk.readData<std::vector<Scalar>>(pvNamesFunc(CellPrimaryVariables::dimension),
                                                           VTKReader::DataType::pointData);
    for (size_t i = 0; i < faceSol.size(); ++i)
        faceSol[i][0] = std::accumulate(&vec[3*i], &vec[3*(i+1)], 0.0);
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primray variable names of a model with privar state
 */
template<class ModelTraits, class FluidSystem>
std::string primaryVariableName(int pvIdx, int state)
{
    static auto numStates = (1 << ModelTraits::numPhases()) - 1;
    const auto paramNameWithState = "LoadSolution.PriVarNamesState" + std::to_string(state);

    if (hasParam("LoadSolution.PriVarNames") && !hasParam(paramNameWithState))
    {
        DUNE_THROW(Dune::NotImplemented, "please provide LoadSolution.PriVarNamesState1..." << numStates
                   << " or remove LoadSolution.PriVarNames to use default names");
    }
    else if (hasParam(paramNameWithState))
    {
        const auto pvNames = getParam<std::vector<std::string>>(paramNameWithState);
        return pvNames[pvIdx];
    }
    else
        return ModelTraits::template primaryVariableName<FluidSystem>(pvIdx, state);
}

/*!
 * \ingroup InputOutput
 * \brief helper function to determine the primray variable names of a model
 */
template<class ModelTraits>
std::string primaryVariableName(int pvIdx)
{
    if (hasParam("LoadSolution.PriVarNames"))
    {
        static auto pvNames = getParam<std::vector<std::string>>("LoadSolution.PriVarNames");
        return pvNames[pvIdx];
    }
    else
        return ModelTraits::primaryVariableName(pvIdx);
}


/*!
 * \ingroup InputOutput
 * \brief load a solution vector from file
 * \note Supports the following file extensions: *.vtu *.vtp
 */
template <class SolutionVector, class PvNamesFunc>
void loadSolution(const std::string& fileName,
                  DiscretizationMethod discMethod,
                  PvNamesFunc&& pvNamesFunc,
                  SolutionVector& sol)
{
    const auto extension = fileName.substr(fileName.find_last_of(".") + 1);

    if (extension == "vtu" || extension == "vtp")
    {
        const auto dataType = discMethod == DiscretizationMethod::box
                              ? VTKReader::DataType::pointData : VTKReader::DataType::cellData;
        loadSolutionFromVtkFile(fileName, dataType, pvNamesFunc, sol);
    }
    else if (extension == fileName && discMethod == DiscretizationMethod::staggered)
        loadStaggeredSolutionFromVtkFiles(fileName, pvNamesFunc, sol);
    else
        DUNE_THROW(Dune::NotImplemented, "loadSolution for extension " << extension);
}

} // namespace Dumux

#endif
