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
 * \ingroup PoroElastic
 * \brief Adds vtk output fields specific to the poro-elastic model
 */
#ifndef DUMUX_POROELASTIC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_POROELASTIC_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup PoroElastic
 * \brief Adds vtk output fields specific to the poro-elastic model
 */
class PoroElasticVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.displacement(); }, "u");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.porosity(); }, "porosity");
    }
};

} // end namespace Dumux

#endif