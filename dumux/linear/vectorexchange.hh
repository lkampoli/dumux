// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Linear
 * \brief Contains a class to exchange entries of a vector
 */
#ifndef DUMUX_VECTOR_EXCHANGE_HH
#define DUMUX_VECTOR_EXCHANGE_HH

#warning "This header is deprecated and will be removed after release 3.2. Use linearsolver/linearsolvertraits.hh"

#include <dumux/parallel/vectorcommdatahandle.hh>

namespace Dumux {
    template<class Mapper, class Vector>
    using VectorExchange [[deprecated("Use VectorCommDataHandleEqual<Mapper, Vector, 0> instead. Will be removed after 3.2!")]] = VectorCommDataHandleEqual<Mapper, Vector, 0/*elementCodim*/>;
} // end namespace Dumux

#endif
