/*
 * biosub.hh
 *
 *  Created on: 20.04.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the Substrate component properties
 */
#ifndef DUMUX_BIOSUB_HH
#define DUMUX_BIOSUB_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {
/*!
 * \brief A class for the Substrate component properties
 */
template <class Scalar>
class Substrate
: public Components::Base<Scalar, Substrate<Scalar> >
{
public:
    /*!
     * \brief A human readable name for Substrate (bacteria food).
     */
    static std::string name()
    { return "Substrate"; }

    /*!
     * \brief The mass in [kg] of one mole of Substrate.
     */
    static Scalar molarMass()   //Molar mass of Glucose is assumed.
    { return 0.18016; } // kg/mol
};

} // end namespace Components
} // end namespace Dumux

#endif