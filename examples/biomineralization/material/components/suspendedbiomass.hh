/*
 * biosusp.hh
 *
 *  Created on: 20.04.2011
 *      Author: hommel
 */

/*!
 * \file
 *
 * \brief A class for the suspended biomass component properties
 */
#ifndef DUMUX_BIOSUSP_HH
#define DUMUX_BIOSUSP_HH

#include <dumux/material/components/base.hh>

namespace Dumux {
namespace Components {

/*!
 * \brief A class for the suspended biomass component properties
 */
template <class Scalar>
class SuspBiomass
: public Components::Base<Scalar, SuspBiomass<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the Biosusp.
     */
    static std::string name()
    { return "Suspended_Biomass"; }

    /*!
     * \brief The mass in [kg] of one mole of Biosusp.
     */
    static Scalar molarMass()  // TODO what is the molar Mass of suspended biomass???
    { return 1; } // kg/mol
        //based on a cell mass of 2.5e-16, the molar mass of cells would be 1.5e8 kg/mol.
};

} // end namespace Components
} // end namespace Dumux

#endif