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
 *
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
#define DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(ThermalConductivityModel);
NEW_PROP_TAG(CompositionalFluidState);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(Scalar);
}

/*!
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
template<class TypeTag>
class PlotThermalConductivityModel
{
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef CompositionalFluidState<Scalar, FluidSystem> FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

public:
    /*!
     * \brief Constructor
     *
     * Initializes the fluid system.
     *
     * \param temperature temperature in \f$\mathrm{[K]}\f$
     * \param pressure reference pressure in \f$\mathrm{[Pa]}\f$
     */
    PlotThermalConductivityModel(Scalar temperature = 283.15,
                                 Scalar pressure = 1e5)
    : numIntervals_(1000)
    {
        FluidState fluidstate;
        fluidstate.setTemperature(temperature);
        fluidstate.setPressure(wPhaseIdx, pressure);
        fluidstate.setPressure(nPhaseIdx, pressure);
        lambdaW_ = FluidSystem::template thermalConductivity<FluidState>(fluidstate, wPhaseIdx);
        lambdaN_ = FluidSystem::template thermalConductivity<FluidState>(fluidstate, nPhaseIdx);
    }

    /*!
     * \brief Add a effective thermal conductivity-saturation curve to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param params The material law parameters
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addlambdaeffcurve(GnuplotInterface<Scalar> &gnuplot,
                           Scalar porosity,
                           Scalar rhoSolid,
                           Scalar lambdaSolid,
                           Scalar lowerSat = 0.0,
                           Scalar upperSat = 1.0,
                           std::string curveName = "lambdaeff",
                           std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> lambda(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            lambda[i] = ThermalConductivityModel::effectiveThermalConductivity(sw[i], lambdaW_,
                                                                               lambdaN_, lambdaSolid,
                                                                               porosity, rhoSolid);
        }

        gnuplot.setXlabel("wetting phase saturation [-]");
        gnuplot.setYlabel("effective thermal conductivity [W/(m K)]");
        gnuplot.addDataSetToPlot(sw, lambda, curveName, curveOptions);
    }

private:
    int numIntervals_;
    Scalar lambdaN_;
    Scalar lambdaW_;
};

} // end namespace Dumux

#endif // DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
