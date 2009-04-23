// $Id: brineco2problem.hh 609 2008-09-18 17:01:25Z melanie $
#ifndef DUNE_BRINECO2PROBLEM_HH
#define DUNE_BRINECO2PROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

//#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/minc2p2cni/minc2p2cniproblem.hh>
#include "mincco2_soilproperties.hh"


/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch
 */

namespace Dune
{
//! base class that defines the parameters of a diffusion equation
/*! An interface for defining parameters for the stationary diffusion equation
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 *
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar, int numEq>
class BrineCO2Problem : public TwoPTwoCNIProblem<Grid, Scalar, numEq> {
    typedef typename Grid::ctype DT;
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    enum {dim=Grid::dimension};
    enum {numCont = numEq/3};
    //    enum {wPhase = 0, nPhase = 1};
    //    enum {pWIdx = 0, switchIdx = 1, teIdx = 2};


public:

    virtual FieldVector<Scalar,numEq> q (const FieldVector<DT,dim>& x, const Entity& e,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<DT,dim>& x, const Entity& e,
                                                                  const IntersectionIterator& intersectionIt,
                                                                  const FieldVector<DT,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::neumann);

        if(x[0] >= 10-1e-2 )
            values = Dune::BoundaryConditions::dirichlet;
        if(x[0] < 1e-2 && x[1] < 5)
            values[2] = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual void dirichletIndex(const FieldVector<DT,dim>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<DT,dim>& xi, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual FieldVector<Scalar,numEq> g (const FieldVector<DT,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0);
        for (int nC=0; nC < numEq; nC+=(numEq/numCont)){

            if (nC<3){
                values[0]=1.013e5 + (depthBOR_ - x[1])*1045*9.81;
                values[1]=0.0;
                values[2]=283.15 + (depthBOR_ - x[1])*0.03;
            }
            else {
                values[nC]=values[0]*0.99;
                values[nC+1]=0.0;
                values[nC+2]=values[2]*0.99;
            }
        }

        return values;
    }

    virtual FieldVector<Scalar,numEq> J (const FieldVector<DT,dim>& x, const Entity& e,
                                         const IntersectionIterator& intersectionIt, const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0.0);
        if(x[0] <= 1.e-2 && x[1] <= 5.)
        {
            values[0] = 0.0;//-4.046e-5;
            values[1] = -0.2;
        }
        return values;
    }

    // Initial Conditions for global vector x, element e and local vector xi
    virtual FieldVector<Scalar,numEq> initial (const FieldVector<DT,dim>& x, const Entity& e,
                                               const FieldVector<DT,dim>& xi) const
    {
        FieldVector<Scalar,numEq> values(0.0);
        for (int nEq=0; nEq < numEq; nEq+=(numEq/numCont)){

            if (nEq<3){
                values[0]=1.013e5 + (depthBOR_ - x[1])* 1045 * 9.81;
                values[1]=0.0;
                values[2]=283.15 + (depthBOR_ - x[1])*0.03;
            }
            else {
                values[nEq]=values[0]*0.99;
                values[nEq+1]=0.0;
                values[nEq+2]=values[2]*0.99;
            }
        }

        return values;
    }
    //  MINC implementation
    FieldVector<Scalar, numCont> initialPhaseState (const FieldVector<DT,dim>& x, const Entity& e,
                                                    const FieldVector<DT,dim>& xi) const
    {
        enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
        FieldVector<Scalar, numCont> state;
        for (int nC=0; nC<numCont; nC++){
            state[nC] = waterPhase;
        }
        return state;
    }

    //    int initialPhaseState (const FieldVector<DT,dim>& x, const Entity& e,
    //                      const FieldVector<DT,dim>& xi) const
    //        {
    //            enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
    //            int stateF;
    //            FieldVector<Scalar, numCont> state;
    //            stateF = waterPhase;
    //            for (int nC=0; nC<numCont; nC++){
    //                state[nC] = waterPhase;
    //            }
    //            Scalar stateF_test = state[0];
    //            return stateF;
    //        }





    FieldVector<Scalar,dim> gravity () const
    {
        FieldVector<Scalar,dim> values(0.0);

        values[1] = -9.81;

        return values;
    }


    BrineCO2Problem(Liquid_GL& liq, Gas_GL& gas, MincCO2Soil<Grid, Scalar>& soil,
                    TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>),
                    MultiComp& multicomp = *(new CWaterAir), Scalar depthBOR = 0.0)
        : TwoPTwoCNIProblem<Grid, Scalar, numEq>(liq, gas, soil, multicomp, law)
    {
        depthBOR_ = depthBOR;
    }

private:
    Scalar depthBOR_;
    //        Scalar soilDens_, soilHeatCp_, soilLDry_, soilLSw_;
};

}
#endif
