// $Id$

#ifndef DUNE_FVNONWETTINGVELOCITY2P_HH
#define DUNE_FVNONWETTINGVELOCITY2P_HH

#include "dumux/diffusion/fv/fvpressure2p.hh"
#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
/** \todo Please doc me! */

template<class GridView, class Scalar, class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> > class FVNonWettingPhaseVelocity2P: public FVPressure2P<
        GridView, Scalar, VC, Problem>
{

typedef    typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
    {   dim = GridView::dimension};
    enum
    {   dimWorld = GridView::dimensionworld};
    enum
    {
        pw = 0, pn = 1, pglobal = 2
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    FVNonWettingPhaseVelocity2P(GridView& gridView, Problem& problem, std::string pressureType, std::string satType = "Sw")
    : FVPressure2P<GridView,Scalar,VC, Problem>(gridView, problem, pressureType, satType)
    {}

    FVNonWettingPhaseVelocity2P(GridView& gridView, Problem& problem, std::string pressureType, std::string satType, std::string solverName,
            std::string preconditionerName)
    : FVPressure2P<GridView,Scalar,VC, Problem>(gridView, problem, pressureType, satType, solverName, preconditionerName)
    {}

    void calculateVelocity(const Scalar t=0) const
    {
        // phase densities
        Scalar densityNW = this->diffProblem.nonWettingPhase().density();

        // compute update vector
        ElementIterator eItEnd = this->gridView.template end<0>();
        for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // cell center in reference element
            const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            //
            GlobalPosition globalPos = eIt->geometry().global(localPos);

            // cell index
            int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

            Scalar pressI = this->diffProblem.variables().pressure()[globalIdxI];
            Scalar pcI = this->diffProblem.variables().capillaryPressure()[globalIdxI];
            Scalar lambdaNWI = this->diffProblem.variables().mobilityNonWetting()[globalIdxI];

            // run through all intersections with neighbors and boundary
            IntersectionIterator
            isItEnd = this->gridView.template iend(*eIt);
            for (IntersectionIterator
                    isIt = this->gridView.template ibegin(*eIt); isIt
                    !=isItEnd; ++isIt)
            {
                // local number of facet
                int indexInInside = isIt->indexInInside();

                // get geometry type of face
                Dune::GeometryType faceGT = isIt->geometryInInside().type();

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>&
                faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                // center of face inside volume reference element
                const LocalPosition&
                localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(indexInInside,1);

                Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                // get absolute permeability
                FieldMatrix permeabilityI(this->diffProblem.soil().K(globalPos, *eIt, localPos));
                // compute directed permeability vector permeabilityI.n
                FieldVector<Scalar,dim> normalPermeabilityI(0);
                permeabilityI.umv(unitOuterNormal, normalPermeabilityI);

                // handle interior face
                if (isIt->neighbor())
                {
                    // access neighbor
                    ElementPointer neighborPointer = isIt->outside();
                    int globalIdxJ = this->diffProblem.variables().indexDiffusion(*neighborPointer);

                    // compute factor in neighbor
                    Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                    const LocalPosition&
                    localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                    // cell center in global coordinates
                    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                    // neighbor cell center in global coordinates
                    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // get absolute permeability
                    FieldMatrix permeabilityJ(this->diffProblem.soil().K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

                    // compute vectorized permeabilities
                    FieldVector<Scalar,dim> normalPermeabilityJ(0);
                    permeabilityJ.umv(unitOuterNormal, normalPermeabilityJ);
                    // compute permeability normal to intersection and take harmonic mean
                    Scalar normalComponentPermeabilityI = normalPermeabilityI * unitOuterNormal;
                    Scalar normalComponentPermeabilityJ = normalPermeabilityJ * unitOuterNormal;
                    Scalar meanNormalPermeability = 2 * normalComponentPermeabilityI * normalComponentPermeabilityJ / (normalComponentPermeabilityI + normalComponentPermeabilityJ);
                    // compute permeability tangential to intersection and take arithmetic mean
                    FieldVector<Scalar,dim> normalComponentVector = unitOuterNormal;
                    FieldVector<Scalar,dim> tangentialPermeabilityI = normalPermeabilityI - (normalComponentVector *= normalComponentPermeabilityI);
                    normalComponentVector = unitOuterNormal;
                    FieldVector<Scalar,dim> tangentialPermeabilityJ = normalPermeabilityJ - (normalComponentVector *= normalComponentPermeabilityJ);
                    FieldVector<Scalar,dim> meanTangentialPermeability = (tangentialPermeabilityI += tangentialPermeabilityJ);
                    meanTangentialPermeability *= 0.5;
                    FieldVector<Scalar,dim> meanNormalPermeabilityVector = unitOuterNormal;
                    // Build vectorized averaged permeability
                    FieldVector<Scalar,dim> permeability = (meanTangentialPermeability += (meanNormalPermeabilityVector *= meanNormalPermeability));

                    Scalar pressJ = this->diffProblem.variables().pressure()[globalIdxJ];
                    Scalar pcJ = this->diffProblem.variables().capillaryPressure()[globalIdxJ];
                    Scalar lambdaNWJ = this->diffProblem.variables().mobilityNonWetting()[globalIdxJ];

                    //determine upwind direction
                    Scalar potentialNW = this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaNW = (potentialNW >= 0.) ? lambdaNWI : lambdaNWJ;

                    FieldVector<Scalar,dimWorld> velocity(permeability);
                    FieldVector<Scalar,dimWorld> gravityTerm(this->gravity);
                    for (int i=0;i<dim;i++)
                    {
                        gravityTerm[i] *= permeability[i]*unitOuterNormal[i];
                    }

                    if (this->pressureType == pw)
                    {
                        velocity *= lambdaNW * (pressI - pressJ)/dist + 0.5 * (lambdaNWI + lambdaNWJ) * (pcI - pcJ) / dist;
                        velocity += (gravityTerm *= (lambdaNW * densityNW));
                    }
                    if (this->pressureType == pn)
                    {
                        velocity *= lambdaNW * (pressI - pressJ) / dist;
                        velocity += (gravityTerm *= (lambdaNW * densityNW));
                    }
                    if (this->pressureType == pglobal)
                    {
                        DUNE_THROW(NotImplemented, "Pressure type not supported for v_w!");
                    }
                    this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;

                }

                // handle boundary face
                if (isIt->boundary())
                {
                    // center of face in global coordinates
                    GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                    //get boundary type
                    BoundaryConditions::Flags bcTypeSat = this->diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);
                    BoundaryConditions::Flags bcTypePress = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);

                    // cell center in global coordinates
                    GlobalPosition globalPos = eIt->geometry().global(localPos);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    Scalar residualSatW = this->diffProblem.soil().Sr_w(globalPos, *eIt, localPos);
                    Scalar residualSatNW = this->diffProblem.soil().Sr_n(globalPos, *eIt, localPos);

                    Scalar satBound = 0;
                    if (bcTypeSat == BoundaryConditions::dirichlet)
                    {
                        satBound = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    }
                    else
                    {
                        satBound = this->diffProblem.variables().saturation()[globalIdxI];
                    }

                    if (bcTypePress == BoundaryConditions::dirichlet)
                    {
                        Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                        Scalar pcBound = this->diffProblem.materialLaw().pC(1-satBound, globalPosFace, *eIt, localPosFace);

                        Scalar lambdaNWBound = this->diffProblem.materialLaw().mobN(satBound,globalPosFace, *eIt, localPosFace);

                        Scalar potentialNW = this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside];

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaNW = (potentialNW >= 0.) ? lambdaNWI : lambdaNWBound;

                        FieldVector<Scalar,dimWorld> velocity(normalPermeabilityI);
                        FieldVector<Scalar,dimWorld> gravityTerm(this->gravity);
                        for (int i=0;i<dim;i++)
                        {
                            gravityTerm[i] *= normalPermeabilityI[i]*unitOuterNormal[i];
                        }

                        if (this->pressureType == pw)
                        {
                            velocity *= lambdaNW * (pressI - pressBound)/dist + 0.5 * (lambdaNWI + lambdaNWBound) * (pcI - pcBound) / dist;
                            velocity += (gravityTerm *= (lambdaNW * densityNW));
                        }
                        if (this->pressureType == pn)
                        {
                            velocity *= lambdaNW * (pressI - pressBound) / dist;
                            velocity += (gravityTerm *= (lambdaNW * densityNW));
                        }
                        if (this->pressureType == pglobal)
                        {
                            DUNE_THROW(NotImplemented, "Pressure type not supported for v_w!");
                        }

                        this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;
                    }
                    else
                    {
                        Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                        FieldVector<Scalar,dimWorld> velocity(unitOuterNormal);

                        if (J != 0)
                        {
                            if (satBound <= residualSatNW)
                            {
                                velocity *= 0;
                            }
                            else if (satBound >= (1-residualSatW))
                            {
                                velocity *= J;
                            }
                            else
                            {
                                DUNE_THROW(NotImplemented, "v_w can not be determined from v_t if both phases are present!");
                            }
                        }
                        else
                        {
                            velocity *= 0;
                        }
                        this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;
                    }
                }
            }// end all intersections
        }// end grid traversal
//                printvector(std::cout, this->diffProblem.variables().velocity(), "velocity", "row", 4, 1, 3);
        return;
    }
};
}
#endif
