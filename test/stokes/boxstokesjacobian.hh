#ifndef DUNE_BOXSTOKESJACOBIAN_HH
#define DUNE_BOXSTOKESJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/operators/boxjacobian.hh>
#include "dumux/stokes/stokesproblem.hh"

namespace Dune
{
  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
    diffusion equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses the box method.

    Template parameters are:

    - G     a DUNE grid type
    - RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, G::dimension+1> >
  class BoxStokesJacobian
    : public BoxJacobian<BoxStokesJacobian<G,RT,BoxFunction>,G,RT,G::dimension+1,BoxFunction>
  {
        enum {dim=G::dimension};
        enum {numEq = dim+1};

    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxStokesJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,numEq>::VBlockType VBlockType;
    typedef Dune::FVElementGeometry<G> FVElementGeometry;

  public:

    //! Constructor
    BoxStokesJacobian (StokesProblem<G,RT>& params,
                  bool levelBoundaryAsDirichlet_, const G& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,numEq,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
        VBlockType result(0);

        return result;
    }

    VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
    	VBlockType result = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);

    	RT flux = boundaryFlux(e, sol, node)/this->fvGeom.subContVol[node].volume;
    	result[dim] -= flux;

       return (result);
    }

    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
    	VBlockType flux(0);

    	RT pressValue = 0;
    	FieldVector<RT, dim> velocityValue(0);
		for (int k = 0; k < this->fvGeom.numVertices; k++) {
			pressValue += sol[k][dim]*this->fvGeom.subContVolFace[face].shapeValue[k];
			for (int comp = 0; comp < dim; comp++)
				velocityValue[comp] += sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
		}

		// momentum balance:
		for (int comp = 0; comp < dim; comp++)
    	{
    		FieldVector<RT, dim> gradVComp(0);
    		for (int k = 0; k < this->fvGeom.numVertices; k++) {
    			FieldVector<RT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
    			grad *= sol[k][comp];
    			gradVComp += grad;
    		}

    		gradVComp *= -elData.mu;

    		FieldVector<RT, dim> pComp(0);
    		pComp[comp] = pressValue;

    		gradVComp += pComp;

    		flux[comp] = gradVComp*this->fvGeom.subContVolFace[face].normal;
    	}

		// mass balance:
		flux[dim] = velocityValue*this->fvGeom.subContVolFace[face].normal;

    	return flux;
    }

    void computeElementData (const Entity& e)
    {
        // ASSUMING element-wise constant viscosity, evaluate mu at the cell center
        elData.mu = problem.mu(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
    };

     virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
         return;
    }

    void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
    {
        return;
    }

    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
        return;
    }

    RT boundaryFlux(const Entity& e, const VBlockType* sol, int node) {
    	RT result = 0;

        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
        &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(e);

        const typename ReferenceElementContainer<DT,dim>::value_type
        &referenceElement = ReferenceElements<DT, dim>::general(gt);

        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<G, LeafTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G, LeafTag>::begin(e); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (this->levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                    continue;
                if (!this->levelBoundaryAsDirichlet)
                    continue;
            }

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->numberInSelf();

                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                    if (node != nodeInElement)
                    	continue;

                    int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);

                    // get geometry type of face
                    GeometryType faceGT = it->intersectionSelfLocal().type();

                    // center in face's reference element
                    const FieldVector<RT,dim-1>& faceLocal = ReferenceElements<RT,dim-1>::general(faceGT).position(0,0);

                    FieldVector<DT,dim>  velocityValue(0);
                    for (int vert = 0; vert < this->fvGeom.numVertices; vert++)
                    	for (int comp = 0; comp < dim; comp++)
                    		velocityValue[comp] += sol[vert][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];

                    result += velocityValue*it->outerNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                }
            }

        }

        return result;
    }


    template<class TypeTag> void assembleBC(const Entity& e) {
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
        &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(e);

        const typename ReferenceElementContainer<DT,dim>::value_type
        &referenceElement = ReferenceElements<DT, dim>::general(gt);

        for (int i = 0; i < sfs.size(); i++) {
            this->bctype[i].assign(BoundaryConditions::neumann);
            this->b[i] = 0;
            this->dirichletIndex[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<G, TypeTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G, TypeTag>::begin(e); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (this->levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                    continue;
                if (!this->levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
            FieldVector<int,numEq> dirichletIdx(0);

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->numberInSelf();
                //                 std::cout << "faceIdx = " << faceIdx << ", beginning: " << std::endl;
                //                 for (int i = 0; i < 4; i++)
                //                   std::cout << "bctype[" << i << "] = " << this->bctype[i] << std::endl;

                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                    for (int equationNumber = 0; equationNumber < numEq; equationNumber++) {
                        if (this->bctype[nodeInElement][equationNumber] == BoundaryConditions::neumann) {
                            int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);
                            FieldVector<DT,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                            FieldVector<DT,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                            bctypeface = this->getImp().problem.bctype(global, e, it, local); // eval bctype
                            this->getImp().problem.dirichletIndex(global, e, it, local, dirichletIdx); // eval bctype
                            //                                                     std::cout << "faceIdx = " << faceIdx << ", nodeInElement = " << nodeInElement
                            //                                                           << ", bfIdx = " << bfIdx << ", local = " << local << ", global = " << global
                            //                                                           << ", bctypeface = " << bctypeface << std::endl;
                            if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                break;
                            FieldVector<DT,dim> J = this->getImp().problem.J(global, e, it, local);
                            if (equationNumber < dim) {
                            	J[equationNumber] *= this->fvGeom.boundaryFace[bfIdx].area;
                            	this->b[nodeInElement][equationNumber] += J[equationNumber];
                            }
                        }
                    }
                }

                bool nface(true); // check if face is a neumann face
                for(int i=0; i<numEq; i++)
                {
                    if(bctypeface[i] != BoundaryConditions::neumann)
                        nface = false; // was not a neumann face
                }
                if(nface == true)
                    continue; // was a neumann face, go to next face
            }

            // If we are here, then it is
            // (i)   an exterior boundary face with Dirichlet condition, or
            // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
            // (iii) a level boundary in case of level-wise assemble
            // How processor boundaries are handled depends on the processor boundary mode

            bool pface(false);  // check if face is a process boundary
            for(int i=0; i<numEq; i++)
            {
                if (bctypeface[i]==BoundaryConditions::process
                        && this->procBoundaryAsDirichlet==false
                        && this->levelBoundaryAsDirichlet==false)
                {
                    pface = true;
                    break;
                }
            }
            if(pface == true)
                continue;   // if face is a process boundary it acts like homogeneous Neumann


            for (int equationNumber=0; equationNumber<dim; equationNumber++) {
                for (int i=0; i<sfs.size(); i++) // loop over test function number
                {
                    //this->dirichletIndex[i][equationNumber] = equationNumber;

                    //std::cout<<"i = "<<i<<std::endl;
                    if (sfs[i].codim()==0)
                        continue; // skip interior dof
                    if (sfs[i].codim()==1) // handle face dofs
                    {
                        if (sfs[i].entity()==it->numberInSelf()) {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

                                if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                        continue;
                    }
                    // handle subentities of this face
                    for (int j=0; j<ReferenceElements<DT,dim>::general(gt).size(it->numberInSelf(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity()==ReferenceElements<DT,dim>::general(gt).subEntity(it->numberInSelf(), 1, j, sfs[i].codim()))
                        {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
                                if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                }
            }
        }

    }


    struct ElementData {
        RT mu;
       };

       ElementData elData;
    StokesProblem<G,RT>& problem;
  };
}
#endif
