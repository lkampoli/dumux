// $Id: 2p2cnimodel.hh 736 2008-10-27 14:46:55Z melanie $

#ifndef DUNE_TWOPHASEHEATMODEL_HH
#define DUNE_TWOPHASEHEATMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/io/exporttodgf.hh"
#include <boost/format.hpp>

namespace Dune {
template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
		class FunctionType, class OperatorAssembler> class TwoPhaseHeatModel :
	public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
	{
public:
	typedef NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian,
	FunctionType, OperatorAssembler> ThisNonlinearModel;

	TwoPhaseHeatModel(const Grid& g, ProblemType& prob) :
		ThisNonlinearModel(g, prob), uOldTimeStep(g) {
	}

	TwoPhaseHeatModel(const Grid& g, ProblemType& prob, int level) :
		ThisNonlinearModel(g, prob, level), uOldTimeStep(g, level) {
	}

	virtual void initial() = 0;

	virtual void update(double& dt) = 0;

	virtual void solve() = 0;

	FunctionType uOldTimeStep;
};

template<class Grid, class Scalar, class ProblemType, class LocalJac, int numEq> class LeafP1TwoPhaseModel :
	public TwoPhaseHeatModel<Grid, Scalar, ProblemType, LocalJac,
		LeafP1Function<Grid, Scalar, numEq>, LeafP1OperatorAssembler<Grid, Scalar, numEq> > {
public:
	// define the function type:
	typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

	// define the operator assembler type:
	typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

	typedef TwoPhaseHeatModel<Grid, Scalar, ProblemType, LocalJac,
	FunctionType, OperatorAssembler> ThisTwoPhaseHeatModel;

	typedef LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJac, numEq> ThisType;

	typedef LocalJac LocalJacobian;

	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
	typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
	typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator
			IntersectionIterator;

	LeafP1TwoPhaseModel(const Grid& g, ProblemType& prob) :
		ThisTwoPhaseHeatModel(g, prob), problem(prob), _grid(g), vertexmapper(g,	g.leafIndexSet()), size((*(this->u)).size())
		{
	}

        virtual void update(double &dt) {
            DUNE_THROW(NotImplemented, "This method is obsolete. Use updateModel()!");
        }

	virtual void initial() {
		typedef typename Grid::Traits::template Codim<0>::Entity Entity;
		typedef typename Grid::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = Grid::dimension};
		enum {dimworld = Grid::dimensionworld};

		const GV& gridview(grid().leafView());

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
							1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				// initialize cell concentration
				(*(this->u))[globalId] = this->problem.initial(
						global, entity, local);
			}
			this->localJacobian().clearVisited();
			this->localJacobian().initiateStaticData(entity);
		}

		// set Dirichlet boundary conditions
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
							1);
			int size = sfs.size();

			// set type of boundary conditions
			this->localJacobian().template assembleBC<LeafTag>(entity);

			IntersectionIterator
					endit = IntersectionIteratorGetter<Grid, LeafTag>::end(entity);
			for (IntersectionIterator is = IntersectionIteratorGetter<Grid,
					LeafTag>::begin(entity); is!=endit; ++is)
				if (is->boundary()) {
					for (int i = 0; i < size; i++)
						// handle subentities of this face
						for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
							if (sfs[i].entity()
									== ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
											j, sfs[i].codim())) {
								for (int equationNumber = 0; equationNumber<numEq; equationNumber++) {
									if (this->localJacobian().bc(i)[equationNumber]
											== BoundaryConditions::dirichlet) {
										// get cell center in reference element
										Dune::FieldVector<DT,dim>
												local = sfs[i].position();

										// get global coordinate of cell center
										Dune::FieldVector<DT,dimworld>
												global = it->geometry().global(local);

										int
												globalId = vertexmapper.template map<dim>(
														entity, sfs[i].entity());
										FieldVector<int,numEq> dirichletIndex;
										FieldVector<BoundaryConditions::Flags, numEq>
												bctype = this->problem.bctype(
														global, entity, is,
														local);
												this->problem.dirichletIndex(global, entity, is,
														local, dirichletIndex);

										if (bctype[equationNumber]
												== BoundaryConditions::dirichlet) {
											FieldVector<Scalar,numEq>
													ghelp = this->problem.g(
															global, entity, is,
															local);
											(*(this->u))[globalId][dirichletIndex[equationNumber]]
													= ghelp[dirichletIndex[equationNumber]];
										}
									}
								}
							}
				}
		}

		*(this->uOldTimeStep) = *(this->u);
		return;
	}


    virtual double computeFlux ()
     {
		  typedef typename Grid::Traits::template Codim<0>::Entity Entity;
		  typedef typename Grid::ctype DT;
		  typedef typename GV::template Codim<0>::Iterator Iterator;
		  enum{dim = Grid::dimension};
		  enum{dimworld = Grid::dimensionworld};
	   	  double sign;
	   	  const GV& gridview(_grid.leafView());
		  Iterator eendit = gridview.template end<0>();
		  FieldVector<Scalar,numEq> flux(0);
		  double Flux(0);

		  for (Iterator it = gridview.template begin<0>(); it != eendit; ++it) // loop over all entities
		  {

			  	// get geometry type
			  	Dune::GeometryType gt = it->geometry().type();

			  	// get entity
			  	const Entity& entity = *it;

				FVElementGeometry<Grid> fvGeom;
			    fvGeom.update(entity);

				for (int k = 0; k < fvGeom.numEdges; k++)
			     {
				    int i = fvGeom.subContVolFace[k].i;

				    int j = fvGeom.subContVolFace[k].j;

				    int flag_i, flag_j;

				    // 2D case: give y or x value of the line over which flux is to be
				    //			calculated.
				    // up to now only flux calculation to lines or planes (3D) parallel to
				    // x, y and z axis possible

				    // Flux across plane with z = 80 numEq
				     if(fvGeom.subContVol[i].global[2] < 80.)
			   			  flag_i = 1;
			   		  else flag_i = -1;

			   		  if(fvGeom.subContVol[j].global[2] < 80.)
			   			  flag_j = 1;
			   		  else flag_j = -1;

			   		  if(flag_i == flag_j)
			   		   {
			   			  sign = 0;
			   		   }
			   		  else {
			   			  	if(flag_i > 0)
			   			  		sign = -1;
			   			  	else sign = 1; }

			   			  // get variables

			   		  if(flag_i != flag_j)
			   		  {
						this->localJacobian().setLocalSolution(entity);
						this->localJacobian().computeElementData(entity);
						this->localJacobian().updateVariableData(entity, this->localJacobian().u);


						flux = this->localJacobian().computeA(entity, this->localJacobian().u, k);
						Flux += sign*flux[1];
			   		  }
			     }

		  }
		  return Flux; // Co2 flux
     }


	virtual double totalCO2Mass() {
		typedef typename Grid::Traits::template Codim<0>::Entity Entity;
		typedef typename Grid::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = Grid::dimension};
		enum {dimworld = Grid::dimensionworld};
		enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};	// Phase state
		enum {numCont = numEq/3};
		const GV& gridview(_grid.leafView());
		double totalMass = 0;
		double minSat = 1e100;
		double maxSat = -1e100;
		double minP  = 1e100;
		double maxP = -1e100;
		double minTe = 1e100;
		double maxTe = -1e100;
		double minX = 1e100;
		double maxX = -1e100;

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			FVElementGeometry<Grid> fvGeom;
			fvGeom.update(entity);

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
							1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				int stateF;
				stateF = this->localJacobian().sNDat[globalId].phaseStateF;
				FieldVector<Scalar, numCont> state;
				for (int nC=0; nC<numCont; nC++){
					state[nC]=this->localJacobian().sNDat[globalId].phaseState[nC];
				}
				Scalar vol = fvGeom.subContVol[i].volume;
//				Scalar poro = this->problem.soil().porosity(global, entity, local);
				Scalar poro = this->problem.soil().porosity(global, entity, local);
				Scalar rhoN = (*(this->localJacobian().outDensityN))[globalId];
				Scalar rhoW = (*(this->localJacobian().outDensityW))[globalId];
				Scalar satN = (*(this->localJacobian().outSaturationN))[globalId];
				Scalar satW = (*(this->localJacobian().outSaturationW))[globalId];
				Scalar xAW = (*(this->localJacobian().outMassFracAir))[globalId];
				Scalar xWN = (*(this->localJacobian().outMassFracWater))[globalId];
				Scalar xAN = 1 - xWN;
				Scalar pW = (*(this->u))[globalId][0];
				Scalar Te = (*(this->u))[globalId][2];
				Scalar mass = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);



				minSat = std::min(minSat, satN);
				maxSat = std::max(maxSat, satN);
				minP = std::min(minP, pW);
				maxP = std::max(maxP, pW);
				minX = std::min(minX, xAW);
				maxX = std::max(maxX, xAW);
				minTe = std::min(minTe, Te);
				maxTe = std::max(maxTe, Te);

				totalMass += mass;
			}

		}

		// print minimum and maximum values
		std::cout << "nonwetting phase saturation: min = "<< minSat
				<< ", max = "<< maxSat << std::endl;
		std::cout << "wetting phase pressure: min = "<< minP
				<< ", max = "<< maxP << std::endl;
		std::cout << "mole fraction CO2: min = "<< minX
				<< ", max = "<< maxX << std::endl;
		std::cout << "temperature: min = "<< minTe
				<< ", max = "<< maxTe << std::endl;

		return totalMass;
	}

	virtual void globalDefect(FunctionType& defectGlobal) {
		typedef typename Grid::Traits::template Codim<0>::Entity Entity;
		typedef typename Grid::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = Grid::dimension};
		typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

		const GV& gridview(_grid.leafView());
		(*defectGlobal)=0;
		// allocate flag vector to hold flags for essential boundary conditions
		std::vector<BCBlockType> essential(this->vertexmapper.size());
		for (typename std::vector<BCBlockType>::size_type i=0; i
				<essential.size(); i++)
			essential[i].assign(BoundaryConditions::neumann);
		// iterate through leaf grid
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;
			this->localJacobian().fvGeom.update(entity);
			int size = this->localJacobian().fvGeom.numVertices;
			this->localJacobian().setLocalSolution(entity);
			this->localJacobian().computeElementData(entity);
			bool old = true;
			this->localJacobian().updateVariableData(entity, this->localJacobian().uold, old);
			this->localJacobian().updateVariableData(entity, this->localJacobian().u);
			this->localJacobian().template localDefect<LeafTag>(entity, this->localJacobian().u);

			// begin loop over vertices
			for (int i=0; i < size; i++) {
				int globalId = this->vertexmapper.template map<dim>(entity,i);
				for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
					if (this->localJacobian().bc(i)[equationnumber] == BoundaryConditions::neumann)
						(*defectGlobal)[globalId][equationnumber]
								+= this->localJacobian().def[i][equationnumber];
					else
						essential[globalId].assign(BoundaryConditions::dirichlet);
				}
			}
		}

		for (typename std::vector<BCBlockType>::size_type i=0; i
				<essential.size(); i++)
			for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
			if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
				(*defectGlobal)[i][equationnumber] = 0;
			}
	}


	virtual void vtkout(const char* name, int k) {
	}

	void writerestartfile(int restartNum=0)
	{
		enum {dim = Grid::dimension};
		typedef typename GV::template Codim<dim>::Iterator Iterator;

//		exportToDGF(_grid.leafView(), *(this->u), numEq, "primvar", false);

		const int size = vertexmapper.size();
		BlockVector<FieldVector<double, numEq+1> > data(size);
		data=0;

		Iterator endIt = _grid.leafView().template end<dim>();
		for (Iterator it = _grid.leafView().template begin<dim>(); it != endIt;	++it)
		{
			int index = vertexmapper.map(*it);
			for (int i = 0; i < numEq;i++)
			{
				data[index][i]=(*(this->u))[index][i];
			}
			data[index][numEq]=this->localJacobian().sNDat[index].phaseStateF;
		}
		exportToDGF(_grid.leafView(), data, (numEq+1), "data", false);
	}
    const Grid &grid() const
        { return _grid; }

protected:
  ProblemType& problem;
  const Grid& _grid;
  VertexMapper vertexmapper;
  int size;
};

}
#endif
