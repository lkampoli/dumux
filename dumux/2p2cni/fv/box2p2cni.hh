// $Id$

#ifndef DUNE_BOX2P2CNI_HH
#define DUNE_BOX2P2CNI_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

//#include <dune/common/array.hh>        // defines simple array class
#include <dune/common/fixedarray.hh>   // defines simple array classes
#include <dune/common/geometrytype.hh>
#include <dune/grid/sgrid.hh>          // a complete structured grid
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/universalmapper.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/vbvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/gsetc.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/disc/functions/functions.hh>
#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/2p2cni/2p2cnimodel.hh"
#include "dumux/2p2cni/2p2cniproblem.hh"
#include "dumux/2p2cni/fv/box2p2cnijacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/new_models/2p2cni/2p2cninewtoncontroller.hh"
#include "dumux/io/importfromdgf_leaf.hh"
#include <boost/format.hpp>

namespace Dune
{
/**
   \brief Non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns

   This implements a non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns
*/
template<class ThisGrid, class ThisScalar, class VtkMultiWriter>
class Box2P2CNI: public LeafP1TwoPhaseModel<ThisGrid, ThisScalar,
                                            TwoPTwoCNIProblem<ThisGrid, ThisScalar> , Box2P2CNIJacobian<ThisGrid,
                                                                                                        ThisScalar> >
{
public:
    // define the problem type (also change the template argument above)
    typedef TwoPTwoCNIProblem<ThisGrid, ThisScalar> ProblemType;
    // define the localPos Jacobian (also change the template argument above)
    typedef Box2P2CNIJacobian<ThisGrid, ThisScalar> LocalJacobian;
    typedef LeafP1TwoPhaseModel<ThisGrid, ThisScalar, ProblemType,
                                LocalJacobian> ThisLeafP1TwoPhaseModel;
    typedef Box2P2CNI<ThisGrid, ThisScalar, VtkMultiWriter> ThisType;

    typedef    typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename ThisGrid::LeafGridView GridView;

    enum
        {   numEq = 3};
    enum
        {   wComp = 0, nComp = 1, temp = 2};
    enum
        {   gasPhase = 0, waterPhase = 1, bothPhases = 2};

    typedef typename ThisLeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;

    //////////////////////
    // Stuff required for the new newton method

    //! The traits class for the new newton method.
    struct NewtonTraits
    {
        typedef ThisScalar Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType Function;
        typedef typename ThisType::LocalJacobian LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits
    {
        typedef ThisScalar Scalar;
        typedef ThisGrid Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod> NewtonController;

    typedef typename NewtonTraits::Function Function;
    Function &currentSolution()
    {   return this->u;};

    LocalJacobian &getLocalJacobian()
    {   return this->localJacobian();}

    typedef typename NewtonTraits::JacobianAssembler JacobianAssembler;
    JacobianAssembler &jacobianAssembler()
    {   return this->A;}
    // End of stuff for new newton method
    //////////////////////


    Box2P2CNI(const ThisGrid& grid, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(grid, prob)// (this->size) vectors

    {}

    void initial()
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum
        {   dim = ThisGrid::dimension};
        enum
        {   dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        // iterate through leaf grid an evaluate c0 at cell center
        ElementIterator eendit = gridview.template end<0>();
        for (ElementIterator eIt = gridview.template begin<0>(); eIt
                 != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            this->localJacobian().fvGeom.update(element);

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,1);
            int size = sfs.size();

            for (int idx = 0; idx < size; idx++)
            {
                // get cell center in reference element
                const Dune::FieldVector<CoordScalar,dim>&localPos = sfs[idx].position();

                // get global coordinate of cell center
                Dune::FieldVector<CoordScalar,dimworld> globalPos = eIt->geometry().global(localPos);

                int globalIdx = this->vertexmapper.template map<dim>(element,
                                                                     sfs[idx].entity());

                // initialize cell concentration
                (*(this->u))[globalIdx] = this->problem.initial(
                                                                globalPos, element, localPos);

                // initialize variable phaseState
                this->localJacobian().sNDat[globalIdx].phaseState =
                    this->problem.initialPhaseState(globalPos, element, localPos);
                // initialize variable oldPhaseState
                this->localJacobian().sNDat[globalIdx].oldPhaseState =
                    this->problem.initialPhaseState(globalPos, element, localPos);

            }
            this->localJacobian().clearVisited();
            this->localJacobian().initiateStaticData(element);
        }

        // set Dirichlet boundary conditions
        for (ElementIterator eIt = gridview.template begin<0>(); eIt
                 != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,
                                                                                         1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(element);

            IntersectionIterator
                endit = element.ileafend();
            for (IntersectionIterator isIt = element.ileafbegin(); isIt!=endit; ++isIt)
                if (isIt->boundary())
                {
                    for (int idx = 0; idx < size; idx++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(isIt->indexInInside(), 1, sfs[idx].codim()); j++)
                            if (sfs[idx].entity()
                                == ReferenceElements<CoordScalar,dim>::general(gt).subEntity(isIt->indexInInside(), 1,
                                                                                             j, sfs[idx].codim()))
                            {
                                for (int equationNumber = 0; equationNumber<numEq; equationNumber++)
                                {
                                    if (this->localJacobian().bc(idx)[equationNumber]
                                        == BoundaryConditions::dirichlet)
                                    {
                                        // get cell center in reference element
                                        Dune::FieldVector<CoordScalar,dim>
                                            localPos = sfs[idx].position();

                                        // get globalPos coordinate of cell center
                                        Dune::FieldVector<CoordScalar,dimworld>
                                            globalPos = eIt->geometry().global(localPos);

                                        int
                                            globalIdx = this->vertexmapper.template map<dim>(
                                                                                             element, sfs[idx].entity());
                                        FieldVector<int,numEq> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, numEq>
                                            bctype = this->problem.bctype(
                                                                          globalPos, element, isIt,
                                                                          localPos);
                                        this->problem.dirichletIndex(globalPos, element, isIt,
                                                                     localPos, dirichletIndex);

                                        if (bctype[equationNumber]
                                            == BoundaryConditions::dirichlet)
                                        {
                                            FieldVector<ThisScalar,numEq>
                                                ghelp = this->problem.g(
                                                                        globalPos, element, isIt,
                                                                        localPos);
                                            (*(this->u))[globalIdx][dirichletIndex[equationNumber]]
                                                = ghelp[dirichletIndex[equationNumber]];
                                        }
                                    }
                                }
                            }
                }
            this->localJacobian().setLocalSolution(element);
            for (int idx = 0; idx < size; idx++)
                this->localJacobian().updateVariableData(element, this->localJacobian().u, idx, false);

        }

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    void restart(int restartNum=0)
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum
        {   dim = ThisGrid::dimension};
        enum
        {   dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, numEq+1> > data(size);
        data=0;

        // initialize primary variables
        std::string restartFileName;
        restartFileName = (boost::format("data-%05d")
                           %restartNum).str();
        importFromDGF<GridView>(data, restartFileName, false);
        this->localJacobian().clearVisited();

        for (int globalIdx=0;globalIdx<size;globalIdx++)
        {
            for (int equationNumber=0;equationNumber<numEq;equationNumber++)
            {
                (*(this->u))[globalIdx][equationNumber]=data[globalIdx][equationNumber];
            }

            //////////////////////////// SEQUENTIAL COUPLING WITH 2pni model ////////////////////
            // determine phase state according to saturation calculated by 2pni model
            // For bothphases case saturation needs to be corrected to prevent mass errors
            // see correctSaturation function in boxco2jacobian
            /////////////////////////////////////////////////////////////////////////////////////

            if(this->problem.sequentialCoupling() == true)
            {
                if((*(this->u))[globalIdx][nComp] <= 1.e-10)
                {
                    // initialize variable phaseState
                    this->localJacobian().sNDat[globalIdx].phaseState = waterPhase;
                    // initialize variable oldPhaseState
                    this->localJacobian().sNDat[globalIdx].oldPhaseState = waterPhase;
                    (*(this->u))[globalIdx][nComp] = 0.0;
                    this->localJacobian().sNDat[globalIdx].primVarSet = true;
                }
                else if((*(this->u))[globalIdx][nComp] >= 1-1.e-10)
                {
                    // initialize variable phaseState
                    this->localJacobian().sNDat[globalIdx].phaseState = gasPhase;
                    // initialize variable oldPhaseState
                    this->localJacobian().sNDat[globalIdx].oldPhaseState = gasPhase;
                    (*(this->u))[globalIdx][nComp] = 0.0;
                    this->localJacobian().sNDat[globalIdx].primVarSet = true;
                }
                else if(this->localJacobian().sNDat[globalIdx].primVarSet == false)
                {
                    // initialize variable phaseState
                    this->localJacobian().sNDat[globalIdx].phaseState = bothPhases;
                    // initialize variable oldPhaseState
                    this->localJacobian().sNDat[globalIdx].oldPhaseState = bothPhases;
                }
            }
        }

        // iterate through leaf grid an evaluate c0 at cell center
        ElementIterator eendit = gridview.template end<0>();

        if(this->problem.sequentialCoupling() == false)
        {
            for (ElementIterator eIt = gridview.template begin<0>(); eIt
                     != eendit; ++eIt)
            {
                // get geometry type
                Dune::GeometryType gt = eIt->geometry().type();

                // get element
                const Element& element = *eIt;

                this->localJacobian().fvGeom.update(element);

                const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
                    &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,1);
                int size = sfs.size();

                for (int idx = 0; idx < size; idx++)
                {
                    // get cell center in reference element
                    const Dune::FieldVector<CoordScalar,dim>&localPos = sfs[idx].position();

                    // get globalPos coordinate of cell center
                    Dune::FieldVector<CoordScalar,dimworld> globalPos = eIt->geometry().global(localPos);

                    int globalIdx = this->vertexmapper.template map<dim>(element,
                                                                         sfs[idx].entity());

                    // initialize variable phaseState
                    this->localJacobian().sNDat[globalIdx].phaseState = data[globalIdx][numEq];
                    // initialize variable oldPhaseState
                    this->localJacobian().sNDat[globalIdx].oldPhaseState = data[globalIdx][numEq];

                }
                this->localJacobian().initiateStaticData(element);
            }
        }

        ////////////////////////////// INCLUDE BOUNDARY CONDITIONS ////////////////////////////
        // Take care! Boundary conditions are only included for normal restart
        // for sequential coupling boundary conditions are taken from restart file
        // in the first time step because boundary saturation values might be corrected
        // accidentially
        ///////////////////////////////////////////////////////////////////////////////////////
        for (ElementIterator eIt = gridview.template begin<0>(); eIt != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get entity
            const Element& element = *eIt;

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,
                                                                                         1);
            int size = sfs.size();

            if(this->problem.sequentialCoupling() == false)
            {
                // set type of boundary conditions for normal restart
                this->localJacobian().assembleBoundaryCondition(element);

                IntersectionIterator
                    endit = element.ileafend();
                for (IntersectionIterator isIt = element.ileafbegin(); isIt!=endit; ++isIt)
                    if (isIt->boundary())
                    {
                        for (int idx = 0; idx < size; idx++)
                            // handle subentities of this face
                            for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(isIt->indexInInside(), 1, sfs[idx].codim()); j++)
                                if (sfs[idx].entity()
                                    == ReferenceElements<CoordScalar,dim>::general(gt).subEntity(isIt->indexInInside(), 1,
                                                                                                 j, sfs[idx].codim()))
                                {
                                    for (int equationNumber = 0; equationNumber<numEq; equationNumber++)
                                    {
                                        if (this->localJacobian().bc(idx)[equationNumber]
                                            == BoundaryConditions::dirichlet)
                                        {
                                            // get cell center in reference element
                                            Dune::FieldVector<CoordScalar,dim>
                                                localPos = sfs[idx].position();

                                            // get globalPos coordinate of cell center
                                            Dune::FieldVector<CoordScalar,dimworld>
                                                globalPos = eIt->geometry().global(localPos);

                                            int
                                                globalIdx = this->vertexmapper.template map<dim>(
                                                                                                 element, sfs[idx].entity());
                                            FieldVector<int,numEq> dirichletIndex;
                                            FieldVector<BoundaryConditions::Flags, numEq>
                                                bctype = this->problem.bctype(
                                                                              globalPos, element, isIt,
                                                                              localPos);
                                            this->problem.dirichletIndex(globalPos, element, isIt,
                                                                         localPos, dirichletIndex);

                                            if (bctype[equationNumber]
                                                == BoundaryConditions::dirichlet)
                                            {
                                                FieldVector<ThisScalar,numEq>
                                                    ghelp = this->problem.g(
                                                                            globalPos, element, isIt,
                                                                            localPos);
                                                (*(this->u))[globalIdx][dirichletIndex[equationNumber]]
                                                    = ghelp[dirichletIndex[equationNumber]];
                                            }
                                        }
                                    }
                                }
                    }
            }

            this->localJacobian().setLocalSolution(element);
            for (int idx = 0; idx < size; idx++)
            {
                // write data into output vector
                // if sequenatialCoupling == true calculate secondary variables
                // needed for correctSaturation function
                this->localJacobian().updateVariableData(element, this->localJacobian().u, idx, false);
            }

            /////////////////////////////// SEQUENTIAL COUPLING - saturation correction /////////
            if(this->problem.sequentialCoupling() == true)
            {
                this->localJacobian().correctSaturation(element, size, this->localJacobian().u);
                this->localJacobian().initiateStaticData(element);

                for (int idx = 0; idx < size; idx++)
                {
                    /// write corrected data into output vector
                    this->localJacobian().updateVariableData(element, this->localJacobian().u, idx, false);
                }
            }
        }

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    void updateState()
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum
        {   dim = ThisGrid::dimension};
        enum
        {   dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());
        // iterate through leaf grid and evaluate c0 at cell center
        ElementIterator eendit = gridview.template end<0>();
        for (ElementIterator eIt = gridview.template begin<0>(); eIt
                 != eendit; ++eIt)
        {

            const Element& element = *eIt;
            this->localJacobian().fvGeom.update(element);
            this->localJacobian().setLocalSolution(element);
            this->localJacobian().computeElementData(element);
            this->localJacobian().updateStaticData(element, this->localJacobian().u);
            this->localJacobian().localToGlobal(element);
        }
        return;
    }

    virtual void globalDefect(FunctionType& defectGlobal)
    {
        ThisLeafP1TwoPhaseModel::globalDefect(defectGlobal);
    }

    void solve()
    {
        Operator op(*(this->A)); // make operator out of matrix
        double red=1E-14;

#ifdef HAVE_PARDISO
        SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
        pardiso.factorize(*(this->A));
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2); // an inverse operator
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1); // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }

    void update(double& dt)
    {
        DUNE_THROW(NotImplemented, "the update method is deprecated. use updateModel()");
    }

    void updateModel(double& dt, double &nextDt)
    {
        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        enum
        {   dim = ThisGrid::dimension};

        ///////////////////////////////////
        // define solver tolerances here
        ///////////////////////////////////

        /*
//////////////
ThisScalar absTol = 1e-1;
ThisScalar relTol = 1e-8;
///////////////////
NewtonMethod<ThisGrid, ThisType> newtonMethod(this->_grid, *this, relTol, absTol);
newtonMethod.execute();
        */
        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian().getDt();
            NewtonMethod newton(*this);
            NewtonController newtonCtl;
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian().setDt(nextDt);
            if(!newtonLoop)
            {
                *this->u = *this->uOldTimeStep;
                this->localJacobian().resetPhaseState();
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        double Flux(0), Mass(0);
        Flux = this->computeFlux();
        Mass = this->totalCO2Mass();
        std::cout << Flux << ", "<< Mass<< "    # flux nPhase, total mass nPhase \n";

        this->localJacobian().updatePhaseState(); // update variable oldPhaseState
        // this->localJacobian().clearVisited();
        // updateState();                            // phase switch after each timestep


        *(this->uOldTimeStep) = *(this->u);
        // update old phase state for computation of ComputeM(..uold..)

        return;
    }

    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {
        //        BlockVector<FieldVector<ThisScalar, 1> > &xWN = *writer.template createField<ThisScalar, 1>(this->size);
        //        BlockVector<FieldVector<ThisScalar, 1> > &xAW = *writer.template createField<ThisScalar, 1>(this->size);
        //        BlockVector<FieldVector<ThisScalar, 1> > &satW = *writer.template createField<ThisScalar, 1>(this->size);

        //        writer.addScalarVertexFunction("nonwetting phase saturation", this->u, 1);
        writer.addScalarVertexFunction("pressure wetting phase", this->u, 0);

        ////        writer.addScalarVertexFunction("nonwetting phase saturation",
        ////                                        this->u,
        ////                                        1);
        //        writer.addScalarVertexFunction("wetting phase pressure",
        //                                        this->u,
        //                                        0);
        //        writer.addVertexData(&satW,"wetting phase saturation");
        writer.addVertexData(this->localJacobian().outPressureN,"pressure non-wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outTemperature,"temperature");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracAir,"massfraction nonwetting component in wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracWater,"massfraction wetting component in non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
        writer.addVertexData(this->localJacobian().outPhaseState,"phase state");

    }

    void vtkout (const char* name, int k)
    {

    }

    void setVtkMultiWriter(VtkMultiWriter *writer)
    {   vtkMultiWriter = writer;}
protected:
    VtkMultiWriter *vtkMultiWriter;

};
}
#endif
