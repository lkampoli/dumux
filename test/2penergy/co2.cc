#include "config.h"
#include <iostream>
#ifdef HAVE_UG
#include <iomanip>
//#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
//#include "dumux/twophase/problems/uniformtwophaseproblem.hh"
#include "dumux/2penergy/fv/boxpwsnte.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/readstarformat.cc"
#include "co2_HeatFlux.hh"

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(20);
    Dune::FieldVector<NumberType, dim> innerLowerLeft(1);
    innerLowerLeft[1] = 2;
    Dune::FieldVector<NumberType, dim> innerUpperRight(4);
    innerUpperRight[1] = 3;
    if (argc != 4) {
      std::cout << "usage: co2 basefilename tEnd dt" << std::endl;
      return 0;
    }
    	std::string arg1(argv[2]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
	std::string arg2(argv[3]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;

    // create a grid object
    typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    //typedef Dune::UGGrid<dim> GridType; 
	//typedef Dune::ALUSimplexGrid<dim,dim> GridType; 
	
    Dune::GridPtr<GridType> gridPointer(argv[1]);
    GridType& grid = *gridPointer;
    //readStarFormat(grid, argv[1]);
    //grid.createLGMGrid(argv[1]);

     Dune::gridinfo(grid);
    
    Brine brine;
    CO2 co2;

//     Dune::LinearLaw law(brine, co2);
//     Dune::CO2HeatFlux<GridType, NumberType> problem(law, 1.e7); 
      Dune::BrooksCoreyLaw law(brine, co2);
      Dune::CO2HeatFlux<GridType, NumberType> problem(law, 1.e7, 0.2, 0.05, 10000, 2); 

    typedef Dune::BoxPwSnTe<GridType, NumberType> TwoPhase;
    TwoPhase twoPhase(grid, problem);


    
    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "co2", 1);
    
    Dune::Timer timer;
    timer.reset();
    timeloop.execute(twoPhase);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
     
	//std::cout << twoPhase.injected() << " kg CO2 injected." << std::endl;

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
#else 

int main (int argc , char **argv) try
{
  std::cout << "Please install the UG library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
