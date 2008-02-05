#include "config.h"
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "fvca5test7problem.hh"
#include "../benchmarkresult.hh"

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::UGGrid<dim> GridType; 

    if (argc != 2 && argc != 3) {
    	std::cout << "Usage: fvca5_test7 dgffilename [refinementsteps]" << std::endl;
    	return (1);
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( argv[1] );

    // grid reference 
    GridType& grid = *gridPtr;

    if (refinementSteps)
    	grid.globalRefine(refinementSteps);

    double delta = 0.2;
	Dune::FVCA5Test7Problem<GridType, NumberType> problem(delta);
    Dune::SimpleProblem<GridType, NumberType> satprob;

    Dune::Timer timer;
    timer.reset();
    //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    //Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem);
    Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem, satprob, grid.maxLevel());
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    //printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    
    Dune::BenchmarkResult result;
    result.evaluate(grid, problem, diffusion);  
    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.setf(std::ios_base::uppercase);
    std::cout.precision(2);
    
    std::cout << "sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf \n        = " 
    	<< result.flux0 << " + " << result.flux1 << " + " 
    	<< result.fluy0 << " + " << result.fluy1 << " - " 
    	<< result.sumf << "\n        = " << result.sumflux << std::endl;
    std::cout << "relative discrete L2 error erl2 = " << result.relativeL2Error << std::endl;
    std::cout << "relative discrete L2 error ergrad = " << result.ergrad << std::endl;
    std::cout << "errflx0 = abs((flux0 + exactflux0)/exactflux0) = abs((" 
    	<< result.flux0 << " + " << result.exactflux0 << ")/" << result.exactflux0 
    	<< ") = " << result.errflx0 << std::endl;
    std::cout << "errflx1 = abs((flux1 + exactflux1)/exactflux1) = abs((" 
    	<< result.flux1 << " + " << result.exactflux1 << ")/" << result.exactflux1 
    	<< ") = " << result.errflx1 << std::endl;
    std::cout << "errfly0 = abs((fluy0 + exactfluy0)/exactfluy0) = abs((" 
    	<< result.fluy0 << " + " << result.exactfluy0 << ")/" << result.exactfluy0 
    	<< ") = " << result.errfly0 << std::endl;
    std::cout << "errfly1 = abs((fluy1 + exactfluy1)/exactfluy1) = abs((" 
    	<< result.fluy1 << " + " << result.exactfluy1 << ")/" << result.exactfluy1 
    	<< ") = " << result.errfly1 << std::endl;
    std::cout << "mean value error erflm = " << result.erflm << std::endl;
    std::cout << "energy ener2 = " << result.ener2 << std::endl;
    std::cout << "umin = " << result.uMin << std::endl;
    std::cout << "umax = " << result.uMax << std::endl;
    
    
    diffusion.vtkout("fvca5_test7", 0);


    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
