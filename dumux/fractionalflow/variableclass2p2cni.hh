// $Id$

#ifndef DUNE_VARIABLECLASS2P2CNI_HH
#define DUNE_VARIABLECLASS2P2CNI_HH

/**
 * @file
 * @brief  class including the variables needed for decoupled nonisothermal 2p2c computations.
 * @author Jochen Fritz
 *
 * The function of this class is to concentrate all variables which are needed for a decoupled
 * two phase two component calculation in one place and to do the output.
 */

namespace Dune {

template<class G, class RT> class VariableClass2p2cni {

	enum {n=G::dimension};
	typedef typename G::ctype DT;
	typedef Dune::BlockVector< Dune::FieldVector<RT,1> > ScalarType;
	typedef Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, n>, 2*n> >
			VelType;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	ScalarType saturation;
	ScalarType pressure;
	VelType velocity;
	ScalarType totalConcentration;
	ScalarType wet_X1, nonwet_X1;
	ScalarType density_wet, density_nonwet;
	ScalarType temperature, enthalpy_l, enthalpy_g;
	ScalarType volErr;

	template<int dim> struct ElementLayout
	{
		bool contains(Dune::GeometryType gt)
		{
			return gt.dim() == dim;
		}
	};

	typedef typename G::LevelGridView::IndexSet IS;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>	ElementMapper;
	G& grid;
	ElementMapper mapper;

private:
	int size;
public:

	VariableClass2p2cni(G& g, int lev = 0) :
		grid(g), mapper(g, g.levelIndexSet(lev)), size(mapper.size())
	{
		saturation.resize(size);
		pressure.resize(size);
		velocity.resize(size);
		totalConcentration.resize(3*size);
		wet_X1.resize(size);
		density_wet.resize(size);
		nonwet_X1.resize(size);
		density_nonwet.resize(size);
		volErr.resize(size);
		temperature.resize(size);
		enthalpy_l.resize(size);
		enthalpy_g.resize(size);

		saturation = 0;
		pressure = 0;
		totalConcentration = 0;
		wet_X1 = 0;
		density_wet = 0;
		nonwet_X1 = 0;
		density_nonwet = 0;
		volErr = 0;
	}

	ScalarType& sat() const {
		return saturation;
	}
	ScalarType& press() const {
		return pressure;
	}
	VelType& vel() const {
		return velocity;
	}

	const Dune::FieldVector<RT,1>& sat(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return saturation[mapper.map(e)];;
	}

	const Dune::FieldVector<RT,1>& press(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return pressure[mapper.map(e)];
	}

	const Dune::FieldVector<DT,n>& vTotal(const Entity& e,
			const int numberInSelf) const {
		int elemId = mapper.map(e);

		return (velocity[elemId][numberInSelf]);
	}

	/*! @brief writes all variables to a VTK File
	 *
	 *  The file name is "<name>-<k>.vtu" where k is an integer number.
	 *  @param name specifies the name of the VTK file
	 *  @param k specifies a number
	 */
	void vtkout(const char* name, int k, int pressurelevel = 0, int satlevel = 0)
	{
		if (pressurelevel == satlevel)
		{
			ScalarType C1, C2, U;
			C1.resize(size); C2.resize(size); U.resize(size);
			for (int i = 0; i < size; i++)
			{
				C1[i] = totalConcentration[i];
				C2[i] = totalConcentration[i + size];
				U[i] = totalConcentration[i + 2*size];
			}
			volErr *= 100.;
			VTKWriter<G> vtkwriter(grid);
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriter.addCellData(saturation, "saturation [-]");
			vtkwriter.addCellData(pressure, "total pressure p~ [Pa]");
			vtkwriter.addCellData(temperature, "temperature [K]");
			vtkwriter.addCellData(C1, "total concentration 1 [kg/m^3]");
			vtkwriter.addCellData(C2, "total concentration 2 [kg/m^3]");
			vtkwriter.addCellData(wet_X1, "Mass fraction 1 in wetting phase [kg/m^3]");
			vtkwriter.addCellData(nonwet_X1, "Mass fraction 1 in non-wetting phase [kg/m^3]");
			vtkwriter.addCellData(density_wet, "Density of wetting phase [kg/m^3]");
			vtkwriter.addCellData(density_nonwet, "Density of non-wetting phase [kg/m^3]");
			vtkwriter.addCellData(U, "Internal energy in cell [J/m^3]");
			vtkwriter.addCellData(enthalpy_l, "hl [J/kg]");
			vtkwriter.addCellData(enthalpy_g, "hg [J/kg]");
			vtkwriter.addCellData(volErr, "volume error [%]");
			vtkwriter.write(fname, VTKOptions::ascii);
			volErr /= 100.;
		}
		else
		{
			Dune::VTKWriter<G, typename G::LevelGridView>
					vtkwriterpressure(grid.levelView(pressurelevel));
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriterpressure.addCellData(pressure, "total pressure p~");
			vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

			VTKWriter<G, typename G::LevelGridView>
					vtkwritersaturation(grid.levelView(satlevel));
			sprintf(fname, "%s-press%05d", name, k);
			vtkwritersaturation.addCellData(saturation, "saturation");
			vtkwritersaturation.write(fname, VTKOptions::ascii);
		}
		return;
	}
	void vtkoutpressure(const char* name, int k) const {
		VTKWriter<G> vtkwriter(grid);
		char fname[128];
		sprintf(fname, "%s-press%05d", name, k);
		vtkwriter.addCellData(pressure, "total pressure p~");
		vtkwriter.write(fname, VTKOptions::ascii);
	}
};
}
#endif
