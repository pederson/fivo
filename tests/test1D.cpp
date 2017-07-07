#include <fivo.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>
#include <functional>

#include "../../fidi/fidi.h"


// using namespace fivo;

// compile with:
// 			clang++ -std=c++14 -I../ test1D.cpp -o test

struct SomeQuantity{
	double mN;

	SomeQuantity():mN(0.0){};

	double & n(){return mN;};
	constexpr double n() const {return mN;};
};


template <class ConservedQuantityFunctor,
		  class NumericalFluxFunctor,
		  class... SourceFunctors>
struct ExplicitUpdate{
	double mdt, mdx, mdtdx;
	std::tuple<ConservedQuantityFunctor, NumericalFluxFunctor, SourceFunctors...> mtpl;

	constexpr ExplicitUpdate(double dt, double dx, std::tuple<ConservedQuantityFunctor, NumericalFluxFunctor, SourceFunctors...> tpl)
	: mdt(dt)
	, mdx(dx)
	, mdtdx(dt/dx)
	, mtpl(tpl){};

	template <typename CellT>
	double operator()(const CellT & cl){
		double flx=0.0;

		flx -= std::get<1>(mtpl).operator()(cl, cl.getNeighborMin(0));
		flx += std::get<1>(mtpl).operator()(cl.getNeighborMax(0), cl);

		return std::get<0>(mtpl).operator()(cl) - mdtdx*flx;
	}
};


template <class ConservedQuantityFunctor,
		  class NumericalFluxFunctor,
		  class... SourceFunctors>
auto make_explicit_update(double dt, double dx, ConservedQuantityFunctor cf, NumericalFluxFunctor ff, SourceFunctors... sf){
	return ExplicitUpdate<ConservedQuantityFunctor,
			   			  NumericalFluxFunctor,
			   			  SourceFunctors...>(dt, dx, std::make_tuple(cf,ff,sf...));
};


template <class ConservedQuantityFunctor, 
		  class FluxFunctor>
struct NumericalFlux{
	ConservedQuantityFunctor mcf;
	FluxFunctor mff;
	double mCc; 	// coefficient of the flux on the conserved quantity value
	double mCn;		// coefficient of the flux on the neighbor conserved quantity value
	double mFc;		// coefficient of the flux on the flux value
	double mFn;		// coefficient of the flux on the neighbor flux value

	NumericalFlux(ConservedQuantityFunctor cf, FluxFunctor ff, double Cc, double Cn, double Fc, double Fn)
	: mcf(cf), mff(ff)
	, mCc(Cc), mCn(Cn)
	, mFc(Fc), mFn(Fn) {};

	template <typename CellT>
	double operator()(const CellT & cl, const CellT & neighb){
		return mCc*mcf.operator()(cl)+mCn*mcf.operator()(neighb)+mFc*mff.operator()(cl)+mFn*mff.operator()(neighb);
	}
};

template <class ConservedQuantityFunctor,
		  class FluxFunctor>
auto make_numerical_flux(ConservedQuantityFunctor cf, FluxFunctor ff, double Cc, double Cn, double Fc, double Fn){
	return NumericalFlux<ConservedQuantityFunctor, FluxFunctor>(cf, ff, Cc, Cn, Fc, Fn);
}


// create an update struct
template <typename CellIterator, 
		  typename SolutionIterator, 
		  typename Functor> 
void for_each_update(CellIterator first, CellIterator last, Functor f, SolutionIterator sfirst){
	auto soln = sfirst;
	for (auto it=first; it!=last; it++){
		(*soln) = f.operator()(*it);
		soln++;
	}
}


int main(int argc, char * argv[]){

	std::size_t ncells = 50;

	// time-stepping parameters
	const double c 	= 1.0;								// velocity [m/s]
	const double cfl = 0.5;								// cfl number
	const double dx = 1.0/static_cast<double>(ncells);	// cell size [m]
	const double dt = cfl*dx/std::fabs(c);							// time step [s]

	// std::cout << typeid(y).name() << std::endl;

	// YeeCell typedefs
	typedef fidi::ObjectStencil<1,1,SomeQuantity> 		CellT;

	// define the domain data structure and connect neighbors
	std::vector<CellT> cells(ncells);
	for (auto i=1; i<ncells-1; i++){
		cells[i].setNeighborMin(cells[i-1], 0);
		cells[i].setNeighborMax(cells[i+1], 0);
	}

	// define the initial density
	for (auto i=20; i<30; i++) cells[i].n() = 1.0-std::fabs(i-25)/5.0;

	// initialize a solution vector
	std::vector<double> soln(ncells, 0.0);

	// make the explicit update struct
	auto cons = [](const CellT & cl){return cl.n();};
	auto flux = [&c](const CellT & cl){return c*cl.n();};

	auto lax_fried = make_numerical_flux(cons, flux, -0.5*dx/dt, 0.5*dx/dt, 0.5, 0.5);
	auto exp_up_LF = make_explicit_update(dt, dx, cons, lax_fried);

	auto lax_wend = [&c, &exp_up_LF](const CellT & cl, const CellT & neighb){return c*exp_up_LF.operator()(cl);};
	auto exp_up = make_explicit_update(dt, dx, cons, lax_wend);
	
	
	// start time-stepping
	for (auto t=0; t<20; t++){
		for_each_update(++cells.begin(), --cells.end(), exp_up, ++soln.begin());
		for(auto i=0; i<ncells; i++) cells[i].n() = soln[i];
		// for_each_update(++soln.begin(), --soln.end(), [](double & d){}, ++cells.begin());

		// std::cout << "n:" ;
		for(auto i=0; i<ncells; i++) std::cout << ", " << cells[i].n() ;
		std::cout << std::endl;
	}

	return 0;
}