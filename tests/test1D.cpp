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


template <std::size_t ndim, 
		  typename CellT,
		  typename SolutionIterator,
		  class ConservedQuantityFunctor,
		  class FaceFluxFunctor,
		  class SourceFunctors...>
struct ExplicitUpdate{
	// SolutionIterator msit;
	double mdt, mdx, mdtdx;
	std::tuple<ConservedQuantityFunctor, FaceFluxFunctor, SourceFunctors...> mtpl;

	constexpr ExplicitUpdate(double dt, double dx, std::tuple<ConservedQuantityFunctor, FaceFluxFunctor, SourceFunctors...> tpl)
	: mdt(dt)
	, mdx(dx)
	, mdtdx(dt/dx)
	, mtpl(tpl){};

	double operator()(CellT & cl){
		double flx=0.0;
		for (auto d=0; d<ndim; d++){
			flx += 0.5*(std::get<1>(mtpl).operator()(cl.getNeighborMin<d>()));
			flx -= 0.5*(std::get<1>(mtpl).operator()(cl.getNeighborMax<d>()));
		}
		return std::get<0>(mtpl).operator()(cl) - mdtdx*flx;
	}
};


template <std::size_t ndim,
		  typename CellT,
		  typename SolutionIterator,
		  class ConservedQuantityFunctor,
		  class FaceFluxFunctor,
		  class SourceFunctors...>
ExplicitUpdate<ndim,CellT,SolutionIterator,
			   ConservedQuantityFunctor,
			   FaceFluxFunctor,
			   SourceFunctors...> make_explicit_update(double dt, double dx, ConservedQuantityFunctor cf, FaceFluxFunctor ff, SourceFunctors... sf){
	return ExplicitUpdate(dt, dx, std::make_tuple(cf,ff,sf));
};


// create an update struct
template <typename CellIterator, 
		  typename SolutionIterator, 
		  typename Functor> 
void for_each_update(CellIterator first, CellIterator last, Functor f, SolutionIterator sfirst){
	auto soln = sfirst;
	for (auto it=first; it!=last; it++){
		(*soln) = f.operator()(it);
		soln++;
	}
}


int main(int argc, char * argv[]){

	std::size_t ncells = 50;

	// time-stepping parameters
	const double c 	= 1.0;								// velocity [m/s]
	const double cfl = 0.9;								// cfl number
	const double dx = 1.0/static_cast<double>(ncells);	// cell size [m]
	const double dt = cfl*dx/c;							// time step [s]

	// std::cout << typeid(y).name() << std::endl;

	// YeeCell typedefs
	typedef fidi::ObjectStencil<1,1,SomeQuantity> 		CellT;

	// define the domain data structure and connect neighbors
	std::vector<CellT> cells(50);
	for (auto i=1; i<49; i++){
		cells[i].setNeighborMin(0, cells[i-1]);
		cells[i].setNeighborMax(0, cells[i+1]);
	}

	// define the initial density
	for (auto i=20; i<30; i++) cells[i].n() = 1.0;

	// define the conserved quantities 
	// std::tuple<std::function<const double(CellT &)>, std::function<const double(CellT &)>> 	conserved = {, };

	// initialize a solution vector
	std::vector<double> soln(ncells, 0.0);

	// make the explicit update struct
	auto exp_up = make_explicit_update(dt, dx, [](CellT & cl){return cl.n();}, [&c](CellT & cl){return cl.n()*c;});
	
	// start time-stepping
	for (auto t=0; t<200; t++){
		for_each_update(++cells.begin(), --cells.end(), exp_up, ++soln.begin());
		for_each_update(++soln.begin(), --soln.end(), [](double & d){}, ++cells.begin());

		// std::cout << "n:" ;
		for(auto i=0; i<50; i++) std::cout << ", " << cells[i].n() ;
		std::cout << std::endl;
	}

	return 0;
}