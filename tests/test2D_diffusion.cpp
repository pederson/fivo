#include <fivo.h>

#include <iostream>
#include <typeinfo>
#include <vector>
#include <algorithm>
#include <functional>

#include "../../fidi/fidi.h"


using namespace fivo;

// compile with:
// 			clang++ -std=c++14 -I../ test2D_diffusion.cpp -o test2D_diff









struct SomeQuantity{
	double mN;

	SomeQuantity():mN(0.0){};

	double & n(){return mN;};
	constexpr double n() const {return mN;};
	constexpr double D() const {return 1.0;};
};


template <typename T> int sgn(T val) {
    return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
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


///////////////// TWO-CELL FLUX FUNCTORS /////////////////////////
template <class InwardVelocityFunctor, 
		  class ConservedQuantityFunctor,
		  class FluxFunctor>
struct UpwindTwoCell{
	double 							mdxdt;
	InwardVelocityFunctor 			mIVF;
	ConservedQuantityFunctor 		mCQF;
	FluxFunctor 					mFF;

	template <typename CellT>
	double operator()(const CellT & c, const CellT & n){
		double ubar = mIVF.operator()(c,n);
		return std::max(ubar, 0.0)*(mFF.operator()(c)-mFF.operator()(n));
	}
};


template <class InwardVelocityFunctor, 
		  class ConservedQuantityFunctor,
		  class FluxFunctor>
struct LaxFriedrichsTwoCell{
	double 							mdxdt;
	InwardVelocityFunctor 			mIVF;
	ConservedQuantityFunctor 		mCQF;
	FluxFunctor 					mFF;

	template <typename CellT>
	double operator()(const CellT & c, const CellT & n){
		double ubar = mIVF.operator()(c,n);
		return 0.5*(mFF.operator()(c)+mFF.operator()(n)) - 0.5*mdxdt*sgn(ubar)*(mCQF.operator()(c)-mCQF.operator()(n));
	}
};



inline double LaxFriedrichsTwoCellFlux(double dt, double dx, double ubar,
									   double fc, double fn, double qc, double qn){
	return 0.5*(fc+fn)-0.5*dx/dt*sgn(ubar)*(qc-qn);
}


template <class InwardVelocityFunctor, 
		  class ConservedQuantityFunctor,
		  class FluxFunctor>
struct LaxWendroffTwoCell{
	double 							mdxdt;
	InwardVelocityFunctor 			mIVF;
	ConservedQuantityFunctor 		mCQF;
	FluxFunctor 					mFF;

	LaxWendroffTwoCell(double dt, double dx, InwardVelocityFunctor ivf, ConservedQuantityFunctor cf, FluxFunctor ff)
	: mdxdt(dx/dt)
	, mIVF(ivf)
	, mCQF(cf)
	, mFF(ff){};

	template <typename CellT>
	double operator()(const CellT & c, const CellT & n){
		double ubar = mIVF.operator()(c,n);
		return 0.5*(mFF.operator()(c)+mFF.operator()(n)) - 0.5/mdxdt*ubar*ubar*sgn(ubar)*(mCQF.operator()(c)-mCQF.operator()(n));
	}
};

template <class InwardVelocityFunctor, 
		  class ConservedQuantityFunctor,
		  class FluxFunctor>
auto make_lax_wendroff_two_cell(double dt, double dx, InwardVelocityFunctor ivf, ConservedQuantityFunctor cf, FluxFunctor ff){
	return LaxWendroffTwoCell<InwardVelocityFunctor, ConservedQuantityFunctor, FluxFunctor>(dt, dx, ivf, cf, ff);
}
///////////////////////////////////////////////////////////////////






template <class ConservedQuantityFunctor,
		  class FluxFunctor>
struct LaxWendroff{
	double mdx, mdtdx;
	ConservedQuantityFunctor mCQF;
	FluxFunctor 			 mFF;

	constexpr LaxWendroff(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff)
	: mdx(dx)
	, mdtdx(dt/dx)
	, mCQF(cf)
	, mFF(ff) {};

	template <typename CellT>
	double operator()(const CellT & cl){
		double flx = 0.0;

		double nr = 0.5*(mCQF.operator()(cl) + mCQF.operator()(cl.getNeighborMax(0))) 
				  - mdtdx*0.5*(mFF.operator()(cl.getNeighborMax(0)) - mFF.operator()(cl));
		double nl = 0.5*(mCQF.operator()(cl) + mCQF.operator()(cl.getNeighborMin(0))) 
				  - mdtdx*0.5*(mFF.operator()(cl) - mFF.operator()(cl.getNeighborMin(0)));

		// double nu = 0.5*(mCQF.operator()(cl) + mCQF.operator()(cl.getNeighborMax(1))) 
		// 		  - mdtdx*0.5*(mFF.operator()(cl.getNeighborMax(1)) - mFF.operator()(cl));
		// double nd = 0.5*(std::get<0>(mtpl).operator()(cl) + std::get<0>(mtpl).operator()(cl.getNeighborMin(1))) 
		// 		  - mdtdx*0.5*(std::get<1>(mtpl).operator()(cl) - std::get<1>(mtpl).operator()(cl.getNeighborMin(1)));

		// flux in
		flx -= mFF.operator()(nl);
		// flx -= mFF.operator()(nd);

		// flux out
		flx += mFF.operator()(nr);
		// flx += mFF.operator()(nu);

		return flx/mdx;
	}
};


template <class ConservedQuantityFunctor,
		  class FluxFunctor>
auto make_lax_wendroff(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff){
	return LaxWendroff<ConservedQuantityFunctor,
			   			  FluxFunctor>(dt, dx, cf, ff);
};








template <class ConservedQuantityFunctor,
		  class DiffusionCoefficientFunctor>
struct DiffusiveFlux{
	double mdx, mdtdx;
	ConservedQuantityFunctor 				mCQF;
	DiffusionCoefficientFunctor 			mDCF;

	constexpr DiffusiveFlux(double dt, double dx, ConservedQuantityFunctor cf, DiffusionCoefficientFunctor dcf)
	: mdx(dx)
	, mdtdx(dt/dx)
	, mCQF(cf)
	, mDCF(dcf) {};

	template <typename CellT>
	double operator()(const CellT & cl){
		double flx = 0.0;

		double diff = mDCF.operator()(cl);
		double nr = mdtdx*diff*(mCQF.operator()(cl.getNeighborMax(0)) - mCQF.operator()(cl));
		double nl = mdtdx*diff*(mCQF.operator()(cl.getNeighborMin(0)) - mCQF.operator()(cl));
		double nu = mdtdx*diff*(mCQF.operator()(cl.getNeighborMax(1)) - mCQF.operator()(cl));
		double nd = mdtdx*diff*(mCQF.operator()(cl.getNeighborMin(1)) - mCQF.operator()(cl));

		// flux out
		flx += nr+nl+nu+nd;

		return flx/mdx;
	}
};


template <class ConservedQuantityFunctor,
		  class DiffusionCoefficientFunctor>
auto make_diffusive_flux(double dt, double dx, ConservedQuantityFunctor cf, DiffusionCoefficientFunctor ff){
	return DiffusiveFlux<ConservedQuantityFunctor,
			   			  DiffusionCoefficientFunctor>(dt, dx, cf, ff);
};






















int main(int argc, char * argv[]){

	std::size_t ncells = 100;

	// time-stepping parameters
	const double cx 	= 1.0/sqrt(2.0);								// velocity [m/s]
	const double cy 	= 1.0/sqrt(2.0);
	const double cfl = 0.5;								// cfl number
	const double dx = 1.0/static_cast<double>(ncells);	// cell size [m]
	const double dt = cfl*dx/std::fabs(sqrt(cx*cx+cy*cy));							// time step [s]

	// std::cout << typeid(y).name() << std::endl;

	// YeeCell typedefs
	typedef fidi::ObjectStencil<2,1,SomeQuantity> 		CellT;

	// define the domain data structure and connect neighbors
	std::vector<std::vector<CellT>> cells(ncells);
	for (auto i=0; i<ncells; i++) cells[i].resize(ncells);
	for (auto i=1; i<ncells-1; i++){
		for (auto j=1; j<ncells-1; j++){
		cells[i][j].setNeighborMin(cells[i-1][j], 0);
		cells[i][j].setNeighborMax(cells[i+1][j], 0);
		cells[i][j].setNeighborMin(cells[i][j-1], 1);
		cells[i][j].setNeighborMax(cells[i][j+1], 1);
		}
	}

	// define the initial density
	for (auto i=0; i<ncells; i++){
		for (auto j=0; j<ncells; j++){
			cells[i][j].n() = exp(-static_cast<double>((i-ncells/2)*(i-ncells/2) + (j-ncells/2)*(j-ncells/2))/(2.0*5.0*5.0));
		}
	} 

	// for(auto i=0; i<ncells; i++){
	// 	for (auto j=0; j<ncells; j++){
	// 		std::cout << ", " << cells[i][j].n() ;
	// 	}
	// } 
	// std::cout << std::endl;

	// initialize a solution vector
	std::vector<std::vector<double>> soln(ncells);
	for (auto i=0; i<ncells; i++) soln[i].resize(ncells);

	// make the explicit update struct
	auto cons = [](const CellT & cl){return cl.n();};
	auto flux = [&cx](const CellT & cl){return cx*cl.n();};

	struct FluxGetter{
		double mCx, mCy;

		FluxGetter(double c): mCx(c){};
		double operator()(const CellT & cl){return mCx*cl.n();};
		double operator()(double n){return mCx*n;};
	};


	auto diffusive = make_diffusive_flux(dt, dx, cons, [](auto & c){return c.D();});
	auto exp_up = make_explicit_update(dt, cons, diffusive);
	
	// start time-stepping
	for (auto t=0; t<50; t++){
		for(auto i=1; i<ncells-1; i++){
			for (auto j=1; j<ncells-1; j++){
				soln[i][j] = exp_up.operator()(cells[i][j]);
			}
		} 
		for(auto i=0; i<ncells; i++){
			for (auto j=0; j<ncells; j++){
				cells[i][j].n() = soln[i][j];
			}
		} 

		// std::cout << "n:" ;
		for(auto i=0; i<ncells; i++){
			for (auto j=0; j<ncells; j++){
				std::cout << ", " << cells[i][j].n() ;
			}
		} 
		std::cout << std::endl;
	}

	return 0;
}