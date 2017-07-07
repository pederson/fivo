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



namespace Detail{
	// compile-time for_each on a tuple
	template <typename TupleType, typename FunctionType>
	void for_each(TupleType&&, FunctionType
	            , std::integral_constant<std::size_t, std::tuple_size<typename std::remove_reference<TupleType>::type >::value>) {}

	template <std::size_t I, typename TupleType, typename FunctionType
	       , typename = typename std::enable_if<I!=std::tuple_size<typename std::remove_reference<TupleType>::type>::value>::type >
	void for_each(TupleType&& t, FunctionType f, std::integral_constant<size_t, I>)
	{
	    f(std::get<I>(t));
	    for_each(std::forward<TupleType>(t), f, std::integral_constant<size_t, I + 1>());
	}

	template <typename TupleType, typename FunctionType>
	void for_each(TupleType&& t, FunctionType f)
	{
	    for_each(std::forward<TupleType>(t), f, std::integral_constant<size_t, 0>());
	}
}






struct SomeQuantity{
	double mN;

	SomeQuantity():mN(0.0){};

	double & n(){return mN;};
	constexpr double n() const {return mN;};
};



template <class ConservedQuantityFunctor,
		  class NumericalFluxInFunctor,
		  class NumericalFluxOutFunctor,
		  class... SourceFunctors>
struct ExplicitUpdate{
	double mdt, mdx, mdtdx;
	ConservedQuantityFunctor mcf;
	NumericalFluxInFunctor mfif;
	NumericalFluxOutFunctor mfof;
	std::tuple<SourceFunctors...> mstpl;

	constexpr ExplicitUpdate(double dt, double dx, ConservedQuantityFunctor cf, NumericalFluxInFunctor fif, NumericalFluxOutFunctor fof, std::tuple<SourceFunctors...> stpl)
	: mdt(dt)
	, mdx(dx)
	, mdtdx(dt/dx)
	, mcf(cf)
	, mfif(fif)
	, mfof(fof)
	, mstpl(stpl){};

	template <typename CellT>
	double operator()(const CellT & cl){
		double flx=0.0;
		// flux in
		flx += mfif.operator()(cl);
		// flux out	
		flx -= mfof.operator()(cl);

		// sources
		double src=0.0;
		Detail::for_each(mstpl, [&cl, &src](const auto & s){src += s.operator()(cl);});

		return mcf.operator()(cl) + mdtdx*flx + mdt*src;
	}
};


template <class ConservedQuantityFunctor,
		  class NumericalFluxInFunctor,
		  class NumericalFluxOutFunctor,
		  class... SourceFunctors>
auto make_explicit_update(double dt, double dx, ConservedQuantityFunctor cf, NumericalFluxInFunctor fif, NumericalFluxOutFunctor fof, SourceFunctors... sf){
	return ExplicitUpdate<ConservedQuantityFunctor,
						  NumericalFluxInFunctor,
						  NumericalFluxOutFunctor,
						  SourceFunctors...>(dt, dx, cf, fif, fof, std::make_tuple(sf...));
}



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




template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class NeighborGetter>
struct LaxFriedrichsFluxIn : public NumericalFlux<ConservedQuantityFunctor, FluxFunctor>{
	typedef NumericalFlux<ConservedQuantityFunctor, FluxFunctor> 		NumFlux;
	NeighborGetter mng;

	LaxFriedrichsFluxIn(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff, NeighborGetter ng)
	: NumericalFlux<ConservedQuantityFunctor, FluxFunctor>(cf, ff, 0.5, 0.5, -0.5*dt/dx, 0.5*dt/dx)
	, mng(ng) {};

	template <typename CellT>
	double operator()(const CellT & cl){return NumFlux::operator()(cl, mng.operator()(cl));};
};


template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class NeighborGetter>
auto make_lax_friedrichs_in(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff, NeighborGetter ng){
	return LaxFriedrichsFluxIn<ConservedQuantityFunctor, FluxFunctor, NeighborGetter>(dt, dx, cf, ff, ng);
}



template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class NeighborGetter>
struct LaxFriedrichsFluxOut : public NumericalFlux<ConservedQuantityFunctor, FluxFunctor>{
	typedef NumericalFlux<ConservedQuantityFunctor, FluxFunctor> 		NumFlux;
	NeighborGetter mng;

	LaxFriedrichsFluxOut(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff, NeighborGetter ng)
	: NumericalFlux<ConservedQuantityFunctor, FluxFunctor>(cf, ff, 0.5, 0.5, 0.5*dt/dx, -0.5*dt/dx)
	, mng(ng) {};

	template <typename CellT>
	double operator()(const CellT & cl){return NumFlux::operator()(cl, mng.operator()(cl));};
};

template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class NeighborGetter>
auto make_lax_friedrichs_out(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff, NeighborGetter ng){
	return LaxFriedrichsFluxOut<ConservedQuantityFunctor, FluxFunctor, NeighborGetter>(dt, dx, cf, ff, ng);
}



// template <class ConservedQuantityFunctor,
// 		  class NumericalFluxFunctor,
// 		  class... SourceFunctors>
// struct ExplicitUpdate{
// 	double mdt, mdx, mdtdx;
// 	std::tuple<ConservedQuantityFunctor, NumericalFluxFunctor> mtpl;
// 	std::tuple<SourceFunctors...> mstpl;

// 	constexpr ExplicitUpdate(double dt, double dx, std::tuple<ConservedQuantityFunctor, NumericalFluxFunctor> tpl, std::tuple<SourceFunctors...> stpl)
// 	: mdt(dt)
// 	, mdx(dx)
// 	, mdtdx(dt/dx)
// 	, mtpl(tpl)
// 	, mstpl(stpl){};

// 	template <typename CellT>
// 	double operator()(const CellT & cl){
// 		double flx=0.0;

// 		flx -= std::get<1>(mtpl).operator()(cl, cl.getNeighborMin(0));
// 		flx += std::get<1>(mtpl).operator()(cl.getNeighborMax(0), cl);


// 		double src=0.0;
// 		Detail::for_each(mstpl, [&cl, &src](const auto & s){src += s.operator()(cl);});
// 		// for (auto s=0; s<sizeof...(SourceFunctors); s++) src += std::get<2+s>(mtpl).operator()(cl);

// 		return std::get<0>(mtpl).operator()(cl) - mdtdx*flx + mdt*src;
// 	}
// };


// template <class ConservedQuantityFunctor,
// 		  class NumericalFluxFunctor,
// 		  class... SourceFunctors>
// auto make_explicit_update(double dt, double dx, ConservedQuantityFunctor cf, NumericalFluxFunctor ff, SourceFunctors... sf){
// 	return ExplicitUpdate<ConservedQuantityFunctor,
// 			   			  NumericalFluxFunctor,
// 			   			  SourceFunctors...>(dt, dx, std::make_tuple(cf,ff), std::make_tuple(sf...));
// };





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









template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class... SourceFunctors>
struct LaxWendroff{
	double mdt, mdx, mdtdx;
	std::tuple<ConservedQuantityFunctor, FluxFunctor> mtpl;
	std::tuple<SourceFunctors...> mstpl;

	constexpr LaxWendroff(double dt, double dx, std::tuple<ConservedQuantityFunctor, FluxFunctor> tpl, std::tuple<SourceFunctors...> stpl)
	: mdt(dt)
	, mdx(dx)
	, mdtdx(dt/dx)
	, mtpl(tpl)
	, mstpl(stpl){};

	template <typename CellT>
	double operator()(const CellT & cl){
		double flx=0.0;

		double nr = 0.5*(std::get<0>(mtpl).operator()(cl) + std::get<0>(mtpl).operator()(cl.getNeighborMax(0))) 
				  - mdtdx*0.5*(std::get<1>(mtpl).operator()(cl.getNeighborMax(0)) - std::get<1>(mtpl).operator()(cl));
		double nl = 0.5*(std::get<0>(mtpl).operator()(cl) + std::get<0>(mtpl).operator()(cl.getNeighborMin(0))) 
				  - mdtdx*0.5*(std::get<1>(mtpl).operator()(cl) - std::get<1>(mtpl).operator()(cl.getNeighborMin(0)));

		// flux in
		flx += std::get<1>(mtpl).operator()(nl);

		// flux out
		flx -= std::get<1>(mtpl).operator()(nr);

		// flx -= std::get<1>(mtpl).operator()(cl, cl.getNeighborMin(0));
		// flx += std::get<1>(mtpl).operator()(cl.getNeighborMax(0), cl);

		double src=0.0;
		Detail::for_each(mstpl, [&cl, &src](const auto & s){src += s.operator()(cl);});

		return std::get<0>(mtpl).operator()(cl) + mdtdx*flx + mdt*src;
	}
};



template <class ConservedQuantityFunctor,
		  class FluxFunctor,
		  class... SourceFunctors>
auto make_lax_wendroff(double dt, double dx, ConservedQuantityFunctor cf, FluxFunctor ff, SourceFunctors... sf){
	return LaxWendroff<ConservedQuantityFunctor,
			   			  FluxFunctor,
			   			  SourceFunctors...>(dt, dx, std::make_tuple(cf,ff), std::make_tuple(sf...));
};


int main(int argc, char * argv[]){

	std::size_t ncells = 100;

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

	struct FluxGetter{
		double mC;

		FluxGetter(double c): mC(c){};
		double operator()(const CellT & cl){return mC*cl.n();};
		double operator()(double n){return mC*n;};
	};

	auto lax_fried = make_explicit_update(dt, dx, cons, 
										  make_lax_friedrichs_in( dt, dx, cons, FluxGetter(c), [](const CellT & cl){return cl.getNeighborMin(0);}),
										  make_lax_friedrichs_out(dt, dx, cons, FluxGetter(c), [](const CellT & cl){return cl.getNeighborMax(0);}));

	// auto lax_fried = make_numerical_flux(cons, flux, -0.5*dx/dt, 0.5*dx/dt, 0.5, 0.5);
	// auto exp_up = make_explicit_update(dt, dx, cons, lax_fried);

	// auto lax = make_numerical_flux(cons, flux, 0.5, 0.5, -0.5*dt/dx, 0.5*dt/dx);
	// auto exp_up = make_explicit_update(dt, dx, cons, [&c, &lax](const CellT & cl, const CellT & neighb){return c*lax.operator()(cl, neighb);});
	// auto lax_wend = [&c, &exp_up_LF](const CellT & cl, const CellT & neighb){return c*exp_up_LF.operator()(cl);};
	// auto lax_wend = [&c, &exp_up_LF](const CellT & cl, const CellT & neighb){return c*exp_up_LF.operator()(cl);};
	// auto exp_up = make_explicit_update(dt, dx, cons, lax_wend);
	
	auto exp_up = make_lax_wendroff(dt, dx, cons, FluxGetter(c));
	
	// start time-stepping
	for (auto t=0; t<50; t++){
		for_each_update(++cells.begin(), --cells.end(), exp_up, ++soln.begin());
		for(auto i=0; i<ncells; i++) cells[i].n() = soln[i];
		// for_each_update(++soln.begin(), --soln.end(), [](double & d){}, ++cells.begin());

		// std::cout << "n:" ;
		for(auto i=0; i<ncells; i++) std::cout << ", " << cells[i].n() ;
		std::cout << std::endl;
	}

	return 0;
}