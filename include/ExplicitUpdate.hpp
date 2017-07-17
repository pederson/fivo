#ifndef _EXPLICITUPDATE_H
#define _EXPLICITUPDATE_H

#include <type_traits>
#include <tuple>


namespace fivo{



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
} // end namespace Detail





template <class ConservedQuantityFunctor,
		  class NumericalFluxInFunctor,
		  class... SourceFunctors>
struct ExplicitUpdate{
	double mdt;
	ConservedQuantityFunctor 			mCQF;
	NumericalFluxInFunctor 				mNFF;
	std::tuple<SourceFunctors...> 		mStpl;



	constexpr ExplicitUpdate(double dt, ConservedQuantityFunctor cf, NumericalFluxInFunctor fif, std::tuple<SourceFunctors...> stpl)
	: mdt(dt)
	, mCQF(cf)
	, mNFF(fif)
	, mStpl(stpl){};



	template <typename CellT>
	double operator()(const CellT & cl){

		// sources
		double src=0.0;
		Detail::for_each(mStpl, [&cl, &src](const auto & s){src += s.operator()(cl);});

		return mCQF.operator()(cl) + mdt*mNFF.operator()(cl) + mdt*src;
	}


	// alter the timestep
	void set_timestep(double dt){mdt = dt;};
};






// convenience function to make an explicit updater
template <class ConservedQuantityFunctor,
		  class NumericalFluxInFunctor,
		  class... SourceFunctors>
auto make_explicit_update(double dt, ConservedQuantityFunctor cf, NumericalFluxInFunctor fif, SourceFunctors... sf){
	return ExplicitUpdate<ConservedQuantityFunctor,
						  NumericalFluxInFunctor,
						  SourceFunctors...>(dt, cf, fif, std::make_tuple(sf...));
}





}// end namespace fivo

#endif