#ifndef SIMPLICIAL_COMPLEX_HPP
#define SIMPLICIAL_COMPLEX_HPP

#include <boost/hana/assert.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/core/make.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/map.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/type.hpp>
#include <boost/hana/fold_left.hpp>
#include <boost/hana/insert.hpp>
#include <boost/hana/remove_if.hpp>
#include <boost/hana/functional/compose.hpp>
#include <boost/hana/functional/iterate.hpp>
#include <boost/hana/functional/placeholder.hpp>
#include <boost/hana/fuse.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/traits.hpp>
#include <boost/core/typeinfo.hpp>
#include <variant>

#include "Simplex.hpp"


namespace hana = boost::hana;
namespace mpl  = boost::mpl;

struct print_type
{
    template <class T>
    void operator() (T) const
    {
        auto const& ti = BOOST_CORE_TYPEID(T);
        std::cout << boost::core::demangled_name(ti) << std::endl;
    }
};

template<typename ...Args>
inline auto make_complex(Args&&... args) {
    constexpr auto complex_type = hana::make<hana::tuple_tag>(hana::type<Args>{}...);
    constexpr auto complex_dim  = hana::make<hana::tuple_tag>(hana::int_c<std::decay_t<Args>::dim>...);
    auto complex = hana::make<hana::tuple_tag>(std::forward<Args>(args)...);
    hana::for_each(complex_type, print_type());
    hana::for_each(complex_dim, [&](auto x) {
        std::cout << "Simplex dim = " << x << '\n';
    });
    std::cout << std::endl;
//    hana::for_each(complex, [&](auto x) {
//        if constexpr(decltype(x)::dim == 2)
//        std::cout << x << ' ';
//    });
    std::cout << std::endl;
    return hana::make_tuple(complex_type, complex_dim, complex);
}
template <class _Ty, _Ty _Val>
struct N_constant {
    static constexpr _Ty value = _Val;

    using value_type = _Ty;
    using type       = N_constant;

    constexpr operator value_type() const noexcept {
        return value;
    }

    _NODISCARD constexpr value_type operator()() const noexcept {
        return value;
    }
};

template <bool _Val>
using bool_constant = N_constant<bool, _Val>;

using true_type  = bool_constant<true>;
using false_type = bool_constant<false>;

template <class _Ty, int N>
inline constexpr bool is_N_v = _Ty::dim == N;

template <class Complex>
class Simplicial_complex {
private:
     Complex& complex;
     template<int N>
     struct is_NN {
         template <class _Ty>
         struct is_N : bool_constant<is_N_v<_Ty, N>> {};
     };

public:
    Simplicial_complex(Complex& complex_) : complex(complex_) {}
    template<int N>
    auto get_complex() const {

        constexpr auto is_N__ = hana::compose(hana::trait<is_NN<N>::template is_N>, hana::decltype_);
        auto ret = hana::filter(hana::at(this->complex, hana::size_c<2>), is_N__);
        hana::for_each(ret, [&](auto x) {
            std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
        });
        std::cout << std::endl;
    }
};

#endif // SIMPLICIAL_COMPLEX_HPP
