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
#include <boost/hana/count.hpp>
#include <variant>
#include <any>
#include <cmath>

#include "Simplex.hpp"
#include "Vector_space.hpp"


namespace hana = boost::hana;
namespace mpl  = boost::mpl;

struct print_type
{
private:
    std::ostream& out;
public:
    print_type(std::ostream& out_) : out(out_) {}
    template <class T>
    void operator() (T) const
    {
        auto const& ti = BOOST_CORE_TYPEID(T);
        out << boost::core::demangled_name(ti) << std::endl;
    }
};

template<class Vect, size_t... S>
inline auto unpack_to_tuple(const Vect& vec, std::index_sequence<S...>) {
//    Simpl simpl(vec[S]...);
    return hana::make_tuple(vec[S]...);
}


template<typename ...Args>
inline auto make_complex(Args&&... args) {
    constexpr auto complex_type = hana::make<hana::tuple_tag>(hana::type<Args>{}...);
    constexpr auto complex_dim  = hana::make<hana::tuple_tag>(hana::int_c<std::decay_t<Args>::dim>...);
    auto complex = hana::make<hana::tuple_tag>(std::forward<Args>(args)...);

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

    constexpr value_type operator()() const noexcept {
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
     constexpr auto get_complex_type() const {
         return hana::at(this->complex, hana::size_c<0>);
     }
     constexpr auto get_complex_dim() const {
         return hana::at(this->complex, hana::size_c<1>);
     }
     constexpr auto get_complex_() const {
         return hana::at(this->complex, hana::size_c<2>);
     }
public:
    Simplicial_complex(Complex& complex_) : complex(complex_) {}
    template<int N>
    auto get_complex() const {
        constexpr auto is_Ng = hana::compose(hana::trait<is_NN<N>::template is_N>, hana::decltype_);
        const auto ret = hana::filter(hana::at(complex, hana::size_c<2>), is_Ng);
//        hana::for_each(ret, [&](auto x) {
//            std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
//        });
        std::cout << std::endl;
        return ret;
    }
    friend std::ostream& operator << (std::ostream& out, const Simplicial_complex& cmpl) {
        hana::for_each(cmpl.get_complex_type(), print_type(out));
        out << std::endl;
        hana::for_each(cmpl.get_complex_(), [&](auto x) {
            out << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
        });
        out << std::endl;
        return out;
    }
};
template<int N, class Complex_N, class Vec_spc>
class boundary {
    const Complex_N& complex;
    using subsimplex_type = Simplex<N - 1, Vec_spc >;
    std::vector<subsimplex_type> vector;
public:
    boundary(const Complex_N& complex_) : complex(complex_) {
        hana::for_each(complex, [&](auto& x) {
            fusion::for_each(x.get_simplex(), [](auto y) {
                std::cout << y << ' ';
            });
            std::cout <<  '\n';
            auto tnp = x.boundary_sub();
            fusion::for_each(tnp, [](auto y) {
                std::cout << y << ' ';
            });
            std::cout <<  '\n';
            int p = 0;
            auto dth = fusion::at<mpl::int_<0>>(tnp);
            const auto dth_tmp = fusion::at<mpl::int_<0>>(tnp);
            fusion::for_each(tnp, [&dth, &dth_tmp, j = p](auto& z) mutable {
                if(dth_tmp != z) {
                    std::cout << dth_tmp << "; " << z << "; " << j << '\n';
                    dth = dth + z*std::pow(-1, j);
                }
                ++j;
            });
            vector.emplace_back(dth);
            std::cout << dth << '\n';
            std::cout <<  '\n';
        });
    }
    auto get() const {
        return unpack_to_tuple(vector, std::make_index_sequence<2>());
    }
};

#endif // SIMPLICIAL_COMPLEX_HPP
