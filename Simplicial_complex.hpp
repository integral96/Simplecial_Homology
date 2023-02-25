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
#include <boost/hana/functional/fix.hpp>
#include <boost/hana/fuse.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/traits.hpp>
#include <boost/core/typeinfo.hpp>
#include <boost/hana/count.hpp>
#include <boost/hana/size.hpp>
#include <variant>
#include <any>
#include <cmath>
#include <set>

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
    return hana::make_tuple(vec[S]...);
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

template<typename ...Args>
inline constexpr auto make_complex(Args&&... args) {
//    typedef boost::mpl::vector<Args...> MyMplVector;
    constexpr auto complex_type = hana::make<hana::tuple_tag>(hana::type<Args>{}...);;
    constexpr auto complex_dim  = hana::make<hana::tuple_tag>(hana::int_c<std::decay_t<Args>::dim>...);
    return hana::make<hana::tuple_tag>(std::forward<Args>(args)...);
}

template <class Complex>
class Simplicial_complex {
    using complex_type = Complex;
    template<int N>
    struct is_NN {
        template <class _Ty>
        struct is_N : bool_constant<is_N_v<_Ty, N>> {};
    };
    complex_type& complex;
public:

    template<int N>
    static constexpr auto get_count() {
        constexpr auto is_Ng = hana::compose(hana::trait<is_NN<N>::template is_N>, hana::decltype_);
        return hana::count_if(Complex{}, is_Ng);
    }

public:
    Simplicial_complex() = delete;
    Simplicial_complex(Complex& complex_) : complex{complex_} {}
    template<int N>
    constexpr auto chain_unit() const {
        constexpr auto is_Ng = hana::compose(hana::trait<is_NN<N>::template is_N>, hana::decltype_);
        return hana::filter(complex, is_Ng);
    }

    friend std::ostream& operator << (std::ostream& out, const Simplicial_complex& cmpl) {
        out << std::endl;
        hana::for_each(cmpl.complex, [&](auto x) {
            out << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
        });
        out << std::endl;
        return out;
    }
};
template<typename Left, typename Right>
inline  auto curry(Left& left, Right& right) {
    return hana::equal(left, right);
}

template<int N, class Complex>
inline constexpr auto boundary(const Complex& chain) {
    constexpr size_t tupleSize = decltype(hana::size(std::declval<Complex>()))::value;
    typedef typename std::remove_cvref_t<
                                        decltype(hana::front(std::declval<Complex>()))
                                        >::value_type Vector_type;
    using ring_type = typename Vector_type::value_type;
//    auto const& ti = BOOST_CORE_TYPEID(Vector_type);
//    std::cout << boost::core::demangled_name(ti) << std::endl;
//    std::cout << "Количество симплексов размерности: " << N << "; равно: " << tupleSize << '\n';
    if constexpr(N == 0) {
        using subsimplex_type = Simplex<0, Vector_type>;
        std::vector<subsimplex_type> vector;
        hana::for_each(chain, [&vector](const auto& x) {
            auto tnp = x.get_simplex();
            std::cout << "Разложение в " << N << "-цепь." << '\n';
            int k = 0;
            fusion::for_each(tnp, [j = k](const auto& y) mutable {
                auto znak = std::pow(-1, j);
                if(znak == -1)
                    std::cout << "-" << y << '\n';
                else
                    std::cout << "+" << y << '\n';
                ++j;
            });
            std::cout <<  '\n';
            int p = 0;
            auto dth = fusion::at<mpl::int_<0>>(tnp);
            fusion::for_each(tnp, [&dth, j = p](auto& z) mutable {
                auto z_tmp = z*ring_type(std::pow(-1, j));
                if(z_tmp != dth) {
                    dth += z_tmp;
                }
                ++j;
            });
            vector.emplace_back(dth);
        });
        return unpack_to_tuple(vector, std::make_index_sequence<tupleSize>());
    } else {
        using subsimplex_type = Simplex<N - 1, Vector_type>;
        std::vector<subsimplex_type> vector;
        hana::for_each(chain, [&vector](const auto& x) {
            auto tnp = x.boundary_sub();
            std::cout << "Разложение в " << N << "-цепь." << '\n';
            int k = 0;
            fusion::for_each(tnp, [j = k](const auto& y) mutable {
                auto znak = std::pow(-1, j);
                if(znak == -1)
                    std::cout << "-" << y << '\n';
                else
                    std::cout << "+" << y << '\n';
                ++j;
            });
            std::cout <<  '\n';
            int p = 0;
            auto dth = fusion::at<mpl::int_<0>>(tnp);
            fusion::for_each(tnp, [&dth, j = p](auto& z) mutable {
                auto z_tmp = z*std::pow(-1, j);
                if(z_tmp != dth) {
                    dth = dth + z_tmp;
                }
                ++j;
            });
            vector.emplace_back(dth);
        });
        return unpack_to_tuple(vector, std::make_index_sequence<tupleSize>());
    }
}
template<int N, class Complex>
inline constexpr auto Kernel(const Complex& chain) {
    return boundary<N>(chain);
}
template<int N, class Complex>
inline constexpr auto Image(const Complex& chain) {
    return boundary<N>(chain);
}
using namespace hana::literals;
template<class Ker, class Im>
inline constexpr auto quotient(const Ker& ker, const Im& im) {
    constexpr size_t tupleSize = decltype(hana::size(std::declval<Im>()))::value;
    typedef typename std::remove_cvref_t<
                                        decltype(hana::front(std::declval<Im>()))
                                        > Simplex_type;
    typedef typename std::remove_cvref_t<
                                        decltype(hana::front(std::declval<Im>()))
                                        >::value_type Vector_type;
    using subsimplex_type = Simplex<Simplex_type::dim, Vector_type>;
    using group_type = std::set<typename Vector_type::value_type>;
//    auto const& ti = BOOST_CORE_TYPEID(Simplex_type);
//    std::cout << boost::core::demangled_name(ti) << "; " << Simplex_type::dim << std::endl;
//    hana::for_each(ker, print_type(std::cout));
//    hana::for_each(im, print_type(std::cout));

    std::vector<group_type> vector_ker;
    hana::for_each(ker, [&vector_ker, &im](const auto& x) {
        hana::for_each(im, [&vector_ker, &x](const auto& y) {
            y.template emplace_vector<std::vector<group_type>>(x, vector_ker);
        });
    });
    std::sort( vector_ker.begin(), vector_ker.end() );
    vector_ker.erase( std::unique( vector_ker.begin(), vector_ker.end() ), vector_ker.end() );
    return vector_ker;
}

#endif // SIMPLICIAL_COMPLEX_HPP
