#ifndef BASE_MATH_HPP
#define BASE_MATH_HPP

#include <stddef.h>
#include <type_traits>
#include <vector>
#include <iostream>
#include <utility>
#include <concepts>
#include <cmath>

#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/iterator/iterator_adapter.hpp>
#include <boost/mpl/int.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/static_assert.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/type_index.hpp>
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
#include <boost/hana/is_empty.hpp>
#include <boost/hana/fuse.hpp>
#include <boost/hana/at.hpp>
#include <boost/hana/traits.hpp>
#include <boost/core/typeinfo.hpp>
#include <boost/hana/count.hpp>
#include <boost/hana/size.hpp>

namespace fusion = boost::fusion;
namespace mpl = boost::mpl;

template<int N, typename T>
struct Simplex;
template<size_t N, typename T>
class Vector_space;

template<int N>
class quaternon_type;

//template<int N>
//using quaternon_type = std::array<int, N>;

template <int N, int I, class Closure> requires (I == N)
constexpr void is_meta_loop(Closure& closure) {}

template <int N, int I, class Closure> requires (I < N)
constexpr void is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <int N, class Closure>
constexpr void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}

template<int N>
struct quaternon_prod_scalar {
    quaternon_prod_scalar(const quaternon_type<N>& arr_, quaternon_type<N>& res_, int val) :
        arr(arr_), res(res_), val(val) {}
    template <int I>
    void apply() {
        res[I] = arr[I]*val;
    }
private:
    const quaternon_type<N>& arr;
    quaternon_type<N>& res;
    int val;
};
template<int N>
struct quaternon_prod_self {
    quaternon_prod_self(quaternon_type<N>& arr_, int val) :
        arr(arr_), val(val) {}
    template <int I>
    void apply() {
        arr[I] *= val;
    }
private:
    quaternon_type<N>& arr;
    int val;
};
template<int N>
struct quaternon_sum {
    quaternon_sum(const quaternon_type<N>& arr_, quaternon_type<N>& res_) :
        arr(arr_), res(res_) {}
    template <int I>
    void apply() {
        res[I] += arr[I];
    }
private:
    const quaternon_type<N>& arr;
    quaternon_type<N>& res;
};
template<int N>
struct quaternon_pow {
    quaternon_pow(const quaternon_type<N>& arr_, typename quaternon_type<N>::value_type& tmp_) :
        arr(arr_), tmp(tmp_) {}
    template <int I>
    void apply() {
        tmp += arr[I]*arr[I];
    }
private:
    const quaternon_type<N>& arr;
    typename quaternon_type<N>::value_type& tmp;
};

template<int N>
class quaternon_type : public std::array<int, N> {
public:
    template<typename ...Args, std::enable_if_t<(sizeof... (Args) == N)>* = nullptr>
    quaternon_type(Args&&... args) : std::array<int, N>{std::forward<Args>(args)...} {}
    quaternon_type() : std::array<int, N>{} {}
    friend inline quaternon_type<N> operator*(const quaternon_type<N>& arr, int val) {
        quaternon_type<N> tmp_qv;
        quaternon_prod_scalar<N> closure(arr, tmp_qv, val);
        meta_loop<N>(closure);
        return tmp_qv;
    }

    friend inline quaternon_type<N>& operator*=(quaternon_type<N>& arr, int val) {
        quaternon_prod_self<N> closure(arr, val);
        meta_loop<N>(closure);
        return arr;
    }

    friend inline quaternon_type<N>& operator+=(quaternon_type<N>& arr, const quaternon_type<N>& other) {
        quaternon_sum<N> closure(other, arr);
        meta_loop<N>(closure);
        return arr;
    }

    friend inline std::ostream& operator << (std::ostream& o, const quaternon_type<N>& s) {
        for(const auto& x : s) {
            if(x != 0) o << x;
        }
        return o;
    }
};



template<typename Simplex, class = bool>
struct IsSimplex : mpl::false_ {};
template<typename Simplex>
struct IsSimplex<Simplex, typename Simplex::value > : mpl::true_ {};

template<class T, typename... Args>
inline constexpr bool is_same_v_ = (std::is_same_v<T, Args> && ...);

template< class T >
struct remove_cvref {
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};
template< class T >
using remove_cvref_t = typename remove_cvref<T>::type;


struct print_ {
    std::ostream& o_;
    print_(std::ostream& o) : o_(o) {}
    template<typename T>
    void operator()(T& t) const {
        o_ << t << " ";
    }
};

template<class Simpl, size_t... S> requires Simpl::space_type::is_vs
inline Simpl unpack_vector(const std::vector<typename Simpl::space_type>& vec, std::index_sequence<S...>) {
    Simpl simpl(vec[S]...);
    return simpl;
}

template <int I, int N1, typename IN, typename T>
struct for_one {
private:
    const IN& simplex_vertex;
public:
    std::vector<T> RESULT;
    for_one(const IN& simplex_) : simplex_vertex(simplex_)  {}
    template <int J>
    void apply() {
        if constexpr(I != J) {
            RESULT.push_back(fusion::at_c<J>(simplex_vertex));
        }
    }
};
template <int N1, typename OUT, typename IN, class Simpl>
struct for_two {
private:
    OUT& RESULT;
    const IN& simplex_vertex;
public:
    for_two(OUT& RESULT_, const IN& simplex_) : RESULT(RESULT_), simplex_vertex(simplex_)  {}
    template <int I>
    void apply() {
        for_one<I, N1, IN, typename Simpl::space_type> closure(simplex_vertex);
        meta_loop<N1>(closure);
        fusion::at_c<I>(RESULT) = unpack_vector<Simpl>(closure.RESULT, std::make_index_sequence<N1 - 1>());
    }
};
template <int N1, int N2, typename OUT, typename IN, class Simpl>
inline void for_all(OUT& RESULT_, const IN& simplex_) {
    for_two<N1, OUT, IN, Simpl> closure(RESULT_, simplex_);
    meta_loop<N2>(closure);
}

template <int I, int M, typename IN>
struct plus_one {
private:
    IN& simplex_vertex1;
    const IN& simplex_vertex2;
public:
    plus_one(IN& simplex1_, const IN& simplex2_) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_)  {}
    template <int J>
    void apply() {
        if constexpr(J == I) {
            fusion::at_c<J>(simplex_vertex1) += fusion::at_c<I>(simplex_vertex2);
        }
    }
};
template <int M, typename IN>
struct plus_two {
private:
    IN& simplex_vertex1;
    const IN& simplex_vertex2;
public:
    plus_two(IN& simplex1_, const IN& simplex2_) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_)  {}
    template <int I>
    void apply() {
        plus_one<I, M, IN> closure(simplex_vertex1, simplex_vertex2);
        meta_loop<M>(closure);
    }
};

template <int M, typename IN>
inline void plus_all(IN& simplex1_, const IN& simplex2_) {
    plus_two<M, IN> closure(simplex1_, simplex2_);
    meta_loop<M>(closure);
}
template <typename IN>
struct predicate_one {
private:
    const IN& simplex_vertex1;
    const IN& simplex_vertex2;
    int& index;
public:
    predicate_one(const IN& simplex1_, const IN& simplex2_, int& index) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), index(index)  {}
    template <int J>
    void apply() {
        if(fusion::at_c<J>(simplex_vertex1) == fusion::at_c<J>(simplex_vertex2)) {
            ++index;
//            std::cout << fusion::at_c<J>(simplex_vertex1) << " == " << fusion::at_c<J>(simplex_vertex2) << std::endl;
        }
//        else
//            std::cout << fusion::at_c<J>(simplex_vertex1) << " != " << fusion::at_c<J>(simplex_vertex2) << std::endl;
    }
};

template <int N, typename IN>
inline void predicate_all(const IN& simplex1_, const IN& simplex2_, int& index) {
    predicate_one<IN> closure(simplex1_, simplex2_, index);
    meta_loop<N>(closure);
}
template <int I, typename IN, typename IN_1>
struct predicateN_one {
private:
    const IN& simplex_vertex1;
    const IN_1& simplex_vertex2;
    int& index;
public:
    predicateN_one(const IN& simplex1_, const IN_1& simplex2_, int& index) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), index(index)  {}
    template <int J>
    void apply() {
        if(fusion::at_c<I>(simplex_vertex1) != fusion::at_c<J>(simplex_vertex2)) {
            ++index;
        }
    }
};
template <int M_2, typename IN, typename IN_1>
struct predicateN_two {
private:
    const IN& simplex_vertex1;
    const IN_1& simplex_vertex2;
    int& index;
public:
    predicateN_two(const IN& simplex1_, const IN_1& simplex2_, int& index) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), index(index)  {}
    template <int I>
    void apply() {
        predicateN_one<I, IN, IN_1> closure(simplex_vertex1, simplex_vertex2, index);
        meta_loop<M_2>(closure);
    }
};

template <int M_1, int M_2, typename IN, typename IN_1>
inline void predicateN_all(const IN& simplex1_, const IN_1& simplex2_, int& index) {
    predicateN_two<M_2, IN, IN_1> closure(simplex1_, simplex2_, index);
    meta_loop<M_1>(closure);
}
template <int I, typename SimN, typename SimN_1, typename Vector>
struct emplace_vector_one {
private:
    using value_type = typename Vector::value_type;
    const SimN& simplex_vertex1;
    const SimN_1& simplex_vertex2;
    Vector& vec;
public:
    emplace_vector_one(const SimN& simplex1_, const SimN_1& simplex2_, Vector& vec_) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), vec(vec_)  {}
    template <int J>
    void apply() {
        if(fusion::at_c<I>(simplex_vertex1) != fusion::at_c<J>(simplex_vertex2)) {
//            auto tmp1 = fusion::at_c<I>(simplex_vertex1).get_vector();
            auto tmp2 = fusion::at_c<J>(simplex_vertex2).get_vector();
            value_type tmp;
            fusion::for_each(tmp2, [this, &tmp](const auto& z) {
//                std::cout << boost::typeindex::type_id_with_cvr<decltype(z)>().pretty_name() << "\n";
                tmp.insert(z);
            });
            vec.emplace_back(tmp);
        }
    }
};
template <int M_2, typename SimN, typename SimN_1, typename Vector>
struct emplace_vector_two {
private:
    const SimN& simplex_vertex1;
    const SimN_1& simplex_vertex2;
    Vector& vec;
public:
    emplace_vector_two(const SimN& simplex1_, const SimN_1& simplex2_, Vector& vec_) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), vec(vec_)  {}
    template <int I>
    void apply() {
        emplace_vector_one<I, SimN, SimN_1, Vector> closure(simplex_vertex1, simplex_vertex2, vec);
        meta_loop<M_2>(closure);
    }
};

template <int M_1, int M_2, typename SimN, typename SimN_1, typename Vector> requires (decltype(fusion::size(std::declval<SimN>()))::value == M_1 &&
                                                                                       decltype(fusion::size(std::declval<SimN_1>()))::value == M_2)
inline void emplace_vector_all(const SimN& simplex1_, const SimN_1& simplex2_, Vector& vec) {
    emplace_vector_two<M_2, SimN, SimN_1, Vector> closure(simplex1_, simplex2_, vec);
    meta_loop<M_1>(closure);
}

namespace my {
    template<class T, int N>
    struct helper {
        static constexpr T pow_(const T x){
            return helper<T, N-1>::pow_(x) * x;
        }
    };

    template<class T>
    struct helper<T, 1> {
        static constexpr T pow_(const T x){
            return x;
        }
    };

    template<class T>
    struct helper<T, 0> {
        static constexpr T pow_(const T x){
            return T(1);
        }
    };
    template<int N, class T>
    T consteval  pow_(T const x) {
        return helper<T, N>::pow_(x);
    }
    template<int N>
    inline double norm(const quaternon_type<N>& arr) {
        typename quaternon_type<N>::value_type tmp{};
        quaternon_pow<N> closure(arr, tmp);
        meta_loop<N>(closure);
        return std::sqrt(tmp);
    }
}



#endif // BASE_MATH_HPP
