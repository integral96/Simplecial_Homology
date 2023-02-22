#ifndef BASE_MATH_HPP
#define BASE_MATH_HPP

#include <stddef.h>
#include <type_traits>
#include <vector>
#include <iostream>
#include <utility>

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

namespace fusion = boost::fusion;
namespace mpl = boost::mpl;

template<class T, typename... Args>
inline constexpr bool is_same_v_ = (std::is_same_v<T, Args> && ...);

template< class T >
struct remove_cvref {
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};
template< class T >
using remove_cvref_t = typename remove_cvref<T>::type;

template <int N, int I, class Closure>
typename std::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

template <int N, int I, class Closure>
typename std::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <int N, class Closure>
void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}
struct print_ {
    std::ostream& o_;
    print_(std::ostream& o) : o_(o) {}
    template<typename T>
    void operator()(T& t) const {
        o_ << t << " ";
    }
};

template<class Simpl, size_t... S>
inline Simpl unpack_vector(const std::vector<typename Simpl::value_type>& vec, std::index_sequence<S...>) {
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
        for_one<I, N1, IN, typename Simpl::value_type> closure(simplex_vertex);
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
            fusion::at_c<J>(simplex_vertex1) += fusion::at_c<J>(simplex_vertex2);
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
//            std::cout << fusion::at_c<I>(simplex_vertex1) << " == " << fusion::at_c<J>(simplex_vertex2) << std::endl;
        }
//        else
//            std::cout << fusion::at_c<J>(simplex_vertex1) << " != " << fusion::at_c<J>(simplex_vertex2) << std::endl;
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
    T constexpr pow_(T const x) {
        return helper<T, N>::pow_(x);
    }
}



#endif // BASE_MATH_HPP
