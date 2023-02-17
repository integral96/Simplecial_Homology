#ifndef VECTOR_SPACE_HPP
#define VECTOR_SPACE_HPP

#include <boost/shared_array.hpp>
#include <boost/smart_ptr/detail/sp_nullptr_t.hpp>
#include "base_math.hpp"

template<typename T>
struct product_scalr_ {
    T& o_;
    product_scalr_(T& o) : o_(o) {}
    template<typename Vec>
    void operator()(Vec& t) const {
        t *= o_;
    }
};


template <typename IN>
struct predicate_ {
private:
    const IN& simplex_vertex1;
    const IN& simplex_vertex2;
    int& index;
public:
    predicate_(const IN& simplex1_, const IN& simplex2_, int& index) :
        simplex_vertex1(simplex1_), simplex_vertex2(simplex2_), index(index)  {}
    template <int J>
    void apply() {
        if(fusion::at_c<J>(simplex_vertex1) == fusion::at_c<J>(simplex_vertex2)) {
            ++index;
        }
    }
};

template <int N, typename IN>
inline void predicate_all_(const IN& simplex1_, const IN& simplex2_, int& index) {
    predicate_<IN> closure(simplex1_, simplex2_, index);
    meta_loop<N>(closure);
}

template<size_t N, typename T>
class Vector_space : public boost::shared_array<T> {
public:
    using value_type = T;
    static constexpr int dim = N;
private:
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, N>;
    gen_vector_t vector_;
public:
    Vector_space(){}
    template<typename... Args, typename = std::enable_if_t<(sizeof... (Args) == N) && is_same_v_<T, remove_cvref_t<Args>...>>>
    explicit Vector_space(Args&& ... args)
        : vector_(fusion::make_vector(std::forward<Args>(args)...)) {}
    /* -------------------------------------------------------------------- */
    template< typename F, typename = std::enable_if_t<std::is_same_v<T, F>>>
    Vector_space& operator *= (F value) {
        fusion::for_each(vector_, product_scalr_(value));
        return *this;
    }

    //operators==================================================================
    Vector_space& operator += (const Vector_space& other) {
        T value = 3;
        int I = 0;
        fusion::for_each(vector_, [&other, P_I = I](auto& x) mutable {
            int J = 0;
            fusion::for_each(other.vector_, [&x, &P_I, P_J = J](auto y) mutable {
                if(P_I == P_J) x += y;
                P_J++;
            });
            P_I++;
        });
        return *this;
    }
    template< typename F, typename = std::enable_if_t<std::is_same_v<T, F>>>
    friend Vector_space operator*(Vector_space lhs, F rhs)
    {
        lhs *= rhs;
        return lhs;
    }
    bool operator == (const Vector_space& other) const {
        int predicate_index{};
        predicate_all_<N, gen_vector_t>(vector_, other.vector_, predicate_index);
        if(predicate_index == N)
            return true;
        else return false;
    }
    friend std::ostream& operator << (std::ostream& o, const Vector_space& s) {
        o << "{ ";
        fusion::for_each(s.vector_, print_(o));
        o << "}";
        return o;
    }
};

#endif // VECTOR_SPACE_HPP
