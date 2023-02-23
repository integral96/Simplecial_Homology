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


template<size_t N, typename T>
class Vector_space {
public:
    using value_type = T;
    static constexpr int dim = N;
private:
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, N>;
    gen_vector_t vector_;
public:
    Vector_space(){}
    template<typename... Args, typename = std::enable_if_t<(sizeof... (Args) == N) && is_same_v_<T, remove_cvref_t<Args>...>>>
    Vector_space(Args&& ... args)
        : vector_(fusion::make_vector(std::forward<Args>(args)...)) {}
    /* -------------------------------------------------------------------- */
    template< typename F, typename = std::enable_if_t<std::is_same_v<T, F>>>
    Vector_space& operator *= (F value) {
        fusion::for_each(vector_, product_scalr_(value));
        return *this;
    }
    bool get_zero() const {
        return vector_ == gen_vector_t();
    }
    gen_vector_t const& get_vector() const {
        return vector_;
    }

    //operators==================================================================
    Vector_space& operator += (const Vector_space& other) {
        plus_all<N, gen_vector_t>(vector_, other.vector_);
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
        predicate_all<N, gen_vector_t>(vector_, other.vector_, predicate_index);
        if(predicate_index == N)
            return true;
        else return false;
    }
    bool operator != (const Vector_space& other) const {
        return !this->operator==( other );
    }
    friend std::ostream& operator << (std::ostream& o, const Vector_space& s) {
        o << "{ ";
        fusion::for_each(s.vector_, print_(o));
        o << "}";
        return o;
    }
};

#endif // VECTOR_SPACE_HPP
