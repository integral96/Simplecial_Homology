#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include "base_math.hpp"

template <typename T>
struct boundary_predicate {
    T& o_;
    int& value;
    boundary_predicate(T& o, int& value) : o_(o), value(value) {}
    template<typename F>
    void operator()(F& t) const {
        value = t == o_ ? 0 : 1;
    }
};

template <typename IN>
struct boundary_one {
private:
    const IN& simplex_vertex;
public:
    boundary_one(const IN& simplex_) : simplex_vertex(simplex_)  {}
    template <int J>
    void apply() {
        auto tmp = fusion::at_c<J>(simplex_vertex);
        int value;
        if constexpr(my::pow_<J, int>(-1) == -1)
            std::cout << " - " << tmp;
        else
            std::cout << " + " << tmp;
//        fusion::for_each(simplex_vertex, boundary_predicate(tmp, value));
//        std::cout << "VALUE = " << value;
    }
};

template <int N, typename IN>
inline void boundary_all(const IN& simplex_) {
    boundary_one<IN> closure(simplex_);
    meta_loop<N>(closure);
    std::cout << std::endl;
}

template <int N, typename Simpl>
struct Chain
{
    Simpl simplex1;
    Simpl simplex2;
    char operat;
    Chain() {}
    Chain(const Simpl& simplex1) : simplex1(simplex1) {
    }
    Chain(const Simpl& simplex1, char operat, const Simpl& simplex2) : simplex1(simplex1), operat(operat), simplex2(simplex2) {
    }
    friend std::ostream& operator << (std::ostream& o, const Chain& s) {
        Simpl dff;
        o << "[ ";
        if(s.simplex2 == dff) o << s.simplex1;
        else o << s.simplex1 << " " << s.operat << " " << s.simplex2;
        o << " ]";
        return o;
    }
};

template <int N, typename Simpl>
struct boundary_chain
{
    const Simpl& chain;
    boundary_chain(const Simpl& chain) : chain(chain) {
        std::cout << "BOUNDARY = " << N << std::endl;
        boundary_all<N + 1, Simpl>(chain);
    }
    /* -------------------------------------------------------------------- */
    template <typename Simpl_sub, typename = std::enable_if_t<boost::mpl::size<Simpl_sub>::value == N>>
    int operator * (const boundary_chain<N - 1, Simpl_sub>& value) {
        return 1;
    }
};

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

template<int N, typename T>
struct Simplex {
    using value_type = T;
    using simpl_type = Simplex<N - 1, T>;
    static constexpr int dim = N;
private:
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, N + 1>;
    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, N + 1>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex() {};
    template<typename ...Args, typename = std::enable_if_t<(sizeof... (Args) == N + 1) && is_same_v_<T, remove_cvref_t<Args>...>>>
    Simplex(Args&&... args) : simplex_vertex(fusion::make_vector(std::forward<Args>(args)...)) {
        for_all<N + 1, N + 1, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
        std::cout << "SIMPLEX DIM = " << N << std::endl;
        fusion::for_each(simplex_vertex, print_(std::cout));
        std::cout << std::endl;
        std::cout << "SUB SIMPLEX DIM = " << N - 1 << std::endl;
        fusion::for_each(sub_simplex, print_(std::cout));
        std::cout << std::endl;

    }
    gen_vector_t const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<N, gen_vector_simpl_t> boundary() {
        return boundary_chain<N, gen_vector_simpl_t>(sub_simplex);
    }
    //operators==================================================================
    Chain<N, Simplex> operator + (const Simplex& other) {
        if(*this == other) {
            return Chain<N, Simplex>(*this);
        } else {
            return Chain<N, Simplex>(*this, '+', other);
        }
    }
    //operators==================================================================
    Chain<N, Simplex> operator - (const Simplex& other) {
        if(*this == other) {
            return Chain<N, Simplex>();
        } else {
            return Chain<N, Simplex>(*this, '-', other);
        }
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<N, gen_vector_t>(simplex_vertex, other.simplex_vertex, predicate_index);
        std::cout << "Operator = " << predicate_index << "; size = " << N + 1 << std::endl;
        if(predicate_index == N)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend std::ostream& operator << (std::ostream& o, const Simplex& s) {
        o << "{ ";
        fusion::for_each(s.simplex_vertex, print_(o));
        o << "}";
        return o;
    }
};
template<typename T>
struct Simplex<2, T> {
    using value_type = T;
    using simpl_type = Simplex<1, T>;
    static constexpr int dim = 2;
private:
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, 3>;
    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, 3>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex(){};
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>>
    Simplex(F&& arg1, F&& arg2, F&& arg3) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1),
                                                                            std::forward<F>(arg2),
                                                                            std::forward<F>(arg3))) {
        sub_simplex = fusion::make_vector(simpl_type(arg1, arg2), simpl_type(arg1, arg3), simpl_type(arg2, arg3));
//        std::cout << "DIM = " << 2 << std::endl;
//        fusion::for_each(simplex_vertex, print_(std::cout));
//        std::cout << std::endl;
    }
    gen_vector_t const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<2, gen_vector_simpl_t>  boundary() {
        return boundary_chain<2, gen_vector_simpl_t>(sub_simplex);
    }
    //operators==================================================================
    Chain<2, Simplex> operator + (const Simplex& other) {
        if(*this == other) {
            return Chain<2, Simplex>(*this);
        } else {
            return Chain<2, Simplex>(*this, '+', other);
        }
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<2, gen_vector_t>(simplex_vertex, other.simplex_vertex, predicate_index);
        if(predicate_index == 2)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend std::ostream& operator << (std::ostream& o, const Simplex& s) {
        o << "{ ";
        fusion::for_each(s.simplex_vertex, print_(o));
        o << "}";
        return o;
    }
};
template<typename T>
struct Simplex<1, T> {
    using value_type = T;
    using simpl_type = Simplex<0, T>;
    static constexpr int dim = 1;
private:
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, 2>;
    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, 2>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex(){};
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>>
    Simplex(F&& arg1, F&& arg2) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1),
                                                                  std::forward<F>(arg2))) {
        sub_simplex = fusion::make_vector(simpl_type(arg1), simpl_type(arg2));
//        std::cout << "DIM = " << 1 << std::endl;
//        fusion::for_each(simplex_vertex, print_(std::cout));
//        std::cout << std::endl;
    }
    gen_vector_t const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<1, gen_vector_simpl_t>  boundary() {
        return boundary_chain<1, gen_vector_simpl_t>(sub_simplex);
    }
    //operators==================================================================
    Chain<1, Simplex> operator + (const Simplex& other) {
        if(*this == other) {
            return Chain<1, Simplex>(*this);
        } else {
            return Chain<1, Simplex>(*this, '+', other);
        }
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<1, gen_vector_t>(simplex_vertex, other.simplex_vertex, predicate_index);
        if(predicate_index == 1)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend std::ostream& operator << (std::ostream& o, const Simplex& s) {
        o << "{ ";
        fusion::for_each(s.simplex_vertex, print_(o));
        o << "}";
        return o;
    }
};
template<typename T>
struct Simplex<0, T> {
    using value_type = T;
    using simpl_type = Simplex<-1, T>;
    static constexpr int dim = 0;
private:
    fusion::vector<T> simplex_vertex;
    fusion::vector<simpl_type> sub_simplex;
public:
    Simplex(){};
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>>
    Simplex(F&& arg1) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1))) {
        sub_simplex = fusion::make_vector(simpl_type(arg1));
//        std::cout << "DIM = " << 0 << std::endl;
//        fusion::for_each(simplex_vertex, print_(std::cout));
//        std::cout << std::endl;
    }
    fusion::vector<T> const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<0, fusion::vector<simpl_type>>  boundary() {
        return boundary_chain<0, fusion::vector<simpl_type>>(sub_simplex);
    }
    //operators==================================================================
    Chain<0, Simplex> operator + (const Simplex& other) {
        if(*this == other) {
            return Chain<0, Simplex>(*this);
        } else {
            return Chain<0, Simplex>(*this, '+', other);
        }
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<0, fusion::vector<T>>(simplex_vertex, other.simplex_vertex, predicate_index);
        if(predicate_index == 0)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend std::ostream& operator << (std::ostream& o, const Simplex& s) {
        o << "{ ";
        fusion::for_each(s.simplex_vertex, print_(o));
        o << "}";
        return o;
    }
};
template<typename T>
struct Simplex<-1, T> {
    using value_type = T;
    static constexpr int dim = -1;
    Simplex(){};
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, F>>>
    Simplex(F arg1 = F())  {
//        std::cout << "AUGUMENTACIA" << std::endl;
    }
};



#endif // SIMPLEX_HPP
