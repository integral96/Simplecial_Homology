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

template<int N, typename T>
struct Simplex;

template <int N, typename T>
struct Chain
{
    using simpl_type = Simplex<N - 1, T>;
    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, N + 1>;
private:
    gen_vector_simpl_t& sub_simpleces;
public:
    Chain() {}
    Chain(gen_vector_simpl_t& sub_simpleces_) : sub_simpleces(sub_simpleces_) {
    }
    friend std::ostream& operator << (std::ostream& o, const Chain& s) {
        o << "[ ";
        fusion::for_each(s.sub_simpleces, [&o](auto x) {
            o << x << ' ';
        });
        o << " ] ";
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


template<int N, typename T>
struct Simplex {
    using value_type = T;
    using ring_type  = typename T::value_type;
    using simpl_type = Simplex<N - 1, T>;
    static constexpr int dim = N;
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, N + 1>;
private:

    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, N + 1>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex() {}
    template<typename ...Args, std::enable_if_t<(sizeof... (Args) == N + 1) && is_same_v_<T, remove_cvref_t<Args>...>>* = nullptr>
    Simplex(Args&&... args) : simplex_vertex(fusion::make_vector(std::forward<Args>(args)...)) {
        for_all<N + 1, N + 1, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    Simplex(Simplex const& simpl) : simplex_vertex(simpl.simplex_vertex) {
        for_all<N + 1, N + 1, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    gen_vector_t const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<N, gen_vector_simpl_t> boundary() {
        return boundary_chain<N, gen_vector_simpl_t>(sub_simplex);
    }
    gen_vector_simpl_t  boundary_sub() const {
        return sub_simplex;
    }
    bool empty() const {
        int predicate_index{};
        fusion::for_each(simplex_vertex, [&predicate_index](const auto& x) {
            if(x == T()) ++predicate_index;
        });
        if(predicate_index == N + 1)
            return true;
        else return false;
    }
    //operators==================================================================

    Simplex& operator + (Simplex& other) {
        plus_all<N + 1, gen_vector_t>(simplex_vertex, other.simplex_vertex);
        return *this;
    }
    Simplex& operator * (ring_type other) {
        fusion::for_each(simplex_vertex, [&other](auto& x) {
            x *= other;
        });
        return *this;
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<N + 1, gen_vector_t>(simplex_vertex, other.simplex_vertex, predicate_index);
        if(predicate_index == N + 1)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend bool operator != (const Simplex& left, const Simplex<N - 1, T>& right) {
        int predicate_index{};
        Simplex<N, T> sd;
        predicateN_all<N + 1, N, typename Simplex::gen_vector_t, typename Simplex<N - 1, T>::gen_vector_t>(left.simplex_vertex, right.get_simplex(), predicate_index);
        if(predicate_index >= N*(N + 1))
            return true;
        else return false;
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
    using ring_type  = typename T::value_type;
    using simpl_type = Simplex<1, T>;
    static constexpr int dim = 2;
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, 3>;
private:

    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, 3>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex(){}
    template<typename F, std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>* = nullptr>
    Simplex(F&& arg1, F&& arg2, F&& arg3) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1),
                                                                            std::forward<F>(arg2),
                                                                            std::forward<F>(arg3))) {
        for_all<3, 3, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    Simplex(Simplex const& simpl) : simplex_vertex(simpl.simplex_vertex) {
        for_all<3, 3, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    gen_vector_t const& get_simplex() const {
//        std::cout << "BINGO DIM = " << 2 << std::endl;
        return simplex_vertex;
    }
    boundary_chain<2, gen_vector_simpl_t>  boundary() {
        return boundary_chain<2, gen_vector_simpl_t>(sub_simplex);
    }
    gen_vector_simpl_t  boundary_sub()  const {
        return sub_simplex;
    }
    bool empty() const {
        int predicate_index{};
        fusion::for_each(sub_simplex, [&predicate_index](const auto& x) {
            if(x.empty()) ++predicate_index;
        });
        if(predicate_index == 3)
            return true;
        else return false;
    }
    //operators==================================================================

    Simplex& operator + (const Simplex& other) {
        plus_all<3, gen_vector_t>(simplex_vertex, other.simplex_vertex);
        return *this;
    }
    Simplex& operator * (ring_type other) {
        fusion::for_each(simplex_vertex, [&other](auto& x) {
            x *= other;
        });
        return *this;
    }
    //comparison operators==================================================================
    bool operator == (const Simplex& other) const {
        int predicate_index{};
        predicate_all<3, gen_vector_t>(simplex_vertex, other.simplex_vertex, predicate_index);
        if(predicate_index == 3)
            return true;
        else return false;
    }
    bool operator != (const Simplex& other) const {
        return !this->operator==( other );
    }
    friend bool operator != (const Simplex& left, const Simplex<1, T>& right) {
        int predicate_index{};
        predicateN_all<3, 2, typename Simplex::gen_vector_t, typename Simplex<1, T>::gen_vector_t>(left.simplex_vertex, right.get_simplex(), predicate_index);
        if(predicate_index >= 2*3)
            return true;
        else return false;
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
    using ring_type  = typename T::value_type;
    using simpl_type = Simplex<0, T>;
    static constexpr int dim = 1;
    using gen_vector_t = boost::mp11::mp_repeat_c<fusion::vector<T>, 2>;
private:

    using gen_vector_simpl_t = boost::mp11::mp_repeat_c<fusion::vector<simpl_type>, 2>;
    gen_vector_t simplex_vertex;
    gen_vector_simpl_t sub_simplex;
public:
    Simplex(){}
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>>
    Simplex(F&& arg1, F&& arg2) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1),
                                                                  std::forward<F>(arg2))) {
        for_all<2, 2, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    Simplex(Simplex const& simpl) : simplex_vertex(simpl.simplex_vertex) {
        for_all<2, 2, gen_vector_simpl_t, gen_vector_t, simpl_type>(sub_simplex, simplex_vertex);
    }
    gen_vector_t const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<1, gen_vector_simpl_t>  boundary() {
        return boundary_chain<1, gen_vector_simpl_t>(sub_simplex);
    }
    gen_vector_simpl_t  boundary_sub() const {
        return sub_simplex;
    }
    bool empty() const {
        int predicate_index{};
        fusion::for_each(sub_simplex, [&predicate_index](const auto& x) {
            if(x.empty()) ++predicate_index;
        });
        if(predicate_index == 2)
            return true;
        else return false;
    }
    //operators==================================================================

    Simplex& operator + (Simplex& other) {
        plus_all<2, gen_vector_t>(simplex_vertex, other.simplex_vertex);
        return *this;
    }
    Simplex& operator * (ring_type other) {
        fusion::for_each(simplex_vertex, [&other](auto& x) {
            x *= other;
        });
        return *this;
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
    friend bool operator != (const Simplex& left, const Simplex<0, T>& right) {
        int predicate_index{};
        predicateN_all<2, 1, typename Simplex::gen_vector_t, typename Simplex<0, T>::gen_vector_t>(left.simplex_vertex, right.get_simplex(), predicate_index);
        if(predicate_index == 1)
            return true;
        else return false;
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
    using ring_type  = typename T::value_type;
    using simpl_type = Simplex<-1, T>;
    static constexpr int dim = 0;
    using gen_vector_t = fusion::vector<T>;
private:
    fusion::vector<T> simplex_vertex;
    fusion::vector<simpl_type> sub_simplex;
public:
    Simplex(){}
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, remove_cvref_t<F>>>>
    Simplex(F&& arg1) : simplex_vertex(fusion::make_vector(std::forward<F>(arg1))) {
        for_all<1, 1, fusion::vector<simpl_type>, fusion::vector<T>, simpl_type>(sub_simplex, simplex_vertex);
//        sub_simplex = fusion::make_vector(simpl_type(arg1));
    }
    Simplex(Simplex const& simpl) : simplex_vertex(simpl.simplex_vertex) {
        for_all<1, 1, fusion::vector<simpl_type>, fusion::vector<T>, simpl_type>(sub_simplex, simplex_vertex);
    }
    fusion::vector<T> const& get_simplex() const {
        return simplex_vertex;
    }
    boundary_chain<0, fusion::vector<simpl_type>>  boundary() {
        return boundary_chain<0, fusion::vector<simpl_type>>(sub_simplex);
    }
    fusion::vector<simpl_type>  boundary_sub() const {
        return sub_simplex;
    }
    bool empty() const {
        int predicate_index{};
        fusion::for_each(simplex_vertex, [&predicate_index](const auto& x) {
//            std::cout << typeid(x).name() << '\n';
            if(x.get_zero()) ++predicate_index;
        });
        if(predicate_index == 1)
            return true;
        else return false;
    }
    //operators==================================================================

    Simplex& operator + (Simplex& other) {
        plus_all<1, fusion::vector<T>>(simplex_vertex, other.simplex_vertex);
        return *this;
    }
    Simplex& operator * (ring_type other) {
        fusion::for_each(simplex_vertex, [&other](auto& x) {
            x *= other;
        });
        return *this;
    }
    //comparison operators==================================================================
    friend bool operator != (const Simplex<0, T>& left, const Simplex<-1, T>& right) {
        int predicate_index{};
        predicateN_all<1, 0, typename Simplex<0, T>::gen_vector_t, typename Simplex<-1, T>::gen_vector_t>(left.simplex_vertex, right.get_simplex(), predicate_index);
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
    using gen_vector_t = fusion::vector<T>;
    Simplex(){}
    template<typename F, typename = std::enable_if_t<std::is_same_v<T, F>>>
    Simplex(F arg1 = F())  {
//        std::cout << "AUGUMENTACIA" << std::endl;
    }
    gen_vector_t get_simplex() {
        return gen_vector_t();
    }
};

#endif // SIMPLEX_HPP
