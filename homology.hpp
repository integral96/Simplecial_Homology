#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"
#include "group_comm.hpp"

template <int N, int I, class Closure> requires (I == N)
constexpr void is_chain_loop(Closure& closure) {}

template <int N, int I, class Closure> requires (I > N)
constexpr void is_chain_loop(Closure& closure) {
    closure.template apply<I>();
    is_chain_loop<N, I - 1>(closure);
}
template <int N, class Closure>
constexpr void chain_loop(Closure& closure) {
    is_chain_loop<-1, N>(closure);
}
template <int N, class Closure>
constexpr void bound_loop(Closure& closure) {
    is_chain_loop<0, N>(closure);
}

template <typename Simplcl_cmplx>
struct _each {
private:
    const Simplcl_cmplx& simplcl_cmplx;
public:
    _each(const Simplcl_cmplx& simplcl_cmplx_) :
        simplcl_cmplx(simplcl_cmplx_)  {}
    template <int J>
    void apply() {
        std::cout << "C_" << J << " =\n";
        hana::for_each(simplcl_cmplx.template chain_unit<J>(), [](auto x) {
            std::cout << "Simplex dim = " << decltype(x)::dim << "; " << x << '\n';
        });
        std::cout << "=======>" << '\n';
    }
};
template <typename Simplcl_cmplx>
struct assert_boundary_ {
private:
    const Simplcl_cmplx& simplcl_cmplx;
public:
    assert_boundary_(const Simplcl_cmplx& simplcl_cmplx_) :
        simplcl_cmplx(simplcl_cmplx_)  {}
    template <int J>
    void apply() {
        if constexpr(J - 1 == 0) {
            std::cout << "boundary<" << J << ">(boundary<" << J - 1 << ">) \033[1;32mAUGUMENTACIA\033[0m" << std::endl;
        } else {
            auto chain = simplcl_cmplx.template chain_unit<J>();
            auto chain_N  = boundary<J>(chain); //Kernel
            auto chain_N_1  = boundary<J - 1>(chain_N);
            std::cout << "boundary<" << J << ">(boundary<" << J - 1 << ">) " << '\n';
            hana::for_each(chain_N_1, [&](auto x) {
                BOOST_ASSERT_MSG(x.empty(), "Двойное взятие граничного оператора должно быть равно нулую.");

            });
            std::cout << "\033[1;32mДвойное взятие граничного оператора должно быть равно нулую; Выполнено!\033[0m" << std::endl;
        }
    }
};
template <typename Simplcl_cmplx>
struct homology_group_ {
private:
    const Simplcl_cmplx& simplcl_cmplx;
public:
    homology_group_(const Simplcl_cmplx& simplcl_cmplx_) :
        simplcl_cmplx(simplcl_cmplx_)  {}
    template <int J>
    void apply() {
        if constexpr(J == -1) {
            std::cout << J << "\033[1;32m AUGUMENTACIA\033[0m" << std::endl;
        } else {
            auto chain_Im = simplcl_cmplx.template chain_unit<J + 1>();
            auto Img_N1  = Image<J + 1>(chain_Im);
            auto chain_Ker = simplcl_cmplx.template chain_unit<J>();
            auto Ker_N1  = Kernel<J>(chain_Ker);
            auto homology = quotient(Ker_N1, Img_N1);
            std::cout << "\033[1;32mГруппа гомологий: K = \033[0m" << J << std::endl;
            for(const  auto& x : homology) {
                std::cout << "<";
                for(const auto& y : x) {
                    std::cout << y << " ";
                }
                std::cout << ">\n";
            }
        }
    }
};

template <int N, typename Simplcl_cmplx>
inline void chain_complex(const Simplcl_cmplx& simplcl_cmplx) {
    _each<Simplcl_cmplx> closure(simplcl_cmplx);
    chain_loop<N>(closure);
}
template <int N, typename Simplcl_cmplx>
inline void assert_boundary(const Simplcl_cmplx& simplcl_cmplx) {
    assert_boundary_<Simplcl_cmplx> closure(simplcl_cmplx);
    bound_loop<N>(closure);
}
template <int N, typename Simplcl_cmplx>
inline void homology_group(const Simplcl_cmplx& simplcl_cmplx) {
    homology_group_<Simplcl_cmplx> closure(simplcl_cmplx);
    chain_loop<N>(closure);
}

#endif // HOMOLOGY_HPP
