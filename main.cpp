#include <iostream>
#include <array>

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"
static constexpr int C_N = 5;
static constexpr int N = 4;
static constexpr int M = 0;

using vector_space = Vector_space<4, quaternon_type>;

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



int main()
{
    quaternon_type e0{1, 0, 0, 0};
    quaternon_type e1{0, 1, 0, 0};
    quaternon_type e2{0, 0, 1, 0};
    quaternon_type e3{0, 0, 0, 1};

    vector_space d0( e0*2,    e1*(-1), e2*(-5), e3);
    vector_space d1( e0*(-1), e1*3,    e2*(-2), e3*5);
    vector_space d2( e0*1,    e1*(-5), e2*3,    e3*2);
    vector_space d3( e0*(-1), e1*(-6), e2*2,    e3*4);
    vector_space d4( e0*(-2), e1,      e2*(-3), e3*2);

    Simplex<5, vector_space> simplex5D_0(d0, d1, d2, d3, d4, d2);
    Simplex<5, vector_space> simplex5D_1(d2, d4, d0, d3, d0, d1);
    Simplex<5, vector_space> simplex5D_2(d2, d3, d1, d3, d4, d0);

    Simplex<4, vector_space> simplex4D_0(d0, d1, d2, d3, d4);
    Simplex<4, vector_space> simplex4D_1(d2, d4, d0, d3, d0);

    Simplex<3, vector_space> simplex3D_0(d0, d1, d2, d3);
    Simplex<3, vector_space> simplex3D_1(d0, d1, d2, d3);
    Simplex<3, vector_space> simplex3D_2(d1, d3, d2, d4);

    Simplex<2, vector_space> simplex2D_0(d0, d1, d2);
    Simplex<2, vector_space> simplex2D_1(d0, d1, d2);

    Simplex<1, vector_space> simplex1D_0(d0, d1);
    Simplex<1, vector_space> simplex1D_1(d4, d3);

    Simplex<0, vector_space> simplex0D_0(d0);
    Simplex<0, vector_space> simplex0D_1(d1);
    Simplex<0, vector_space> simplex0D_2(d2);


    std::cout << "\033[1;31mСоздание комплекса до размерностей \033[0m" << N << "-х" << std::endl;
    auto complex = make_complex(simplex5D_0, simplex5D_1, simplex5D_2,
                                simplex4D_0, simplex4D_1, simplex3D_0,
                                simplex3D_1, simplex3D_2, simplex2D_0,
                                simplex2D_1, simplex1D_0, simplex1D_1,
                                simplex0D_0, simplex0D_1, simplex0D_2);

    using Simpl_complex_type = Simplicial_complex<decltype (complex)>;
    std::cout << "\033[1;31mСоздание цепного комплекса из аблевых групп:\033[0m" << std::endl;
    Simpl_complex_type complex_simpl(complex);
    chain_complex<C_N, Simpl_complex_type>(complex_simpl);
    std::cout << "\033[1;31mВыполнено\033[0m" << std::endl;
    std::cout << "\033[1;31mПроверка двойного взятие граничного оператора по цепному комплексу:\033[0m" << std::endl;
    assert_boundary<C_N, Simpl_complex_type>(complex_simpl);
    std::cout << "\033[1;31mВыполнено\033[0m" << std::endl;
    std::cout << "\033[1;31mРасчет гомологий:\033[0m" << std::endl;
    homology_group<N, Simpl_complex_type>(complex_simpl);

    std::cout << "\033[1;31mВыполнено\033[0m" << std::endl;

    return 0;
}
