#include <iostream>
#include <array>

#include "homology.hpp"

static constexpr int C_N = 5;
static constexpr int N = 4;
static constexpr int M = 0;

using vector_space = Vector_space<N, quaternon_type<N>>;

int main()
{
    const vector_space::base_type e0{1, 0, 0, 0};
    const vector_space::base_type e1{0, 1, 0, 0};
    const vector_space::base_type e2{0, 0, 1, 0};
    const vector_space::base_type e3{0, 0, 0, 1};

    vector_space d0( e0*2,    e1*(-1), e2*(-5), e3  ); // 4+1 +25 +1
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
    const auto grp = make_group(e0*3, e1*3,    e2*(-2), e3*5);
    group_comm<N, decltype (grp)> grp_p(grp);
    std::cout << grp_p << std::endl;

//    static_assert(product<quaternion_i, quaternion_i>::value == quaternion_zero::value, "NOT BINGO");
    constexpr auto tmp = mult<quaternion<'j'>, quaternion<'e'>>::value;
    auto const& ti = BOOST_CORE_TYPEID(tmp);
    std::cout << boost::core::demangled_name(ti) << std::endl;

    return 0;
}
