#include <iostream>

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"
static constexpr int N = 4;
static constexpr int M = 3;
using vector_space = Vector_space<4, int>;
int main()
{
    vector_space d0( 2,  4,  3, -2);
    vector_space d1(-1,  3, -2,  5);
    vector_space d2( 1, -5,  3,  2);
    vector_space d3(-1, -6,  2,  4);
    vector_space d4(-2,  1, -3,  2);

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


    std::cout << "Создание комплекса до размерностей " << N << "-х" << std::endl;
    auto complex = make_complex(simplex5D_0, simplex5D_1, simplex5D_2,
                                simplex4D_0, simplex4D_1, simplex3D_0,
                                simplex3D_1, simplex3D_2, simplex2D_0,
                                simplex2D_1, simplex1D_0, simplex1D_1,
                                simplex0D_0, simplex0D_1, simplex0D_2);
    std::cout << "Создание цепного комплекса из аблевых групп." << std::endl;
    std::cout << "GET COMPLEX DIMENTION = " << N << std::endl;
    Simplicial_complex complex_simpl(complex);
    auto chain = complex_simpl.chain_unit<N>();


    std::cout << "END COMPLEX DIMENTION = " << complex_simpl << N << std::endl;
    std::cout << "GET COMPLEX BOUNDARY: " << std::endl;
    auto chain_N  = boundary<N>(chain); //Kernel
    auto chain_N_1  = boundary<N - 1>(chain_N);
    std::cout << "Взятие граничного оператора: " << std::endl;
    hana::for_each(chain_N, [](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
    std::cout << "Двойное взятие граничного оператора должно быть равно нулую: " << std::endl;
    hana::for_each(chain_N_1, [&](auto x) {
        BOOST_ASSERT_MSG(x.empty(), "Двойное взятие граничного оператора должно быть равно нулую.");
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
    std::cout << "Определим цепь для циклов: " << std::endl;
    auto chain_Im = complex_simpl.chain_unit<M + 1>();
    hana::for_each(chain_Im, [](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
//    auto Kern_N  = Kernel<N>(chain_1);
    std::cout << "Определим границы для циклов: " << M + 1 << std::endl;
    auto Img_N1  = Image<M + 1>(chain_Im);
    hana::for_each(Img_N1, [](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
    std::cout << "Определим циклы: " << M << std::endl;
    auto chain_Ker = complex_simpl.chain_unit<M>();
    auto Ker_N1  = Kernel<M>(chain_Ker);
    hana::for_each(Ker_N1, [](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
    auto homology = quotient(Ker_N1, Img_N1);
    std::cout << "Группа гомологий: K = " << M << std::endl;
    for(const  auto& x : homology) {
        std::cout << "<";
        for(const  auto& y : x) {
            std::cout << y << " ";
        }
        std::cout << ">\n";
    }

    return 0;
}
