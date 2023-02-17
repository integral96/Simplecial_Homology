#include <iostream>

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"
static constexpr int N = 3;
using vector_space = Vector_space<4, int>;
int main()
{
    vector_space d0(0, 0, 0, 0);
    vector_space d1(1, 0, 0, 0);
    vector_space d2(0, 1, 0, 0);
    vector_space d3(0, 0, 1, 0);
    vector_space d4(0, 0, 0, 1);


    Simplex<4, vector_space> simplex4D_0(d0*2, d1*(-2), d2*8, d3*(-3), d4*(-2));
    Simplex<4, vector_space> simplex4D_1(d0*4, d2*(-6), d1*8, d3*(-3), d4*8);

    Simplex<3, vector_space> simplex3D_0(d0*2, d1*(-2), d2*8, d3*(-3));
    Simplex<3, vector_space> simplex3D_1(d0*4, d2*(-6), d1*8, d3*(-3));
    Simplex<2, vector_space> simplex2D_0(d0*2, d1*(-2), d2*8);
    Simplex<2, vector_space> simplex2D_1(d0*2, d1*(2), d2*5);
    Simplex<1, vector_space> simplex1D_0(d0*2, d1*(-2));

    Simplex<0, vector_space> simplex0D_0(d1*2);
    Simplex<0, vector_space> simplex0D_1(d2*2);
    Simplex<0, vector_space> simplex0D_2(d3*2);


    std::cout << "START MAKE COMPLEX" << std::endl;
    auto complex = make_complex(simplex4D_0, simplex4D_1, simplex3D_0, simplex3D_1, simplex2D_0, simplex2D_1, simplex1D_0, simplex0D_1, simplex0D_2);
    std::cout << "END MAKE COMPLEX" << std::endl;
    std::cout << "GET COMPLEX DIMENTION = " << N << std::endl;
    Simplicial_complex complex_simpl(complex);
    std::cout << complex_simpl << std::endl;
    const auto compl_N = complex_simpl.get_complex<N>();
    std::cout << "END COMPLEX DIMENTION = " << N << std::endl;
    std::cout << "GET COMPLEX BOUNDARY: " << std::endl;
    boundary<N, decltype(compl_N), vector_space> bndr_N(compl_N);
    const auto B_N = bndr_N.get();
    hana::for_each(B_N, [&](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });
    boundary<N - 1, decltype(B_N), vector_space> bndr_N_1(B_N);
    const auto B_N_1 = bndr_N_1.get();
    hana::for_each(B_N_1, [&](auto x) {
        std::cout << "Simplex dim = " << decltype(x)::dim << "; value = " << x << '\n';
    });

//    std::cout << "dif(N)*dif(N - 1) = " << tmp3 << std::endl;
    return 0;
}
