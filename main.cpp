#include <iostream>

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"
static constexpr int N = 3;
int main()
{
    Vector_space<4, int> d0(0, 0, 0, 0);
    Vector_space<4, int> d1(1, 0, 0, 0);
    Vector_space<4, int> d2(0, 1, 0, 0);
    Vector_space<4, int> d3(0, 0, 1, 0);
    Vector_space<4, int> d4(0, 0, 0, 1);

    Vector_space<3, int> a0(0, 0, 0);
    Vector_space<3, int> a1(1, 0, 0);
    Vector_space<3, int> a2(0, 1, 0);
    Vector_space<3, int> a3(0, 0, 1);

    Vector_space<2, int> b0(0, 0);
    Vector_space<2, int> b1(1, 0);
    Vector_space<2, int> b2(0, 1);

    Vector_space<1, int> c0(1);
    Vector_space<1, int> c1(0);

    Simplex<4, Vector_space<4, int>> simplex4D_0(d0*2, d1*(-2), d2*8, d3*(-3), d4*(-2));
    Simplex<4, Vector_space<4, int>> simplex4D_1(d0*4, d2*(-6), d1*8, d3*(-3), d4*8);

    Simplex<3, Vector_space<3, int>> simplex3D_0(a0*2, a1*(-2), a2*8, a3*(-3));
    Simplex<3, Vector_space<3, int>> simplex3D_1(a0*4, a2*(-6), a1*8, a3*(-3));
    Simplex<2, Vector_space<2, int>> simplex2D_0(b0*2, b1*(-2), b2*8);
    Simplex<2, Vector_space<2, int>> simplex2D_1(b0*2, b1*(2), b2*5);
    Simplex<1, Vector_space<1, int>> simplex1D_0(c0*2, c1*(-2));

    Simplex<0, Vector_space<1, int>> simplex0D_0(c0*2);
    Simplex<0, Vector_space<1, int>> simplex0D_1(c0*2);
    Simplex<0, Vector_space<1, int>> simplex0D_2(c0*2);


    std::cout << "START MAKE COMPLEX" << std::endl;
    auto complex = make_complex(simplex4D_0, simplex4D_1, simplex3D_0, simplex3D_1, simplex2D_0, simplex2D_1, simplex1D_0, simplex0D_1, simplex0D_2);
    std::cout << "END MAKE COMPLEX" << std::endl;
    std::cout << "GET COMPLEX DIMENTION = " << N << std::endl;
    Simplicial_complex complex_simpl(complex);
    std::cout << complex_simpl << std::endl;
    const auto compl_N = complex_simpl.get_complex<N>();
    std::cout << "END COMPLEX DIMENTION = " << N << std::endl;
    std::cout << "GET COMPLEX BOUNDARY: " << std::endl;
    boundary<N, decltype(compl_N)> bndr(compl_N);
//    simplex3D.boundary();
//    simplex2D.boundary();
//    auto tmp3 = simplex3D.boundary()*simplex2D.boundary();
//    std::cout << "dif(N)*dif(N - 1) = " << tmp3 << std::endl;
    return 0;
}
