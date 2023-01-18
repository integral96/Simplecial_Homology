#include <iostream>

#include "Simplicial_complex.hpp"
#include "Vector_space.hpp"

int main()
{
    Vector_space<3, int> a0(0, 0, 0);
    Vector_space<3, int> a1(1, 0, 0);
    Vector_space<3, int> a2(0, 1, 0);
    Vector_space<3, int> a3(0, 0, 1);

    Vector_space<2, int> b0(0, 0);
    Vector_space<2, int> b1(1, 0);
    Vector_space<2, int> b2(0, 1);
//    std::cout << a0 << std::endl;
//    Simplex<6, int> simplex(1, 2, 3, 4, 5, 6, 7);
    Simplex<3, Vector_space<3, int>> simplex3D(a0*2, a1*(-2), a2*8, a3*(-3));
    Simplex<3, Vector_space<3, int>> simplex3D1(a0*4, a2*(-6), a1*8, a3*(-3));
    Simplex<2, Vector_space<2, int>> simplex2D(b0*2, b1*(-2), b2*8);
    auto tmp = (simplex3D - simplex3D);
    auto tmp1 = (simplex3D - simplex3D1);
    std::cout << "CHAIN: " << tmp1 << std::endl;
//    simplex3D.boundary();
//    simplex2D.boundary();
    auto tmp3 = simplex3D.boundary()*simplex2D.boundary();
    std::cout << "dif(N)*dif(N - 1) = " << tmp3 << std::endl;
    return 0;
}
