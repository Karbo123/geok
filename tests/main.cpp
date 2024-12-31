#include <iostream>

#include <Eigen/Dense>
#include <ostream>

#include "geok.hpp"

using namespace Eigen;

void test_1() {
    auto geok = GeoK::GeoK();

    MatrixXd A_ineq(3, 3);
    MatrixXd A_eq(0, 3);

    // -x  - y  +1  <= 0
    // -x  +3y  -3  <= 0
    //  x   -y  -1  <= 0
    A_ineq << -1, -1, 1, -1, 3, -3, 1, -1, -1;

    // expected vertex:
    // (0, 1)
    // (1, 0)
    // (3, 2)

    auto V = geok.H2V(A_ineq, A_eq);
    std::cout << "test_1:\n\n" << V << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_2() {
    auto geok = GeoK::GeoK();

    MatrixXd A_ineq(4, 3);
    MatrixXd A_eq(0, 3);

    // -2x -3y +5 <= 0
    // -7x +5y -3 <= 0
    //  2x  +y -4 <= 0
    //  3x -3y -1 <= 0
    A_ineq << -2, -3, 5, -7, 5, -3, 2, 1, -4, 3, -3, -1;

    // expected vertex (approximate):
    // A=(1, 2)
    // B=(0.51, 1.32)
    // C=(1.20, 0.86)
    // D=(1.44, 1.11)

    auto V = geok.H2V(A_ineq, A_eq);
    std::cout << "test_2:\n\n" << V << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_3() {
    auto geok = GeoK::GeoK();

    MatrixXd A_ineq(4, 3);
    MatrixXd A_eq(1, 3);

    // -2x -3y +5 <= 0
    // -7x +5y -3 <= 0
    //  2x  +y -4 <= 0
    //  3x -3y -1 <= 0
    A_ineq << -2, -3, 5, -7, 5, -3, 2, 1, -4, 3, -3, -1;

    // x -2y +2 = 0
    A_eq << 1, -2, 2;

    // expected vertex (approximate):
    // (0.57, 1.28)
    // (1.20, 1.60)

    auto V = geok.H2V(A_ineq, A_eq);
    std::cout << "test_3:\n\n" << V << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_4() {
    auto geok = GeoK::GeoK();

    MatrixXd V(4, 2);
    V << 1, 1, 1, 2, 2, 1, 2, 2;

    // expected:
    //   a square with min=(1,1) and max=(2,2)

    auto [H_ineq, H_eq] = geok.V2H(V);
    std::cout << "test_4:\n\nH_ineq =\n" << H_ineq << "\n\nH_eq=\n" << H_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_5() {
    auto geok = GeoK::GeoK();

    MatrixXd V(5, 2);
    V << -1, 2, 3, 5, 7, 2, 1, 3, 4, 5;

    // expected:
    //   (1, 3) is inactive

    auto [H_ineq, H_eq] = geok.V2H(V);
    std::cout << "test_5:\n\nH_ineq =\n" << H_ineq << "\n\nH_eq=\n" << H_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_6() {
    auto geok = GeoK::GeoK();

    MatrixXd V(2, 2);
    V << -1, 2, 3, 5;

    // expected:
    //   an equality, and two inequalities

    auto [H_ineq, H_eq] = geok.V2H(V);
    std::cout << "test_6:\n\nH_ineq =\n" << H_ineq << "\n\nH_eq=\n" << H_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_7() {
    auto geok = GeoK::GeoK();

    MatrixXd V(1, 2);
    V << -1, 2;

    // expected:
    //   two equalities
    // NOTE:
    //   result has one redundancy: "-1 <= 0"

    auto [H_ineq, H_eq] = geok.V2H(V);
    std::cout << "test_7:\n\nH_ineq =\n" << H_ineq << "\n\nH_eq=\n" << H_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_8() {
    auto geok = GeoK::GeoK();

    MatrixXd V(1, 2);
    V << -1, 2;

    // expected:
    //   no solution

    auto [H_ineq, H_eq] = geok.V2H(V);
    auto ball = geok.inscribed(H_ineq);

    std::cout << "test_8:\n\n";
    if (ball.has_value()) {
        std::cout << "ball->ctr = \n" << ball->ctr << "\n\n";
        std::cout << "ball->rad = \n" << ball->rad << "\n\n";
    } else {
        std::cout << "no solution" << "\n\n";
    }

    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_9() {
    auto geok = GeoK::GeoK();

    MatrixXd V(3, 2);
    V << 1, 0, 0, 1, 3, 2;

    // expected:
    //   (1, 0.763932)   with r = 0.540182

    auto [H_ineq, H_eq] = geok.V2H(V);
    auto ball = geok.inscribed(H_ineq);

    std::cout << "test_9:\n\n";
    if (ball.has_value()) {
        std::cout << "ball->ctr = \n" << ball->ctr << "\n\n";
        std::cout << "ball->rad = \n" << ball->rad << "\n\n";
    } else {
        std::cout << "no solution" << "\n\n";
    }

    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_10() {
    auto geok = GeoK::GeoK();

    MatrixXd H_ineq(3, 3);

    H_ineq << 1, 1, -1, -1, -1, 1, 1, -1, 0;

    auto ball = geok.inscribed(H_ineq);

    // expected:
    //   (0.5, 0.5)  with r = 0

    std::cout << "test_10:\n\n";
    if (ball.has_value()) {
        std::cout << "ball->ctr = \n" << ball->ctr << "\n\n";
        std::cout << "ball->rad = \n" << ball->rad << "\n\n";
    } else {
        std::cout << "no solution" << "\n\n";
    }

    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_11() {
    auto geok = GeoK::GeoK();

    MatrixXd H_ineq(3, 3);

    H_ineq << 2, 3, 5, -3, -1, -9, 1, -1, -5;

    auto ball = geok.inscribed(H_ineq);

    // expected:
    //

    std::cout << "test_11:\n\n";
    if (ball.has_value()) {
        std::cout << "ball->ctr = \n" << ball->ctr << "\n\n";
        std::cout << "ball->rad = \n" << ball->rad << "\n\n";
    } else {
        std::cout << "no solution" << "\n\n";
    }

    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_12() {
    auto geok = GeoK::GeoK();

    MatrixXd A_ineq(3, 3);
    MatrixXd A_eq(0, 3);

    // -x  - y  +1  <= 0
    // -x  +3y  -3  <= 0
    //  x   -y  -1  <= 0
    A_ineq << -1, -1, 1, -1, 3, -3, 1, -1, -1;

    // expected vertex:
    // (0, 1)
    // (1, 0)
    // (3, 2)

    auto [V, Conn] = geok.H2V<true>(A_ineq, A_eq);
    std::cout << "test_12:\n\n" << V << "\n\nConn:\n" << Conn << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_13() {
    auto geok = GeoK::GeoK();

    MatrixXd A_ineq(4, 3);
    MatrixXd A_eq(1, 3);

    // -2x -3y +5 <= 0
    // -7x +5y -3 <= 0
    //  2x  +y -4 <= 0
    //  3x -3y -1 <= 0
    A_ineq << -2, -3, 5, -7, 5, -3, 2, 1, -4, 3, -3, -1;

    // x -2y +2 = 0
    A_eq << 1, -2, 2;

    // expected vertex (approximate):
    // (0.57, 1.28)
    // (1.20, 1.60)

    auto [V, Conn] = geok.H2V<true>(A_ineq, A_eq);
    std::cout << "test_13:\n\n" << V << "\n\nConn:\n" << Conn << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<int>> &vec) {
    os << "[\n";
    for (const auto &row : vec) {
        os << "  [ ";
        for (const auto &elem : row) {
            os << elem << " ";
        }
        os << "]\n";
    }
    os << "]";
    return os;
}

void test_14() {
    auto geok = GeoK::GeoK();

    MatrixXd V(5, 2);
    V << -1, 2, 3, 5, 7, 2, 1, 3, 4, 5;

    // expected:
    //   (1, 3) is inactive

    auto [H_ineq, H_eq, Conn_ineq, Conn_eq] = geok.V2H<true>(V);
    std::cout << "test_14:\n\nH_ineq =\n"
              << H_ineq << "\n\nH_eq=\n"
              << H_eq << "\n\nConn_ineq:\n"
              << Conn_ineq << "\n\nConn_eq:\n"
              << Conn_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

void test_15() {
    auto geok = GeoK::GeoK();

    MatrixXd V(6, 4);
    V << 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1;

    // expected:
    //   3D CUBE

    auto [H_ineq, H_eq, Conn_ineq, Conn_eq] = geok.V2H<true>(V);
    std::cout << "test_15:\n\nH_ineq =\n"
              << H_ineq << "\n\nH_eq=\n"
              << H_eq << "\n\nConn_ineq:\n"
              << Conn_ineq << "\n\nConn_eq:\n"
              << Conn_eq << "\n\n";
    std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = =\n\n";
}

int main() {
    // H2V
    test_1();
    test_2();
    test_3();
    // V2H
    test_4();
    test_5();
    test_6();
    test_7();
    // inscribed
    test_8();
    test_9();
    test_10();
    test_11();
    // connectivity for H2V
    test_12();
    test_13();
    // connectivity for V2H
    test_14();
    test_15();
}