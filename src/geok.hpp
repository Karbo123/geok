//
// each row is for one constraint
//     A_ineq * homo(x) <= 0
//     A_eq   * homo(x) == 0
// assume: double precision

#include <cassert>
#include <optional>

#include <Eigen/Dense>
#include <tuple>
#include <vector>

namespace GeoK {

#include "cddlib/cdd.h"
#include "cddlib/setoper.h"

using namespace Eigen;

struct Ball { // inscribed ball
    double rad;
    VectorXd ctr;
    Ball(int dim) : rad(0.0), ctr(dim) {}
};

struct GeoK {

    GeoK() {
        // some global variables
        dd_set_global_constants();
    }

    inline double read(mytype x) { return x[0]; }

    void copy_in_(dd_MatrixPtr M, const MatrixXd &mat, int row_start = 0, bool is_eq = false) {
        const int num_rows = mat.rows();
        const int num_cols = mat.cols();

        // copy one by one
        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < num_cols; ++j) {
                dd_set_d(M->matrix[i + row_start][j], mat(i, j));
            }
        }

        // if they are equalities
        if (is_eq) {
            for (int n = 1; n <= num_rows; ++n) {
                set_addelem(M->linset, row_start + n);
            }
        }
    }

    void copy_in_(dd_MatrixPtr M, const VectorXd &vec) {
        const int num_elems = vec.size();

        // copy one by one
        for (int i = 0; i < num_elems; ++i) {
            dd_set_d(M->rowvec[i], vec[i]);
        }
    }

    MatrixXd copy_out(dd_MatrixPtr M) {

        const int num_rows = M->rowsize;
        const int num_cols = M->colsize;

        MatrixXd out(num_rows, num_cols);
        dd_Amatrix ptr = M->matrix;

        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < num_cols; ++j) {
                out(i, j) = read(ptr[i][j]);
            }
        }
        return out;
    }

    MatrixXd shift_to_first(const MatrixXd &mat) {
        // shift the last column to the first
        MatrixXd newMat(mat.rows(), mat.cols());
        newMat.col(0) = mat.col(mat.cols() - 1);
        newMat.rightCols(mat.cols() - 1) = mat.leftCols(mat.cols() - 1);
        return newMat;
    }

    MatrixXd shift_to_last(const MatrixXd &mat) {
        // shift the first column to the last
        MatrixXd newMat(mat.rows(), mat.cols());
        newMat.col(mat.cols() - 1) = mat.col(0);
        newMat.leftCols(mat.cols() - 1) = mat.rightCols(mat.cols() - 1);
        return newMat;
    }

    MatrixXd pad_first_one(const MatrixXd &mat) {
        // pad to the first with all ones column
        MatrixXd newMat(mat.rows(), 1 + mat.cols());
        newMat.rightCols(mat.cols()) = mat;
        newMat.col(0).fill(1);
        return newMat;
    }

    MatrixXd hconcat(const MatrixXd &A, const MatrixXd &B) {
        Eigen::MatrixXd C(A.rows(), A.cols() + B.cols());
        C << A, B;
        return C;
    }

    template <bool return_connectivity = false> auto H2V(const MatrixXd &A_ineq, const MatrixXd &A_eq) {
        const int dim = A_ineq.cols() - 1;    // dimension of space
        const int num_ineq = A_ineq.rows();   // number of inequalities
        const int num_eq = A_eq.rows();       // number of equalities
        assert(A_ineq.cols() == A_eq.cols()); // dimensions must match

        dd_PolyhedraPtr poly;
        dd_MatrixPtr H, V;
        dd_ErrorType err;

        dd_rowrange m = num_ineq + num_eq;
        dd_colrange d = dim + 1;
        H = dd_CreateMatrix(m, d);

        //
        MatrixXd A_ineq_shifted = shift_to_first(A_ineq) * (-1.0);
        MatrixXd A_eq_shifted = shift_to_first(A_eq) * (-1.0);

        // copy constraint data
        copy_in_(H, A_ineq_shifted);
        copy_in_(H, A_eq_shifted, num_ineq, true);

        // this is H-repre
        H->representation = dd_Inequality;

        // compute the polyhedron from H-repre
        poly = dd_DDMatrix2Poly(H, &err);
        assert(err == dd_NoError);

        // copy result out as a matrix
        V = dd_CopyGenerators(poly);

        // convert to Eigen matrix
        MatrixXd verts = copy_out(V).rightCols(dim);

        MatrixXd conn(0, dim);
        if constexpr (return_connectivity) {
            const int num_verts = verts.rows();
            conn = {num_verts, dim};

            dd_SetFamilyPtr I = dd_CopyIncidence(poly);
            assert(I != NULL);
            assert(I->famsize == num_verts); // number of sets

            for (int i = 0; i < num_verts; i++) { // for each vertex
                int card = set_card(I->set[i]);   // the actual number in the set
                assert(card == dim);
                assert(I->setsize >= 2 * card); // a set stores two things, but we only use one ???
                auto set = I->set[i];
                int k = 0;
                assert(I->setsize >= set[0]); // the maximal allocated size >= capacity size
                for (int j = 0; j < set[0]; ++j) {
                    if (set_member(j + 1, set)) {
                        // NOTE : if the "value >= num_ineq", this means "value - num_ineq" is the index for "A_eq"
                        conn(i, k++) = j;
                    }
                }
                assert(k == dim);
                assert(k == card);
            }

            dd_FreeSetFamily(I);
        }

        // free memory
        dd_FreeMatrix(H);
        dd_FreeMatrix(V);
        dd_FreePolyhedra(poly);

        if constexpr (!return_connectivity) {
            return verts;
        } else {
            return std::make_tuple(verts, conn);
        }
    }

    template <bool return_connectivity = false> auto V2H(const MatrixXd &pts) {
        const int num_pts = pts.rows(); // number of inequalities
        const int dim = pts.cols();     // dimension of space

        if (num_pts == 0) {
            if constexpr (!return_connectivity) {
                return std::make_tuple(MatrixXd(0, dim + 1), MatrixXd(0, dim + 1));
            } else {
                return std::make_tuple(MatrixXd(0, dim + 1), MatrixXd(0, dim + 1), std::vector<std::vector<int>>{},
                                       std::vector<std::vector<int>>{});
            }
        }

        dd_PolyhedraPtr poly;
        dd_MatrixPtr V, H;
        dd_ErrorType err;

        dd_rowrange m = num_pts;
        dd_colrange d = 1 + dim;
        V = dd_CreateMatrix(m, d);

        // copy constraint data
        copy_in_(V, pad_first_one(pts));

        // this is V-repre
        V->representation = dd_Generator;

        // compute the polyhedron from V-repre
        poly = dd_DDMatrix2Poly(V, &err);
        assert(err == dd_NoError);

        // copy result out as a matrix
        H = dd_CopyInequalities(poly);

        // convert to Eigen matrix
        MatrixXd hs = copy_out(H);
        hs = shift_to_last(hs) * (-1.0);

        // separate equalities and inequalities
        const int num_cons = hs.rows();
        std::vector<int> idx_ineq, idx_eq;
        for (int i = 1; i <= num_cons; ++i) {
            int is_eq = set_member(i, H->linset);
            (is_eq ? idx_eq : idx_ineq).push_back(i - 1);
        }
        MatrixXd A_ineq = hs(idx_ineq, placeholders::all);
        MatrixXd A_eq = hs(idx_eq, placeholders::all);

        // find vertices that form the vertex
        std::vector<std::vector<int>> conn_ineq, conn_eq;
        if constexpr (return_connectivity) {
            conn_ineq.resize(A_ineq.rows());
            conn_eq.resize(A_eq.rows());

            dd_SetFamilyPtr I = dd_CopyIncidence(poly);
            assert(I != NULL);
            assert(I->famsize == num_cons); // number of sets

            dd_WriteIncidence(stdout, poly);

            int k_ineq = 0, k_eq = 0;
            for (int i = 0; i < num_cons; i++) { // for each vertex
                int is_eq = set_member(i + 1, H->linset);
                auto &conn = is_eq ? conn_eq : conn_ineq;
                auto &k = is_eq ? k_eq : k_ineq;

                int card = set_card(I->set[i]); // the actual number in the set
                if (I->setsize < 2 * card) {
                    continue; // we don't know why this is happening in this case
                };
                auto set = I->set[i];
                assert(I->setsize >= set[0]); // the maximal allocated size >= capacity size

                for (int j = 0; j < set[0]; ++j) {
                    if (set_member(j + 1, set)) {
                        conn[i].push_back(j);
                    }
                }
                assert(conn[i].size() == card);
                k++;
            }
            assert(k_ineq == A_ineq.rows());
            assert(k_eq <= A_eq.rows()); // because we skip some cases

            dd_FreeSetFamily(I);
        }

        // free memory
        dd_FreeMatrix(V);
        dd_FreeMatrix(H);
        dd_FreePolyhedra(poly);

        if constexpr (!return_connectivity) {
            return std::make_tuple(A_ineq, A_eq);
        } else {
            return std::make_tuple(A_ineq, A_eq, conn_ineq, conn_eq);
        }
    }

    std::optional<Ball> inscribed(const MatrixXd &A_ineq) {
        // NOTE : not allow any equality as input, otherwise the radius must be zero

        const int num_ineq = A_ineq.rows(); // number of inequalities
        const int dim = A_ineq.cols() - 1;  // dimension of space

        dd_MatrixPtr H;
        dd_rowrange m = num_ineq;
        dd_colrange n = dim + 2;
        H = dd_CreateMatrix(m, n);

        // copy constraint data
        MatrixXd A_ineq_aug = hconcat(shift_to_first(A_ineq), A_ineq.leftCols(dim).rowwise().norm()) * (-1.0);
        copy_in_(H, A_ineq_aug);

        // coefficients of the objective
        // maximize [obj_vec^T @ x]
        //  x = [...coords; radius]     length = dim + 1
        VectorXd obj_vec(dim + 2);
        obj_vec.fill(0);
        obj_vec[dim + 1] = 1;
        copy_in_(H, obj_vec);
        H->objective = dd_LPmax;

        //
        dd_ErrorType err = dd_NoError;
        dd_LPPtr lp;
        lp = dd_Matrix2LP(H, &err);
        assert(err == dd_NoError);
        assert(lp != NULL);

        // solve
        dd_LPSolve(lp, dd_DualSimplex, &err);
        assert(err == dd_NoError);

        // write
        Ball ball(dim);
        if (lp->LPS == dd_Optimal) {
            for (int j = 1; j < lp->d - 1; j++) {
                ball.ctr[j - 1] = read(lp->sol[j]);
            }
            // NOTE : radius might be zero
            ball.rad = read(lp->sol[lp->d - 1]);
        } else {
            return std::nullopt;
        }

        //
        dd_FreeLPData(lp);
        dd_FreeMatrix(H);

        return ball;
    }

    ~GeoK() { // release global variables
        dd_free_global_constants();
    }
};

} // namespace GeoK
