#include <Eigen/Sparse>
#include "cgeneric_defs.h"


extern "C" {

double nChoosek( int n, int k );

void compute_Q_dim1(
    double kappa, double sigma, double gamma, double rho, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list
);
void compute_Q_dim2(
    double kappa, double sigma, double gamma, double rho_1, double rho_2, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list, 
    inla_cgeneric_smat_tp*** M2list2
);

} // extern "C"

// Define a type alias for SparseMatrix with RowMajor storage
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrixRowMajor;

// Helper function to convert inla_cgeneric_smat_tp to Eigen::SparseMatrix with RowMajor storage
SparseMatrixRowMajor convertInlaToEigen(const inla_cgeneric_smat_tp* inlaMat) {
    if (!inlaMat || inlaMat->n <= 0) {
        throw std::invalid_argument("Invalid INLA sparse matrix");
    }
    SparseMatrixRowMajor eigenMat(inlaMat->nrow, inlaMat->ncol);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(inlaMat->n);

    for (int idx = 0; idx < inlaMat->n; ++idx) {
        triplets.emplace_back(inlaMat->i[idx], inlaMat->j[idx], inlaMat->x[idx]);
    }

    eigenMat.setFromTriplets(triplets.begin(), triplets.end());
    return eigenMat;
}

// Function to create matrix L based on Glist, kappa, and n, using RowMajor storage
SparseMatrixRowMajor make_L(int n, double kappa, inla_cgeneric_smat_tp **Glist) {
    SparseMatrixRowMajor L = convertInlaToEigen(Glist[0]);
    L.setZero();  // Start with L as a zero matrix of the same size

    for (int k = 0; k <= n; k++) {
        double coeff = nChoosek(n, k) * pow(kappa, 2 * (n - k));
        SparseMatrixRowMajor Gk = convertInlaToEigen(Glist[k]);

        L += coeff * Gk;
    }

    return L;
}

// Main function to compute Q for the 1D case (d = 1)
void compute_Q_dim1(
    double kappa, double sigma, double gamma, double rho, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list
) {
    // Initialize Q with the first term using RowMajor storage
    SparseMatrixRowMajor Q = make_L(beta, kappa, Gtlist) + 2 * gamma * make_L(beta + alpha, kappa, B0list);



    // Loop over alpha to add the contributions from each term
    for (int k = 0; k <= alpha; ++k) {
        double gamma_sq_rho = gamma * gamma * nChoosek(alpha, k) * pow(rho, 2 * k);

        // Compute the gamma^2 term
        SparseMatrixRowMajor Ct_sum = make_L(beta + 2 * (alpha - k), kappa, &Ctlist[k]);

        // Add this term to Q
        Q += gamma_sq_rho * Ct_sum;

        // Calculate the M2 term
        SparseMatrixRowMajor M2_term = make_L(beta + alpha - k, kappa, M2list[k]);

        // Add the M2 term to Q
        double factor = -0.5 * pow(-1, k / 2) * gamma * nChoosek(alpha, k) * (1 - pow(-1, k)) * pow(rho, k);

        // Ensure that M2_term and its transpose have the same storage order
        SparseMatrixRowMajor M2_term_transpose = M2_term.transpose();

        Q += factor * (M2_term + M2_term_transpose);
    }

    Q /= sigma * sigma;

    // Extract the values from Q into the result array
    Eigen::SparseMatrix<double> Q_triang(Q.rows(), Q.rows());
    Q_triang = Q.triangularView<Eigen::Lower>();    
    int count = 0;
    for (int k = 0; k < Q_triang.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang, k); it; ++it) {
            result[count] = it.value();
            count++;
        }
    }
}

// Main function to compute Q for the 2D case (d = 2)
void compute_Q_dim2(
    double kappa, double sigma, double gamma, double rho_1, double rho_2, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list,
    inla_cgeneric_smat_tp*** M2list2
) {
    // Initialize Q with the first term using RowMajor storage
    SparseMatrixRowMajor Q = make_L(beta, kappa, Gtlist) + 2 * gamma * make_L(beta + alpha, kappa, B0list);

    for (int k = 0; k <= alpha; ++k) {
        double gamma_sq_rho = gamma * gamma * nChoosek(alpha, k) * pow(rho_1, 2 * k);

        SparseMatrixRowMajor Ct_sum = make_L(beta + 2 * (alpha - k), kappa, &Ctlist[k]);

        Q += gamma_sq_rho * Ct_sum;

        // Calculate the M2 terms for 2D (M2list and M2list2)
        SparseMatrixRowMajor M2x_term = make_L(beta + alpha - k, kappa, M2list[k]);
        SparseMatrixRowMajor M2y_term = make_L(beta + alpha - k, kappa, M2list2[k]);

        // Add the M2x and M2y terms to Q
        double factor = -0.5 * gamma * nChoosek(alpha, k) * (1 - pow(-1, k));

        // Ensure that the storage orders match
        SparseMatrixRowMajor M2x_term_transpose = M2x_term.transpose();
        SparseMatrixRowMajor M2y_term_transpose = M2y_term.transpose();

        Q += factor * pow(rho_1, k) * (M2x_term + M2x_term_transpose);
        Q += factor * pow(rho_2, k) * (M2y_term + M2y_term_transpose);
    }

    Q /= sigma * sigma;

    // Extract the values from Q into the result array
    Eigen::SparseMatrix<double> Q_triang(Q.rows(), Q.rows());
    Q_triang = Q.triangularView<Eigen::Lower>();    

    int count = 0;
    for (int k = 0; k < Q_triang.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang, k); it; ++it) {
            result[count] = it.value();
            count++;
        }
    }
}
