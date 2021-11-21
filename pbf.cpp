#include <iostream>

#include <vector>
#include <algorithm>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::VectorXd;

double W_poly6(Eigen::Ref<const VectorXd> r, double h) {
    double r_norm = r.norm();

    if (r_norm <= h) {
        return 315.0 / (64.0 * M_PI * std::pow(h, 9.0)) * std::pow(h*h - r_norm * r_norm, 3.0);
    } else {
        return 0.0;
    }
}

VectorXd dW_spiky(Eigen::Ref<const VectorXd> r, double h) {
    double r_norm = r.norm();

    double c = 0.0;
    if (r_norm <= h)
        c = 45.0 / (M_PI * std::pow(h, 6.0)) * std::pow(h - r_norm, 2.0) / r_norm;

    return c * r;
}

std::tuple<RowMatrixXd, RowMatrixXd> step(Eigen::Ref<const RowMatrixXd> x,
                                          Eigen::Ref<const RowMatrixXd> v0,
                                          Eigen::Ref<const RowMatrixXd> force,
                                          Eigen::Ref<const VectorXd> rho_init,
                                          double dt) {
    RowMatrixXd v = v0 + force * dt;
    RowMatrixXd p = x + v * dt;

    double h = 100.0;
    int N = x.rows();

    // find neighboring particles (should be optimized later)
    std::vector<std::vector<int>> neighbor;
    for(int i=0; i<N; i++) {
        std::vector<std::tuple<double, int>> dist_idx;
        for (int j=0; j<N; j++) {
            double distance = (p.row(i) - p.row(j)).norm();
            dist_idx.push_back(std::make_tuple(distance, j));
        }
        std::sort(dist_idx.begin(), dist_idx.end());

        std::vector<int> n;
        for(int k=0; k<8; k++) {
            n.push_back(std::get<1>(dist_idx[k]));
        }
        neighbor.push_back(n);
    }

    // the simulation loop (Algorithm 1)
    for (int loop = 0; loop < 5; loop++) {
        // calculate lambda
        std::vector<double> lam(N);
        for (int i=0; i<N; i++) {
            double rho = 0;
            double rho0 = rho_init(i);
            for (int j : neighbor[i]) {
                VectorXd r = p.row(i) - p.row(j);
                rho += W_poly6(r, h);
            }
            double C = rho / rho0 - 1.0;


            double Z = 0.0; // denominator for formula (11)
            // formula (8)
            for (int k : neighbor[i]) {
                if (k == i) {
                    VectorXd grad = VectorXd::Zero(3);
                    for (int j : neighbor[i]) {
                        if (i != j) {
                            VectorXd r = p.row(i) - p.row(j);
                            grad += dW_spiky(r, h) / rho0;
                        }
                    }
                    Z += std::pow(grad.norm(), 2.0);
                } else {
                    VectorXd r = p.row(i) - p.row(k);
                    VectorXd grad = -dW_spiky(r, h) / rho0;
                    Z += std::pow(grad.norm(), 2.0);
                }
            }
            lam[i] = C / (Z + 1e-5);
        }
        // calculate delta p_i
        RowMatrixXd delta_p(N, 3);
        VectorXd surface_normal(3);
        surface_normal << 0, 0, 1;

        for (int i=0; i<N; i++) {
            delta_p.row(i) = VectorXd::Zero(3);
            double rho0 = rho_init(i);
            for (int j : neighbor[i]) {
                if (i != j) {
                    VectorXd r = p.row(i) - p.row(j);
                    delta_p.row(i) += (lam[i]+lam[j]) * dW_spiky(r, h) / rho0;
                }
            }
            if (p(i, 2) < 0) {
                VectorXd p_i= p.row(i);
                VectorXd x_i= x.row(i);
                VectorXd q;
                if (x(i, 2) > 0) {
                    q = (x_i * x(i, 2) - p_i * p(i, 2)) / (x(i, 2) - p(i, 2));
                } else {
                    q = p_i - surface_normal.dot(p_i) * surface_normal;
                }
                double s = (p_i - q).dot(surface_normal);
                delta_p.row(i) -= s * surface_normal;
            }
        }

        for (int i=0; i<N; i++) {
            p.row(i) += delta_p.row(i);
        }
    }

    for (int i=0; i<N; i++) {
        v.row(i) = (p.row(i) - x.row(i)) / dt;
    }
    return std::make_tuple(p, v);
}

PYBIND11_MODULE(pbf, m) {
    m.def("W_poly6", &W_poly6);
    m.def("dW_spiky", &dW_spiky);
    m.def("step", &step);
}
