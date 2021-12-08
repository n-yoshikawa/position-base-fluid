#include <iostream>

#include <vector>
#include <algorithm>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <pybind11/eigen.h>
#include <eigen3/Eigen/Dense>

namespace py = pybind11;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorXd = Eigen::VectorXd;


// clip function from https://stackoverflow.com/questions/9323903/most-efficient-elegant-way-to-clip-a-number
float clip(float n, float lower, float upper) {
  return std::max(lower, std::min(n, upper));
}

// cross product of two columns vectors from https://en.wikipedia.org/wiki/Cross_product
VectorXd cross_prod(VectorXd a, VectorXd b) {
    VectorXd s = VectorXd::Zero(3);
    s(0) = a(1)*b(2) - a(2)*b(1);
    s(1) = a(2)*b(0) - a(0)*b(2);
    s(2) = a(0)*b(1) - a(1)*b(0);
    return s;
}

// poly6 kernel from https://matthias-research.github.io/pages/publications/sca03.pdf
double W_poly6(Eigen::Ref<const VectorXd> r, double h) {
    double r_norm = r.norm();
    if (r_norm <= h) {
        return 315.0 / (64.0 * M_PI * std::pow(h, 9.0)) * std::pow(h*h - r_norm * r_norm, 3.0);
    } else {
        return 0.0;
    }
}

// spiky kernel from https://matthias-research.github.io/pages/publications/sca03.pdf
// gradient of kernel computer below wrt r
VectorXd dW_spiky(Eigen::Ref<const VectorXd> r, double h) {
    double r_norm = r.norm();
    double c = 0.0;
    
    if (r_norm <= 0.000001) { // stability issues with small values of r_nor
        c = 0.0;
    } else if (r_norm <= h) {
        c = -(45.0 / (M_PI * std::pow(h, 6.0))) * std::pow(h - r_norm, 2.0) / r_norm;
    }
    return c * r;
}

std::tuple<RowMatrixXd, RowMatrixXd> step(Eigen::Ref<const RowMatrixXd> x,
                                          Eigen::Ref<const RowMatrixXd> v0,
                                          Eigen::Ref<const RowMatrixXd> force,
                                          Eigen::Ref<const VectorXd> rho_init,
                                          double dt,
					  double h, 
					  double boundary,
					  double epsilon) {
    RowMatrixXd v = v0 + force * dt;
    RowMatrixXd p = x + v * dt;

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
    std::vector<double> rho(N);
    // the simulation loop (Algorithm 1)
    for (int loop = 0; loop < 5; loop++) {
        // calculate lambda
        std::vector<double> lam(N);
        for (int i=0; i<N; i++) {
            double rho_i = 0;
            double rho0 = rho_init(i);
            for (int j : neighbor[i]) {
                VectorXd r = p.row(i) - p.row(j);
                rho_i += W_poly6(r, h);
            }
	    rho[i] = rho_i;
            double C = rho[i] / rho0 - 1.0;
    
            double Z = 0.0; // denominator for formula (11)
            // formula (8)
            for (int k : neighbor[i]) {
                if (k == i) {
                    VectorXd grad = VectorXd::Zero(3);
                    for (int j : neighbor[i]) {
                        if (i != j) {
                            VectorXd r = p.row(i) - p.row(j);
                            grad += dW_spiky(r, h);
                        }
                    }
		    grad = grad / rho0;
                    Z += std::pow(grad.norm(), 2.0);

                } else {
                    VectorXd r = p.row(i) - p.row(k);
                    VectorXd grad = -dW_spiky(r, h);
                    grad = grad / rho0;
		    Z += std::pow(grad.norm(), 2.0);
                }
            }
            lam[i] = - C / (Z + epsilon); 
           
	}
        // calculate delta p_i
        RowMatrixXd delta_p(N, 3);
        for (int i=0; i<N; i++) {
            delta_p.row(i) = VectorXd::Zero(3);
            double rho0 = rho_init(i);
            for (int j : neighbor[i]) {
                if (i != j) {
                    VectorXd r = p.row(i) - p.row(j);
		    VectorXd corr(3); 
		    corr << 0.1*h, 0.1*h, 0.1*h;
		    double s_corr = -0.0001 * std::pow(W_poly6(r, h) / W_poly6(corr, h),4);
                    //double s_corr = 0;
		    delta_p.row(i) += (lam[i]+lam[j] + s_corr) * dW_spiky(r, h) / rho0;
                }
            }
        }

        for (int i=0; i<N; i++) {
            p.row(i) += delta_p.row(i);
            // collision detection: clip particle positions to boundary walls
	    p(i, 0) = clip(p(i, 0), -1*boundary, boundary);
            p(i, 1) = clip(p(i, 1), -1*boundary, boundary);
	    if (p(i, 2) < 0){ // TO DO: implement ceiling
	        p(i,2) = 0;
	    }
        }
    }

    for (int i=0; i<N; i++) {
        v.row(i) = (p.row(i) - x.row(i)) / dt;
    }
    
    /*for (int i=0; i<N; i++) {
	// VORTICITY
	VectorXd omega_i = VectorXd::Zero(3);
        for (int j : neighbor[i]) {
	    if (i != j){
	        VectorXd v_ij = v.row(j) - v.row(i);
		VectorXd grad = dW_spiky(p.row(i)-p.row(j), h);
		omega_i += cross_prod(v_ij, grad);
	    }
	}

	VectorXd eta = VectorXd::Zero(3);
	for (int j : neighbor[i]) {
	    eta += dW_spiky(p.row(i)-p.row(j), h) * omega_i.norm();
	}
	// added small constant in case denominator = 0 (source: https://www.cs.ubc.ca/~rbridson/fluidsimulation/fluids_notes.pdf section 5.1)
	VectorXd N_vort = eta / (eta.norm() + 10e-20);

	double epsilon_vort = 2/1000;
	VectorXd f_vort = epsilon_vort * cross_prod(N_vort, omega_i); 
	

	// VISCOSITY 
	VectorXd viscosity = VectorXd::Zero(3);
	double c = 0.01;
        for (int j : neighbor[i]) {
	    viscosity +=(v.row(j) - v.row(i)) * W_poly6(p.row(i) - p.row(j), h);
	}	
	viscosity = c * viscosity;

	v.row(i) += viscosity;
        v.row(i) += dt * f_vort;
    }*/

    return std::make_tuple(p, v);
}

PYBIND11_MODULE(pbf, m) {
    m.def("W_poly6", &W_poly6);
    m.def("dW_spiky", &dW_spiky);
    m.def("step", &step);
}
