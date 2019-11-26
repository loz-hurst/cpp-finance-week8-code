/*
 * Code for week 8 exercises of C++ for Finance.
 *
 * Copyright 2019 Laurence Alexander Hurst
 *
 * This file is part of C++ for Finance.
 *
 *     C++ for Finance is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     C++ for Finance is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * See the file LICENCE in the original source code repository for the
 * full licence.
 */

#include "ImplicitFDM.hpp"
#include "BlackScholes.hpp"
#include <algorithm>
#include <ostream>
#include <vector>

#include <iostream>

namespace ImplicitFdm {
    void Calculate(const BlackScholes::Data &data, const double s_max, const int s_steps, const int t_steps, std::ostream &out_str) {
        // Firstly space to store our values for alpha, beta and gamma

        // Diagonal matrices are mostly 0s, so just store the diagonals for efficiency
        std::vector<double> alpha(s_steps);
        std::vector<double> beta(s_steps);
        std::vector<double> gamma(s_steps);

        // Some useful variables
        double delta_s {s_max/s_steps};
        double delta_t {data.maturity/t_steps};

        // Check the file's good
        if (! out_str.good() ) {
            std::cerr << "Output stream isn't good for ImplicitFDM - aborting" << std::endl;
            return;
        }

        // Now we need the "last" column of values for C_j - payoff initially

        /* Like in the explicit version we are usin the boolean flipping trick to avoid an expensive copy at each
         * iteration (just alternate which is the 'source' and which it the 'destination'.
         */
        std::vector<double> values_time_1 (s_steps, 0.0);
        std::vector<double> values_time_2 (s_steps, 0.0);
        bool reading_payoffs_1 {true};

        /* Our first column of data is to take the final value (at time T - which is the payoff).
         * Conveniently (or by design!) value at 0 is also 0 and is S_max-K at S_max, which happens to match the
         * boundary conditions for this initial case.
         */
        // Need this for alpha, beta and gamma - 1/(1-r*delta_t)
        const double r_delta_t_recip {1/(1-data.rate*delta_t)};

        for (int i {0}; i < s_steps; ++i) {
            // alpha = 1/(1-r*delta_t)((sigma^2j^2delta_t)/2-(r*j*delta_t)/2)
            alpha.at(i) = r_delta_t_recip*(i*delta_t*0.5*(data.sigma*data.sigma*i-data.rate));
            // beta =  1/(1-r*delta_t)(1+(sigma^2j^2delta_t)/2)
            beta.at(i) = r_delta_t_recip*(1+(data.sigma*data.sigma*i*i*delta_t)/2);
            // gamma = 1/(1-r*delta_t)((r*j*delta_t)/2-(sigma^2j^2delta_t)/2)
            gamma.at(i) = r_delta_t_recip*(i*delta_t*0.5*(data.rate-data.sigma*data.sigma*i));

            // For each spot price point, calculate the value at that point
            values_time_1.at(i) = std::max(i * delta_s - data.strike, 0.0); // Standard payoff function at boundary
        }

        // Begin looping
        for (int i {0}; i < (t_steps); ++i) {
            // Calculate out alpha, beta and gamma values

            // Set boundary values (which cannot be calculated due to lack of points outside the bounds of the solution
            const double upper_bound {s_max - data.strike};
            if (reading_payoffs_1) {
                values_time_2.at(0) = 0; // 0 for call
                values_time_2.at(s_steps-1) = upper_bound;
            } else {
                values_time_1.at(0) = 0;
                values_time_1.at(s_steps-1) = upper_bound;
            }

            // Firstly we need to calculate z

            std::vector<double> z (s_steps);
            // z_0 = Cj_0
            z.at(0) = ((reading_payoffs_1) ? values_time_1.at(0) : values_time_2.at(0));

            // We will need to store the values of d
            std::vector<double> d (s_steps);
            // d_0 = beta_0
            d.at(0) = beta.at(0);
            for (int j {1}; j < (s_steps - 1); ++j) {
                // l_i = alpha_i/d_(i-1)
                double l_i {alpha.at(j)/d.at(j-1)};

                // z_i = Cj_i - l_i*Cj_(i-1)
                z.at(j) = ((reading_payoffs_1) ? values_time_1.at(j) : values_time_2.at(j)) - l_i*
                      ((reading_payoffs_1) ? values_time_1.at(j-1) : values_time_2.at(j-1));

                // d_i = beta_i-l_i*u_(i-1)
                // u_i = gamma_i
                // d_i will be d_(i-1) next loop
                d.at(j) = beta.at(j) - l_i*gamma.at(j-1);
            }

            // Now we can calculate the values for C_(j-1)

            // upper boundary (going backwards)
            out_str << (data.maturity-(i*delta_t)) << ',' << s_max << ','<< ((reading_payoffs_1) ? values_time_2.at(s_steps-1) : values_time_1.at(s_steps-1)) << "\n";
            // Remember we have to go backwards for this
            for (int j {s_steps-2}; 0<=j; --j) {
                // C^{j-1}_i=(z_i-u_iC^{j-1}_{i+1})/d_i
                // u_i = gamma_i
                double result {(z.at(j)-gamma.at(j)*((reading_payoffs_1) ? values_time_2.at(j+1) : values_time_1.at(j+1)))/d.at(j)};
                // Store the result in the one not being read from.
                if (reading_payoffs_1) {
                    values_time_2.at(j) = result;
                } else {
                    values_time_1.at(j) = result;
                }

                // Store the result (x (time), y (spot), z (value)
                out_str << (data.maturity-(i*delta_t)) << ',' << j*delta_s << ','<< ((reading_payoffs_1) ? values_time_2.at(j) : values_time_1.at(j)) << "\n";
            }

            // End of loop - the result we've just filled should now be the 'old' one, so flip which one we're using
            reading_payoffs_1 = !reading_payoffs_1;

        }
    }
}

