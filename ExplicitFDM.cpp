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

#include "ExplicitFDM.hpp"
#include "BlackScholes.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ExplicitFdm {

    void Calculate(const BlackScholes::Data & data, const double s_max, const int s_steps, const int t_steps, std::ostream & out_str) {

        const double delta_s {s_max/s_steps};
        const double delta_t {data.maturity/t_steps};

        // Check the file's good
        if (! out_str.good() ) {
            std::cerr << "Output stream isn't good for ExplicitFDM - aborting" << std::endl;
            return;
        }

        // We're going to move column-wise through the space, from time=T to time=0 (i.e. backwards).

        /* We need the values from the previous time step to calculate this one - use the boolean flipping trick to
         * avoid an expensive copy at each iteration (just alternate which is the 'source' and which it the 'destination'.
         */
        std::vector<double> values_time_1 (s_steps, 0.0);
        std::vector<double> values_time_2 (s_steps, 0.0);
        bool reading_payoffs_1 {true};

        /* Our first column of data is to take the final value (at time T - which is the payoff).
         * Conveniently (or by design!) value at 0 is also 0 and is S_max-K at S_max, which happens to match the
         * boundary conditions for this initial case.
         */
        for (int i {0}; i < s_steps; ++i) {
            // For each spot price point, calculate the value at that point
            values_time_1.at(i) = std::max(i * delta_s - data.strike, 0.0); // Standard payoff function at boundary
        }

        // Start looping
        for (int i {0}; i < t_steps; ++i) {
            // For each time step (conceptually starting at T-1 and going back to 0)

            // Set boundary values (which cannot be calculated due to lack of points outside the bounds of the solution
            const double upper_bound {s_max-data.strike};
            if (reading_payoffs_1) {
                values_time_2.at(0) = 0; // 0 for call
                values_time_2.at(s_steps-1) = upper_bound;
            } else {
                values_time_1.at(0) = 0;
                values_time_1.at(s_steps-1) = upper_bound;
            }

            // Output the bottom boundary
            out_str << (data.maturity-(i*delta_t)) << ',' << 0 << ','<< ((reading_payoffs_1) ? values_time_2.at(0) : values_time_1.at(0)) << "\n";

            // For each spot price point, skipping the boundary points (0 and s_max), find values from the previous ones
            for (int j {1}; j < s_steps-1; ++j) {
                // Calculate the value at the next price point

                // Start with co-efficients

                // 0.5*sigma^2*j^2*delta_t-0.5*r*j*delta_t
                double alpha {0.5*delta_t*j*(data.sigma*data.sigma*j-data.rate)};
                // 1-sigma^2*j^2*delta_t-r*delta_t
                double beta {1-delta_t*(data.sigma*data.sigma*j*j-data.rate)};
                // 0.5*sigma^2*j^2*delta_t+0.5*r*j*delta_t
                double gamma {0.5*delta_t*j*(data.sigma*data.sigma*j+data.rate)};

                // Calculate out new value
                if (reading_payoffs_1) {
                    values_time_2.at(j) = alpha * values_time_1.at(j-1) + beta * values_time_1.at(j) + gamma * values_time_1.at(j+1);
                } else {
                    values_time_1.at(j) = alpha * values_time_2.at(j-1) + beta * values_time_2.at(j) + gamma * values_time_2.at(j+1);
                }

                // Store the result (x (time), y (spot price), val)
                out_str << (data.maturity-(i*delta_t)) << ',' << j*delta_s << ','<< ((reading_payoffs_1) ? values_time_2.at(j) : values_time_1.at(j)) << "\n";
            }
            // Output the top boundary
            out_str << (data.maturity-(i*delta_t)) << ',' << s_max << ','<< ((reading_payoffs_1) ? values_time_2.at(s_steps-1) : values_time_1.at(s_steps-1)) << "\n";

            /* Now we have gone down all value of spot price - and flip to using the other vector as 'old'
             * (avoids copying 'new' -> 'old', for efficiency)
             */
            reading_payoffs_1 = !reading_payoffs_1;
        }
    }
}
