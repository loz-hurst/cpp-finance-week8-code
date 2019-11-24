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

#include <iostream>
#include <fstream>
#include <string>
#include "BlackScholes.hpp"
#include "ExplicitFDM.hpp"

int main() {
    const std::string filename {"../explicit.csv"};
    const double s_max {1}; // Upper bound on spot price

    // rate, sigma, maturity, strike, type
    BlackScholes::Data data {0.05, 0.2, 1, 0.5, BlackScholes::OptionType::EurCall};

    std::ofstream explicit_result {filename};
    if (explicit_result.good()) {
        ExplicitFdm::Calculate(data, s_max, explicit_result);
        explicit_result.close();
    } else {
        std::cerr << "Unable to open " << filename << std::endl;
    }
    std::cout << "Hello, World!" << std::endl;
    return 0;
}