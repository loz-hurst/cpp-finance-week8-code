cmake_minimum_required(VERSION 3.12)
project(cpp_finance_week8_code)

set(CMAKE_CXX_STANDARD 14)

add_executable(cpp_finance_week8_code main.cpp BlackScholes.hpp ExplicitFDM.cpp ExplicitFDM.hpp ImplicitFDM.cpp ImplicitFDM.hpp)