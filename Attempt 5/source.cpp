#include <iostream>
#include <vector>
#include <boost\math\special_functions\bessel.hpp>
#include <matplot/matplot.h>
#include <set>
#include <fstream>
#include<string>

//#include "matplotlibcpp.h"
#pragma once;
//namespace plt = matplotlibcpp;
using std::cin;
using std::cout;
using std::endl;
using namespace boost::math;
using std::vector;
using std::string;

double summation_value(double a_value, double T_value) {
	
	
	double summation{ 0 };
	double individual_index{ 0 };
	double zero{ 0 };
	double J_o{ 0 };
	double J_1{ 0 };
	for (unsigned n{ 1 }; n < 1000; n++) {
		 zero = cyl_bessel_j_zero(0.0, n);
		 J_o = cyl_bessel_j(0.0, zero * a_value);
		 J_1 = cyl_bessel_j(1.0, zero);
		individual_index = (J_o / J_1) * (1/pow(zero, 3)) * exp(-1 * pow(zero, 2) * T_value);
		summation = summation + individual_index;

	}
	summation = 8 * summation;
	return summation;
}



int main() {
	string filename = "data.txt";
	
	vector<double> a_value(200);
	double stepwidth = 1.0 / 200.0;
	vector<double> function_value(200);
	vector<double> T{0.05, 0.1, 0.2, 0.5, 1, 1000000};



	std::ofstream outputFile(filename);
	if (outputFile.is_open()) {
		for (unsigned T_steps=0; T_steps < T.size(); T_steps++) {
			outputFile << "Run for timestep: " << T[T_steps] << endl;
			for (unsigned steps{ 0 }; steps < 200; steps++) {
				a_value[steps] = stepwidth * steps;
				function_value[steps] = (1 - pow(a_value[steps], 2) - summation_value(a_value[steps], T[T_steps]));
				
			}
			outputFile << "a value: " << endl;
			for (unsigned j{ 0 }; j < 200; j++) {
				outputFile<< a_value[j] << endl;
			}
			outputFile <<"Function value: " << endl;
			for (unsigned j{ 0 }; j < 200; j++) {
				outputFile << function_value[j] << endl;
			}
			outputFile << endl << endl;
		}
	}
	else {
		cout << "Error opening file." << endl;
	}
	
	
	



	return 0;


		

 
		






	


	


	return 0;
}