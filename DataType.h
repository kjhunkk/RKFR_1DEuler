#pragma once
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <memory>

// Define condition type
// e.g. initial condition, flux, etc
typedef std::string Type;

// Define double vector
typedef std::vector<double> vector_d;

// Define integer type
typedef int int_t;

// Define real type
typedef double real_t;

// Frequently used constants
#define epsilon 1.0e-8

#define GHOST 3

#define QuadDegree 3

#define POST_GRID_NUM 10

#define BOUND_M 10

#define GAMMA 1.4

const real_t CONST13 = 1.0 / 3.0;

const real_t CONST23 = 2.0 / 3.0;

const real_t CONST43 = 4.0 / 3.0;

const real_t CONST59 = 5.0 / 9.0;

const real_t CONST89 = 8.0 / 9.0;

const real_t CONST16 = 1.0 / 6.0;

const real_t CONST1_12 = 1.0 / 12.0;

const real_t CONST1_60 = 1.0 / 60.0;

// Class alert
class Alert
{
protected:
	Alert();
	~Alert();

public:
	// Print error message
	static void error(const std::string& str, const std::string& file, int_t line);
	static void message(const std::string& str);

private:
	// Alert variable
	static Alert _alert;
};

// Message macro
#define ERROR(str) Alert::error(str, __FILE__, __LINE__)
#define MESSAGE(str) Alert::message(str);