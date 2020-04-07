#pragma once
#include <../eigen-3.3.7/Eigen/Dense>
#include <vector>
#include <iostream>

using Eigen::MatrixXf;
using Eigen::VectorXf;


class SerialFunctions
{
public:
	SerialFunctions();
	~SerialFunctions();

	static VectorXf eigenColPivHouseHolderSolve(MatrixXf* coeff, VectorXf* solve);
private:
};

