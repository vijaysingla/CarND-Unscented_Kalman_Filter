#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{


	VectorXd RMSE(4);
	RMSE << 0,0,0,0;
	// Checking validity of estimation and ground truth data
	if (estimations.size() == 0 || estimations.size() != ground_truth.size())
	{
		cout << "Invalid estimation or ground data size " << endl ;
		return RMSE;
	}

	for (int i = 0 ; i < estimations.size();++i)
	{
		VectorXd error = estimations[i] - ground_truth[i];

		// vector element -wise multiplication
		error = error.array()*error.array();
		RMSE += error ;
	}

	// Mean Calculation
	RMSE = RMSE/estimations.size();

	// Root square
	RMSE = RMSE.array().sqrt();
	return RMSE;
}
