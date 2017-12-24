#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	VectorXd resid;
	for(int i=0; i < estimations.size(); ++i){
        
        resid = (estimations[i] - ground_truth[i]);
        resid = resid.array()*resid.array();
		rmse += resid;
	}

	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}