#include "SerialFunctions.h"

SerialFunctions::SerialFunctions()
{

}

SerialFunctions::~SerialFunctions()
{

}

VectorXf  SerialFunctions::eigenColPivHouseHolderSolve(MatrixXf* coeff, VectorXf* solve)
{
	return coeff->colPivHouseholderQr().solve(*solve);
}
