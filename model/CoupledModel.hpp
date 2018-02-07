#ifndef COUPLEDMODEL_HPP_
#define COUPLEDMODEL_HPP_

#include <vector>
#include <memory>

template <class modelType>
class CoupledModel
{
protected:
	std::shared_ptr<FirstModel> fmodel;
	std::vector<std::shared_ptr<SecondModel>> smodels;
public:

	CoupledModel() {};
	virtual ~CoupledModel() {};
};

#endif /* COUPLEDMODEL_HPP_ */