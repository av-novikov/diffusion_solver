#ifndef BASE_TEST_H_
#define BASE_TEST_H_

#define PRES_REL_TOL 1.E-3
#define SAT_REL_TOL 1.E-2
#define TEMP_REL_TOL 1.E-2

template <typename propsType, typename sceneType>
class BaseTest
{
protected:
	virtual propsType* getProps() = 0;

public:
	BaseTest();
	virtual ~BaseTest();

	virtual void run();
	virtual void test() = 0;

};

#endif /* BASE_TEST_H_ */