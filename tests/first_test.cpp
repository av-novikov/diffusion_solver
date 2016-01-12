#include "gtest/gtest.h"

#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>

#include "Scene.h"
#include "model/Gas1D/Gas1D.h"

using std::string;
using std::vector;
using std::ifstream;

/*void setDataFromFile(vector< pair<double,double> >& vec, string fileName)
{
	ifstream file;
	file.open(fileName.c_str(), ifstream::in);
	
	double temp1, temp2;
	while( !file.eof() )
	{
		file >> temp1;
		if( file.eof() )
			break;
		file >> temp2;
		vec.push_back(make_pair(temp1, temp2));
	}

	file.close();
}*/

TEST(FirstTest, First)
{
	/*gas1D::Properties props;
	setDataFromFile(props.z_factor, "props/z.txt");
	Scene<gas1D::Gas1D, gas1D::Gas1DSolver, gas1D::Properties> scene;
	scene.load(props);

	const Gas1D* model = scene->getModel();

	for(int i = 1; i < 600; i++)
		ASSERT_EQ( model->props_gas.z.Solve( (double)(i) ), 1.0 );*/
}