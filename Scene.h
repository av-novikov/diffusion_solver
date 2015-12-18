#ifndef SCENE_H_
#define SCENE_H_

#include <new>
#include <string>

template <class modelType, class methodType, typename propsType>
class Scene
{
protected:
	modelType* model;
	methodType* method;

public:
	Scene();
	~Scene();
	
	void load(propsType& props);
	void setSnapshotterType(std::string type);

	void start();
};

#endif /* SCENE_H_ */