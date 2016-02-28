#ifndef ITERATORS_H_
#define ITERATORS_H_

#include <initializer_list>
#include <iostream>
#include <cassert>

struct Vec3Int
{
	union {
		int ind[3];
		struct {
			int r;
			int phi;
			int z;
		};
	};

	Vec3Int()
	{
		r = phi = z = 0;
	};

	Vec3Int(std::initializer_list<int> list)
	{
		int i = 0;
		for (auto value : list)
			this->ind[i++] = value;
	};
	
	Vec3Int(const Vec3Int& orig) = default;

	Vec3Int& operator=(const Vec3Int rhs)
	{
		int i = 0;
		for (auto value : rhs.ind)
			this->ind[i++] = value;

		return (*this);
	};

};

template <class cellType>
struct Iterator : public Vec3Int
{
private:
	cellType* ptr;
protected:
	Vec3Int left;
	Vec3Int right;
	Vec3Int sizes;
public:
	Iterator() : Vec3Int()
	{
		ptr = nullptr;
	};

	Iterator(cellType* _ptr, Vec3Int _left, Vec3Int _right, Vec3Int _sizes) :
		ptr(_ptr), left(_left), right(_right), sizes(_sizes)
	{
		if (ptr != nullptr)
		{
			this->phi = int(ptr->num / (sizes.r * sizes.z));
			//assert(this->phi <= right.phi && this->phi >= left.phi);
			this->r = int((ptr->num - this->phi * (sizes.r * sizes.z)) / sizes.z);
			//assert(this->r <= right.r && this->r >= left.r);
			this->z = (ptr->num % (sizes.r * sizes.z)) % sizes.z;
			//assert(this->z <= right.z && this->z >= left.z);
		}
	};

	Iterator(const Iterator& it) = default;
	
	~Iterator() {};

	Iterator& operator++()
	{
		this->z = (this->z + 1 - left.z) % (right.z - left.z + 1) + left.z;
		if (this->z == left.z)
		{
			this->r = (this->r + 1 - left.r) % (right.r - left.r + 1) + left.r;
			if (this->r == left.r)
			{
				this->phi = (this->phi + 1 - left.phi) % (right.phi - left.phi + 1) + left.phi;
				if (this->phi == left.phi)
				{
					ptr = nullptr;
					return (*this);
				}
			}
		}

		ptr += (getIdx() - ptr->num);

		return (*this);
	};

	Iterator& operator=(const Iterator& rhs)
	{
		ptr = rhs.ptr;

		Vec3Int::operator=(rhs);
		left = rhs.left;
		right = rhs.right;
		sizes = rhs.sizes;

		return (*this);
	};

	bool operator==(const Iterator& rhs)
	{
		return ptr == rhs.ptr;
	};

	bool operator!=(const Iterator& rhs)
	{
		return ptr != rhs.ptr;
	};

	cellType* operator->() const
	{
		return ptr;
	};

	inline int getIdx()
	{
		return phi * sizes.r * sizes.z + r * sizes.z + z;
	};
};

template <class cellType>
inline std::ostream& operator<<(std::ostream& os, const Iterator<cellType>& it)
{
	os << it.r << "\t" << it.phi << "\t" << it.z << std::endl;

	return os;
};

#endif /* ITERATORS_H_ */