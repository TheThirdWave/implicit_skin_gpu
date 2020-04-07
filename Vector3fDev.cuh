#ifndef VECTOR3FDEV
#define VECTOR3FDEV


class Vector3FDev 
{
private:

public:

	float x;
	float y;
	float z;

	__host__ __device__ Vector3FDev()
	{
		x = 0;
		y = 0;
		z = 0;
	}

	__host__ __device__ Vector3FDev(float in)
	{
		x = in;
		y = in;
		z = in;
	}

	__host__ __device__ Vector3FDev(float xa, float ya, float za)
	{
		x = xa;
		y = ya;
		z = za;
	}

	__host__ __device__ Vector3FDev operator+(const Vector3FDev& a)
	{
		return Vector3FDev(a.x + x, a.y + y, a.z + z);
	}

	__host__ __device__ Vector3FDev operator-(const Vector3FDev& a)
	{
		return Vector3FDev(x - a.x, y - a.y, z - a.z);
	}

	__host__ __device__ float dot(const Vector3FDev& a)
	{
		return x * a.x + y * a.y + z * a.z;
	}
};
#endif