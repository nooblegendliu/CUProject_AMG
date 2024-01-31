#ifndef cuAMG_UTILS_SORT_CUH_
#define cuAMG_UTILS_SORT_CUH_

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

template<typename T, typename U>
void vec_sort_host(thrust::host_vector<T>& vec1, 
	thrust::host_vector<U>& vec2,
	int start = 0, int end = -1)
{
	int n = vec1.size();
	if (end < 0) end = n;
	thrust::sort_by_key(vec1.begin() + start, vec1.begin() + end, vec2.begin() + start);
}

template <typename T, typename U>
void matrix_sort_host(thrust::host_vector<T>& row_ptr, 
	thrust::host_vector<T>& col_ind,
	thrust::host_vector<U>& vals)
{
	int n = row_ptr.size() - 1;
	for (int i = 0; i < n; ++i)
	{
		int start = row_ptr[i];
		int end = row_ptr[i + 1];
		thrust::sort_by_key(col_ind.begin() + start,
			col_ind.begin() + end,
			vals.begin() + start);
	}
}


template <typename T, typename U>
void matrix_sort_device(thrust::device_vector<T>& vec1,
	thrust::device_vector<T>& vec2,
	thrust::device_vector<U>& vec3)
{
	for (size_t i = 0; i < vec1.size() - 1; ++i) {
		size_t begin = vec1[i];
		size_t end = vec1[i + 1];
		thrust::sort_by_key(vec2.begin() + begin,
			vec2.begin() + end, vec3.begin() + begin);
	}
}



#include <numeric>
#include <vector>
#include <algorithm>

template <typename T, typename U>
void vec_sort(std::vector<T>& vec1, std::vector<U>& vec2, int start = 0, int end = -1)
{
	vec1.shrink_to_fit();
	vec2.shrink_to_fit();

	int k, prev_k;
	int n = vec1.size();
	if (end < 0) end = n;
	int size = end - start;

	std::vector<int> p(size);
	std::vector<bool> done(size, false);

	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](const int i, const int j)
		{
			return vec1[i + start] < vec1[j + start];
		});
	for (int i = 0; i < size; i++)
	{
		if (done[i]) continue;
		done[i] = true;
		prev_k = i;
		k = p[i];
		while (i != k)
		{
			std::swap(vec1[prev_k + start], vec1[k + start]);
			std::swap(vec2[prev_k + start], vec2[k + start]);
			done[k] = true;
			prev_k = k;
			k = p[k];
		}
	}
}


#endif