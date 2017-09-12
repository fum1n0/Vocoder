#pragma once

#include"header.hpp"


namespace internale {

	/**
	@brief Serialize Matrix to 1D array
	@param matrix
	@param array1d: pointer of 1D array
	*/
	template <class Matrix, class T>
	inline void matrix_to_array1d(Matrix& matrix, T* array1d)
	{
		const int cols = matrix.size();
		const int rows = matrix[0].size();

		for (int i = 0; i < cols; ++i) {
			std::copy(matrix[i].begin(), matrix[i].end(), array1d);
			array1d += rows;
		}
	}

} // end namespace internale

  /**
  @brief Eigen::Matrix from matrix
  @param matrix
  @return Eigen::Matrix
  */
template <class ValueType = float, class Matrix>
Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic>
eigen_matrix(Matrix& matrix)
{
	const int cols = matrix.size();
	const int rows = matrix[0].size();

	// copy to 1d array
	std::unique_ptr<ValueType[]> array1d(new ValueType[rows*cols]);
	internale::matrix_to_array1d(matrix, array1d.get());

	// eigen matrix from 1darray
	return Eigen::Map<Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic> >(array1d.get(), rows, cols);
}

/**
@brief Eigen::Vector from vector
@param vector
@return Eigen::Vector
*/
template <class Vector>
Eigen::Matrix<typename Vector::value_type, Eigen::Dynamic, 1>
eigen_vector(Vector& vector)
{
	typedef typename Vector::value_type value_type;
	return Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, 1> >(&vector[0], vector.size(), 1);
}
