#pragma once
#include <omp.h>
#include <thread>
#include <vector>
#include <cfloat>

namespace utils
{
	extern unsigned int num_threads;

	template <typename Iterator, typename Function>
	void For(Iterator from, Iterator to, Iterator incrementBy, bool inParallel, Function function)
	{
		if (inParallel)
		{
#pragma omp parallel for
			for (Iterator i = from; i < to; i += incrementBy)
			{
				function(i);
			}
		}
		else
		{
			// Loop for sequential computations.
			// '#pragma omp parallel for if(inParallel)' in case of inParallel==false will execute
			// such loop in one thread, but still with overhead of OpenMP thread creation.
			for (Iterator i = from; i < to; i += incrementBy)
			{
				function(i);
			}
		}
	}

	template <typename T>
	void DeleteJaggedArray(T** jaggedArray, size_t rows, size_t columns)
	{
		for (size_t i = 0; i < rows; i++)
		{
			delete[] jaggedArray[i];
			jaggedArray[i] = nullptr;
		}
		delete[] jaggedArray;
		jaggedArray = nullptr;
	}

	template <typename T>
	T** CreateJaggedArray(size_t rows, size_t columns)
	{
		auto res = new T*[rows];
		for (size_t i = 0; i < columns; i++)
		{
			res[i] = new T[columns];
		}

		return res;
	}

	void SolveCsabaDeboorTridiagonalSystem(
			double main_diagonal_value, double right_side[], unsigned int
	num_equations, double last_main_diagonal_value = DBL_MIN);

	void SolveCsabaDeboorTridiagonalSystemBuffered(double main_diagonal_value,
												   double right_side[], unsigned int num_equations, double L[],
												   double U[], double Y[], double D[], double
												   last_main_diagonal_value= DBL_MIN);

	void SolveDeboorTridiagonalSystemBuffered(double lower_diagonal_value,
		double main_diagonal_value, double upper_diagonal_value,
		double right_side[], size_t num_equations, double buffer[],
		double last_main_diagonal_value = DBL_MIN);

	void solveDeboorTridiagonalSystemBuffered(double main_diagonal_value,
											  double* right_side, size_t num_equations,
											  double* buffer, double
											  last_main_diagonal_value = DBL_MIN);
	
};

