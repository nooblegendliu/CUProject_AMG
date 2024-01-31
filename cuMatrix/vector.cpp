#include "vector.h"
using namespace AMG;

void Vector::set_const_value(data_t alpha)
{
    for (index_t i = 0; i < num_values; i++)
    {
        values[i] = alpha;
    }
}

void Vector::set_rand_values(data_t min, data_t max)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    for (index_t i = 0; i < num_values; i++)
    {
        values[i] = dis(gen);
    }
}

void Vector::axpy(Vector& x, data_t alpha)
{
    for (index_t i = 0; i < num_values; i++)
    {
        values[i] += x.values[i] * alpha;
    }
}

void Vector::copy(const Vector& y)
{
    num_values = y.num_values;
    values.resize(num_values);
    std::copy(y.values.begin(), y.values.end(), values.begin());
}

void Vector::scale(data_t alpha)
{
    for (index_t i = 0; i < num_values; i++)
    {
        values[i] *= alpha;
    }
}

data_t Vector::norm(index_t p)
{
    data_t result = 0.0;
    double val;
    for (index_t i = 0; i < num_values; i++)
    {
        val = values[i];
        if (fabs(val) > zero_tol)
            result += pow(val, p);
    }
    return pow(result, 1.0 / p);
}

void Vector::print(const char* vec_name)
{
    printf("Size = %d\n", num_values);
    for (int i = 0; i < num_values; i++)
    {
        if (fabs(values[i]) > zero_tol)
            printf("%s[%d] = %e\n", vec_name, i, values[i]);
    }
}

data_t& Vector::operator[](const int index)
{
    return values[index];
}


data_t Vector::inner_product(Vector& x)
{
    data_t result = 0.0;

    for (int i = 0; i < num_values; i++)
    {
        result += values[i] * x[i];
    }

    return result;
}