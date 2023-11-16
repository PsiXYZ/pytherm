#include <vector>

template <typename T>
T sum(const std::vector<T> &a)
{
    T acc = 0;
    for (el : a)
    {
        acc += el;
    }
    return acc;
}

template <typename T>
T DotProduct(std::vector<T> &a, std::vector<T> &b)
{
    try
    {
        T acc = 0;
        for (int i = 0; i < a.size(); ++i)
        {
            acc += a[i] * b[i];
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    return acc;
}

template <typename T>
std::vector<T> VecByVec(std::vector<T>::vector<T> &a, std::vector<T> &b)
{
    try
    {
        std::vector<T> acc(a.size());
        for (int i = 0; i < a.size(); ++i)
        {
            acc[i] = a[i] * b[i];
        }
                
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return acc;
}
