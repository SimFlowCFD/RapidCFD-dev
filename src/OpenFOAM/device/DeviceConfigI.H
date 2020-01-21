#pragma once

namespace Foam
{
    template<>
    inline bool hasAtomicAdd<int>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<long>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<float>()
    {
        return true;
    }

    template<>
    inline bool hasAtomicAdd<double>()
    {
        return currentComputeCapability() >= 60;
    }


    template<>
    struct Foam::is_number<int>: std::true_type {};


    template<>
    struct Foam::is_number<long>: std::true_type {};


    template<>
    struct Foam::is_number<float>: std::true_type {};


    template<>
    struct Foam::is_number<double>: std::true_type {};


    template<class T>
    inline bool hasAtomicAdd()
    {
        static_assert(is_number<T>::value, "Number is required");
        return false;
    }
}

