#pragma once

namespace Foam
{

struct luGAMGNonNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x>=0;
    }
};

struct luGAMGNegative
{
    __HOST____DEVICE__
    bool operator()(const label& x)
    {
        return x<0;
    }
};

struct faceToDiagFunctor 
   : public std::unary_function<label,label>
{
    __HOST____DEVICE__
    label operator()(const label& n)
    {
        return -1 - n;
    }
};


struct GAMGSolverAgglomerateAsymFunctor
{
    const scalar* uf;
    const scalar* lf;
    const bool* flip;
    const label* sort;

    GAMGSolverAgglomerateAsymFunctor
    (
        const scalar* _uf,
        const scalar* _lf,
        const bool* _flip,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        flip(_flip),
        sort(_sort)
    {}
    
    template<class Tuple,class Tuple2>
    __HOST____DEVICE__
    thrust::tuple<scalar,scalar> operator()(const Tuple& in,const Tuple2& t)
    {
        scalar uc = thrust::get<0>(in);
        scalar lc = thrust::get<1>(in); 
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            bool flipVal = flip[index];
            scalar ufVal = uf[index];
            scalar lfVal = lf[index];
            if( ! flipVal)
            {
                uc += ufVal;
                lc += lfVal;
            }
            else
            {
                uc += lfVal;
                lc += ufVal;
            }
        }

        return thrust::make_tuple(uc,lc);
    }

};


struct GAMGSolverAgglomerateDiagAsymFunctor
{
    const scalar* uf;
    const scalar* lf;
    const label* sort;

    GAMGSolverAgglomerateDiagAsymFunctor
    (
        const scalar* _uf,
        const scalar* _lf,
        const label* _sort
    ):
        uf(_uf),
        lf(_lf),
        sort(_sort)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            label index = sort[i];
            out += uf[index] + lf[index];
        }

        return out;
    }
};


struct GAMGSolverAgglomerateSymFunctor
{
    const scalar* ff;
    const label* sort;

    GAMGSolverAgglomerateSymFunctor
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}
    
    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += ff[sort[i]];
        }

        return out;
    }

};

struct GAMGSolverAgglomerateDiagSymFunctor
{
    const scalar* ff;
    const label* sort;

    GAMGSolverAgglomerateDiagSymFunctor
    (
        const scalar* _ff,
        const label* _sort
    ):
        ff(_ff),
        sort(_sort)
    {}

    template<class Tuple>
    __HOST____DEVICE__
    scalar operator()(const scalar& s,const Tuple& t)
    {
        scalar out = s;
        label start = thrust::get<0>(t);
        label end = thrust::get<1>(t);

        for(label i = start; i<end; i++)
        {
            out += 2*ff[sort[i]];
        }

        return out;
    }
};


}

