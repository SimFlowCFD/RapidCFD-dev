#ifndef lduMatrixSolverFunctors_H
#define lduMatrixSolverFunctors_H

namespace Foam
{

struct wAPlusBetaPAFunctor
{
    const scalar beta;

    wAPlusBetaPAFunctor(scalar _beta): beta(_beta) {}

    __HOST____DEVICE__
    scalar operator()(const scalar& wA, const scalar& pA)
    {
            return wA + beta*pA;
    }
};

struct psiPlusAlphaPAFunctor
{
    const scalar alpha;

    psiPlusAlphaPAFunctor(scalar _alpha): alpha(_alpha) {}

    __HOST____DEVICE__
    scalar operator()(const scalar& psi, const scalar& pA)
    {
            return psi + alpha*pA;
    }
};

struct rAMinusAlphaWAFunctor
{
    const scalar alpha;

    rAMinusAlphaWAFunctor(scalar _alpha): alpha(_alpha) {}

    __HOST____DEVICE__
    scalar operator()(const scalar& rA, const scalar& wA)
    {
            return rA - alpha*wA;
    }
};

}

#endif
