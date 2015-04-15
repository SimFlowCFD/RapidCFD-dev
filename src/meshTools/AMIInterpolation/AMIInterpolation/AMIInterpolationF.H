#pragma once

namespace Foam
{
    template<class Type, class CombineOp>
    struct AMIInterpolationInterpolateFunctor
    {
        const scalar lowWeightCorrection;
        const CombineOp cop;
        const Type* defaultValues;
        const Type* fld;
        const label* address;
        const label* addressStart;
        const scalar* weights;
        const scalar* weightsSum;
        const Type zero;

        AMIInterpolationInterpolateFunctor
        (
            const scalar _lowWeightCorrection,
            const CombineOp _cop,
            const Type* _defaultValues,
            const Type* _fld,
            const label* _address,
            const label* _addressStart,
            const scalar* _weights,
            const scalar* _weightsSum
        ):
            lowWeightCorrection(_lowWeightCorrection),
            cop(_cop),
            defaultValues(_defaultValues),
            fld(_fld),
            address(_address),
            addressStart(_addressStart),
            weights(_weights),
            weightsSum(_weightsSum),
            zero(pTraits<Type>::zero)
        {}

        __HOST____DEVICE__
        Type operator()(const label& faceI)
        {
            if (weightsSum[faceI] < lowWeightCorrection)
            {
                return defaultValues[faceI];
            }
            else
            {
                Type out = zero;

                label start = addressStart[faceI];
                label end = addressStart[faceI+1];

                for(label i = start; i< end; i++)
                {
                    cop(out, faceI, fld[address[i]], weights[i]);
                }

                return out;
            }
        }
    };
}
