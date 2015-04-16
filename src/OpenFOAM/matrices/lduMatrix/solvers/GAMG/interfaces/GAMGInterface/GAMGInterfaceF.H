#pragma once

namespace Foam
{
    struct GAMGInterfaceAgglomerateCoeffs
    {
        const scalar* ff;
        const label* sort;

        GAMGInterfaceAgglomerateCoeffs
        (
            const scalar* _ff,
            const label* _sort
        ):
            ff(_ff),
            sort(_sort)
        {}

        __HOST____DEVICE__
        scalar operator()(const label& start, const label& end)
        {
            scalar out = 0;

            for(label i = start; i<end; i++)
            {
                out += ff[sort[i]];
            }

            return out;
        }
        
    };
}
