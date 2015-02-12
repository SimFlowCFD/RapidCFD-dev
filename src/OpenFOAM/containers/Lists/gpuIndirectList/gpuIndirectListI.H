/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::gpuIndirectList<T>::gpuIndirectList
(
    const gpuList<T>& completeList,
    const labelgpuList& addr
)
:
    completeList_(const_cast<gpuList<T>&>(completeList)),
    addressing_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::gpuIndirectList<T>::size() const
{
    return addressing_.size();
}


template<class T>
inline bool Foam::gpuIndirectList<T>::empty() const
{
    return addressing_.empty();
}

template<class T>
inline const Foam::gpuList<T>& Foam::gpuIndirectList<T>::completeList() const
{
    return completeList_;
}


template<class T>
inline const Foam::gpuList<Foam::label>& Foam::gpuIndirectList<T>::addressing() const
{
    return addressing_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const Foam::gpuIndirectList<T>::iterator Foam::gpuIndirectList<T>::begin()
{
	 return thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin());
}

template<class T>
inline const Foam::gpuIndirectList<T>::iterator Foam::gpuIndirectList<T>::begin() const
{
	 return thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin());
}

template<class T>
inline const Foam::gpuIndirectList<T>::iterator Foam::gpuIndirectList<T>::end()
{
	 return thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin());
}

template<class T>
inline const Foam::gpuIndirectList<T>::iterator Foam::gpuIndirectList<T>::end() const
{
	 return thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin());
}

template<class T>
inline void Foam::gpuIndirectList<T>::operator=(const gpuList<T>& ae)
{
    if (addressing_.size() != ae.size())
    {
        FatalErrorIn("gpuIndirectList<T>::operator=(const gpuList<T>&)")
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addressing_.size() << " " << ae.size()
            << abort(FatalError);
    }

    thrust::copy(ae.begin(),ae.end(),
    		thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin()));
}


template<class T>
inline void Foam::gpuIndirectList<T>::operator=(const gpuIndirectList<T>& ae)
{
    if (addressing_.size() != ae.size())
    {
        FatalErrorIn("gpuIndirectList<T>::operator=(const gpuIndirectList<T>&)")
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addressing_.size() << " " << ae.size()
            << abort(FatalError);
    }

    thrust::copy(thrust::make_permutation_iterator(ae.completeList_.begin(),ae.addressing_.begin()),
    		     thrust::make_permutation_iterator(ae.completeList_.begin(),ae.addressing_.end()),
        		 thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin()));
}


template<class T>
inline void Foam::gpuIndirectList<T>::operator=(const T& t)
{
	thrust::fill(thrust::make_permutation_iterator(completeList_.begin(),addressing_.begin()),
			     thrust::make_permutation_iterator(completeList_.begin(),addressing_.end()),
			     t);
}


// ************************************************************************* //
