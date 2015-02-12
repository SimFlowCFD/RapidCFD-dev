#include "Ostream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

template<class T>
Foam::gpuList<T>::gpuList(Istream& is)
:
    v_(0)
{
    operator>>(is, *this);
}

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, gpuList<T>& gL)
{
    List<T> L(is);

    gL.operator=(L);

    return is;
}


template<class T>
Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const Foam::gpuList<T>& gL)
{
    List<T> L(gL.size());

    gpu_api::copy(gL.begin(),gL.end(),L.begin());

    os << L;

    return os;
}


template<class T>
void Foam::gpuList<T>::writeEntry(Ostream& os) const
{
    if
    (
        size()
        && token::compound::isCompound
        (
            "List<" + word(pTraits<T>::typeName) + '>'
        )
    )
    {
        os  << word("List<" + word(pTraits<T>::typeName) + '>') << " ";
    }

    os << *this;
}


template<class T>
void Foam::gpuList<T>::writeEntry(const word& keyword, Ostream& os) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


// ************************************************************************* //
