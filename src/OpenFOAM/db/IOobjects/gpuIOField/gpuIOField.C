#include "gpuIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::gpuIOField<Type>::gpuIOField(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("gpuIOField::gpuIOField(const IOobject&)")
            << "gpuIOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but gpuIOField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class Type>
Foam::gpuIOField<Type>::gpuIOField(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("gpuIOField::gpuIOField(const IOobject&, const label)")
            << "gpuIOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but gpuIOField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        gpuField<Type>::setSize(size);
    }
}


template<class Type>
Foam::gpuIOField<Type>::gpuIOField(const IOobject& io, const gpuField<Type>& f)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("gpuIOField::gpuIOField(const IOobject&, const Field<Type>&)")
            << "gpuIOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but gpuIOField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        gpuField<Type>::operator=(f);
    }
}


template<class Type>
Foam::gpuIOField<Type>::gpuIOField(const IOobject& io, const Xfer<gpuField<Type> >& f)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "gpuIOField::gpuIOField(const IOobject&, const Xfer<Field<Type> >&)"
        )   << "gpuIOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but gpuIOField does not support automatic rereading."
            << endl;
    }

    gpuField<Type>::transfer(f());

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::gpuIOField<Type>::~gpuIOField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::gpuIOField<Type>::writeData(Ostream& os) const
{
    return (os << static_cast<const gpuField<Type>&>(*this)).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::gpuIOField<Type>::operator=(const gpuIOField<Type>& rhs)
{
    gpuField<Type>::operator=(rhs);
}


template<class Type>
void Foam::gpuIOField<Type>::operator=(const gpuField<Type>& rhs)
{
    gpuField<Type>::operator=(rhs);
}


// ************************************************************************* //
