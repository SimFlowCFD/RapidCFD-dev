
#ifndef gpuIOField_H
#define gpuIOField_H

#include "regIOobject.H"
#include "gpuField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gpuIOField
:
    public regIOobject,
    public gpuField<Type>
{

public:

    TypeName("gpuField");


    // Constructors

        //- Construct from IOobject
        gpuIOField(const IOobject&);

        //- Construct from IOobject and size (does not set values)
        gpuIOField(const IOobject&, const label size);

        //- Construct from components
        gpuIOField(const IOobject&, const gpuField<Type>&);

        //- Construct by transferring the Field contents
        gpuIOField(const IOobject&, const Xfer<gpuField<Type> >&);


    //- Destructor
    virtual ~gpuIOField();


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const gpuIOField<Type>&);

        void operator=(const gpuField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gpuIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
