// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
__HOST____DEVICE__
inline Foam::Polynomial<PolySize>::Polynomial()
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    for (int i = 0; i < PolySize; ++i)
    {
        this->v_[i] = 0.0;
    }
}


template<int PolySize>
__HOST____DEVICE__
inline Foam::Polynomial<PolySize>::Polynomial
(
    const Polynomial<PolySize>& poly
)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(poly),
    logActive_(poly.logActive_),
    logCoeff_(poly.logCoeff_)
{}


template<int PolySize>
__HOST____DEVICE__
inline Foam::Polynomial<PolySize>::Polynomial(const scalar coeffs[PolySize])
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    for (int i=0; i<PolySize; i++)
    {
        this->v_[i] = coeffs[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
__HOST____DEVICE__
inline bool Foam::Polynomial<PolySize>::logActive() const
{
    return logActive_;
}


template<int PolySize>
__HOST____DEVICE__
inline Foam::scalar Foam::Polynomial<PolySize>::logCoeff() const
{
    return logCoeff_;
}


template<int PolySize>
__HOST____DEVICE__
inline Foam::scalar Foam::Polynomial<PolySize>::value(const scalar x) const
{
    scalar val = this->v_[0];

    // avoid costly pow() in calculation
    scalar powX = 1;
    for (label i=1; i<PolySize; ++i)
    {
        powX *= x;
        val += this->v_[i]*powX;
    }

    if (logActive_)
    {
        val += logCoeff_*log(x);
    }

    return val;
}


template<int PolySize>
__HOST____DEVICE__
inline Foam::scalar Foam::Polynomial<PolySize>::derivative(const scalar x) const
{
    scalar deriv = 0;

    if (PolySize > 1)
    {
        // avoid costly pow() in calculation
        deriv += this->v_[1];

        scalar powX = 1;
        for (label i=2; i<PolySize; ++i)
        {
            powX *= x;
            deriv += i*this->v_[i]*powX;
        }
    }

    if (logActive_)
    {
        deriv += logCoeff_/x;
    }

    return deriv;
}

template<int PolySize>
__HOST____DEVICE__
inline Foam::scalar Foam::Polynomial<PolySize>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    // avoid costly pow() in calculation
    scalar powX1 = x1;
    scalar powX2 = x2;

    scalar integ = this->v_[0]*(powX2 - powX1);
    for (label i=1; i<PolySize; ++i)
    {
        powX1 *= x1;
        powX2 *= x2;
        integ += this->v_[i]/(i + 1)*(powX2 - powX1);
    }

    if (logActive_)
    {
        integ += logCoeff_*((x2*log(x2) - x2) - (x1*log(x1) - x1));
    }

    return integ;
}


template<int PolySize>
__HOST____DEVICE__
inline typename Foam::Polynomial<PolySize>::intPolyType
Foam::Polynomial<PolySize>::integral(const scalar intConstant) const
{
    intPolyType newCoeffs;

    newCoeffs[0] = intConstant;
    forAll(*this, i)
    {
        newCoeffs[i+1] = this->v_[i]/(i + 1);
    }

    return newCoeffs;
}


template<int PolySize>
__HOST____DEVICE__
inline typename Foam::Polynomial<PolySize>::polyType
Foam::Polynomial<PolySize>::integralMinus1(const scalar intConstant) const
{
    polyType newCoeffs;

    if (this->v_[0] > VSMALL)
    {
        newCoeffs.logActive_ = true;
        newCoeffs.logCoeff_ = this->v_[0];
    }

    newCoeffs[0] = intConstant;
    for (label i=1; i<PolySize; ++i)
    {
        newCoeffs[i] = this->v_[i]/i;
    }

    return newCoeffs;
}


