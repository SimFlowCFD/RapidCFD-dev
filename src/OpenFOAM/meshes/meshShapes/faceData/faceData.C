#include "faceData.H"
#include "facegpuList.H"
#include "gpuList.C"


namespace Foam 
{

template class gpuList<faceData>;

const char* const faceData::typeName = "faceDataData";
const faceData Foam::faceData::zero(0,0);

Istream& operator>>(Istream& is, faceData& cd)
{
    is.readBegin("faceData");
    is >> cd.start_ >> cd.size_;
    is.readEnd("cellData");

    return is;
}

Ostream& operator<<(Ostream&os , const faceData& cd)
{
    os  << token::BEGIN_LIST
        << cd.start_ << token::SPACE << cd.size_
        << token::END_LIST;

    return os;
}

}

