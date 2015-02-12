#ifndef coordinateRotationFunctor_H
#define coordinateRotationFunctor_H

namespace Foam
{

__HOST____DEVICE__
inline symmTensor transformPrincipal
(
    const tensor& tt,
    const vector& st
)
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),

        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );

}

template<class Rotation>
struct coordinateRotationFunctor
{
    const Rotation rot;

    coordinateRotationFunctor(const Rotation _rot):rot(_rot){}

    __HOST____DEVICE__
    tensor operator()(const tensor& t)
    {
        return rot.transformTensor(t);
    } 
};

struct transformPrincipalFunctor
{
    const tensor R;
 
    transformPrincipalFunctor(const tensor _R):R(_R){}

    __HOST____DEVICE__
    symmTensor operator()(const vector& v)
    {
        return transformPrincipal(R,v);
    } 
};

struct transformPrincipalFieldFunctor
{
    __HOST____DEVICE__
    symmTensor operator()(const tensor& R, const vector& v)
    {
        return transformPrincipal(R,v);
    } 
};

template<class Rotation>
inline 
tmp<tensorgpuField> transformTensorField(const Rotation& rot,const tensorgpuField& st)
{
    tmp<tensorgpuField> tt(new tensorgpuField(st.size()));
    tensorgpuField& t = tt();

    thrust::transform
    (
        st.begin(),
        st.end(),
        t.begin(),
        coordinateRotationFunctor<Rotation>(rot)
    );
    
    return tt;
}

template<class Rotation>
inline 
tmp<tensorgpuField> transformTensorField(const Rotation& rot, const tensorgpuField& st, const labelgpuList& cellMap)
{
    tmp<tensorgpuField> tt(new tensorgpuField(cellMap.size()));
    tensorgpuField& t = tt();

    thrust::transform
    (
        thrust::make_permutation_iterator
        (
            st.begin(),
            cellMap.begin()
        ),
        thrust::make_permutation_iterator
        (
            st.begin(),
            cellMap.end()
        ),
        t.begin(),
        coordinateRotationFunctor<Rotation>(rot)
    );

    return tt;
}

}

#endif
