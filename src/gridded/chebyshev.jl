export Chebyshev

struct Chebyshev{BC<:Union{Throw{OnGrid}},N} <: DegreeBC{N}
    bc::BC
end
Chebyshev(N) = Chebyshev(Throw(OnGrid()),Val(N))
Chebyshev(bc::BC,v::VN) where {BC,VN} = Chebyshev{BC,VN}(bc)
positions(deg::Chebyshev, knotvec, x) = @inbounds (axes(ax,1)[1], (x-first(ax))/(last(ax)-first(ax)))
value_weights(T::Val{N},z::F) where {N,F} = NTuple{N,F}(ChebyshevPolynomial(N,z))

mutable struct ChebyshevPolynomial{T}
    n::Int
    z::T
    _v::T
    __v::T
end
ChebyshevPolynomial(n,z::T) where T = ChebyshevPolynomial(n,z,convert(T,0.5),zero(T))
function Base.iterate(cheb::ChebyshevPolynomial{T}) where T
    return one(T), 2
end
function Base.iterate(cheb::ChebyshevPolynomial{T}, degree) where T
    degree > cheb.n && return nothing
    v = 2*cheb.z*cheb._v-cheb.__v
    cheb.__v = cheb._v
    cheb._v = v
    return v, degree+1
end
Base.length(cheb::ChebyshevPolynomial)=cheb.n
