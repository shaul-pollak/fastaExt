import Base: length, string, ==

struct Header
    value::Vector{UInt8}
    l::Vector{UInt16}
    gl::Vector{UInt16}
    flg::Vector{Bool}
end

Base.length(x::Header) = x.l[1]
value(x::Header) = x.value[1:length(x)]
Base.string(x::Header) = String(value(x))
gen(x::Header) = x.value[1:x.gl[1]]
genstring(x::Header) = String(gen(x))

==(x::Header, y::Header) = value(x) == value(y)
==(x::Header, y::AbstractString) = string(x) == y
==(y::AbstractString, x::Header) = string(x) == y

function push!(x::Header, v::UInt8, delim)
    @inbounds x.value[length(x)+1] = v
    @inbounds x.l[1] += 1
    if v == delim
        @inbounds x.flg[1] = true
    end
    if !x.flg[1]
        @inbounds x.gl[1] += 1
    end
end
function Base.empty!(x::Header)
    setindex!(x.l, 0, 1)
    setindex!(x.gl, 0, 1)
    setindex!(x.flg, false, 1)
    return nothing
end
Base.isempty(x::Header) = x.gl[1] == 0 ? true : false

Base.show(io::IO, ::MIME"text/plain", g::Header) = print(io, string(g))
