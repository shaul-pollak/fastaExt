import Mmap: mmap
import .Threads: @threads, threadid, nthreads
import ProgressMeter: Progress, next!, finish!
import NaturalSort
import StringViews

function build_index()::Cint
    if length(ARGS) == 0
        error("Need input file path")
    else
        f = ARGS[1]
    end
    gen_delim = 0x23
    record_delim = endswith(f, r"fna|faa|nuc|fasta") ? 0x3e : 0x0a
    fsz = filesize(f)
    mm = open(f, "r") do io
        mmap(io, Vector{UInt8}, (fsz,))
    end
    println("finding where records start")
    seps = findall(mm .== record_delim)
    new_mm(mm, seps, f, gen_delim, fsz)
    return 0
end

function new_mm(mm, seps, f, gen_delim, fsz)
    rio = open("$f.idx1", "w")
    gio = open("$f.idx2", "w")
    i2 = findfirst(view(mm, 1:min(500, fsz)) .== gen_delim) - 1
    ref_gen = view(mm, 2:i2) |> String
    g1 = 1;
    nl = UInt8('\n')
    prg = Progress(length(seps); desc="writing index")
    @inbounds for i in eachindex(seps)
        i1 = seps[i]
        i2 = findfirst(view(mm, i1:min(i1 + 500, fsz)) .== gen_delim) + i1 - 2
        i3 = findfirst(view(mm, i1:min(i1 + 1000, fsz)) .== nl) + i1 - 2
        i4 = i==length(seps) ? fsz : seps[i+1]-1
        rec_gen = view(mm, i1+1:i2) |> StringViews.StringView
        if rec_gen != ref_gen
            # write to gio
            println(gio, ref_gen, '\t', g1, '\t', position(rio)-1)
            g1 = position(rio)
            ref_gen = String(rec_gen)
        end
        # write to rio
        hdr = view(mm, i1+1:i3) |> StringViews.StringView
        println(rio, hdr, '\t', i1, '\t', i4)
        next!(prg)
    end
    println(gio, ref_gen, '\t', g1, '\t', position(rio))
    finish!(prg)
    close(rio)
    close(gio)
end

function build_header_index()::Cint
    if length(ARGS) < 1
        error("Need input file path")
    else
        p = ARGS[1]
    end
    hdrs = get_headers(p)
    open("$p.idx0", "w") do io
        for (i, hdr) in enumerate(hdrs)
            println(io, hdr, '\t', i)
        end
    end
    return 0
end

function step_through_mm(mm, rpos0, gpos0, h0, h, gen_delim, record_delim, rio, gio)
    prg = Progress(size(mm, 1); dt=1, desc="Building indexes")
    nl = UInt8('\n')
    tb = UInt8('\t')
    hstop = record_delim == nl ? tb : nl
    flg = true
    for (pos, v) in enumerate(mm)
        if v == record_delim
            if pos > 1
                p0 = position(rio)
                print(rio, string(h), '\t', rpos0[1], '\t', pos - 1, '\n')
                rpos0[1] = pos
                if isempty(h0)
                    h0 = deepcopy(h)
                elseif gen(h) != gen(h0)
                    print(gio, genstring(h0), '\t', gpos0[1], '\t', p0, '\n')
                    gpos0[1] = p0
                    h0 = deepcopy(h)
                end
            end
            empty!(h)
            flg = true
        elseif flg && v != hstop
            push!(h, v, gen_delim)
        else
            flg = false
        end
        next!(prg)
    end
    print(rio, string(h), '\t', rpos0[1], '\t', size(mm, 1))
    print(gio, genstring(h), '\t', gpos0[1], '\t', position(rio))
end

hdrgen(x, delim) = split(x, delim)[1]

function fastaext()::Cint
    if length(ARGS) == 0
        # f = "/home/user/pollak/scratch/gtdb_hq_gens/output/pyrodigal/prots.faa"
        # hv = ["GCF_000011905.1#NC_002936.3_1017#987763#988491#-1#00","GCF_000011905.1#NC_002936.3_1#261#1598#1#00","GCF_015351475.1#NZ_JADKPR010000065.1_2#906#1205#-1#00","GCF_015351475.1#NZ_JADKPR010000090.1_1#3#386#-1#11"]
        error("Input path needed")
    else
        f = ARGS[1]
        hv = ARGS[2:end]
    end
    gen_delim = '#'
    record_delim = endswith(f, r"fna|faa|nuc|fasta") ? 0x3e : 0x0a
    gend = read_idx2(f)
    rio = open("$f.idx1", "r")
    io = open(f, "r")
    cg = hdrgen(hv[1],"#")
    recd = read_idx1(rio, gend[cg])
    for h in hv
        hg = hdrgen(h, gen_delim)
        if hg != cg
            recd = read_idx1(rio, gend[hg])
            cg = String(hg)
        end
        hidx = recd[h]
        read_record(io, hidx, record_delim)
    end
    close(io)
    close(rio)
    return 0
end


"""
    extgen()::Cint

This function depends on ripgrep. Make sure it is installed before using

This function takes two arguments:
(1) the path of the concatenated fasta file
(2) genome names separated by spaces
"""
function extgen()::Cint

    # stop if not enough input arguments
    if length(ARGS) < 2
        error("need at least two arguments: path to concatenated fasta file and genome to extract")
    end

    # parse input arguments
    inpath = ARGS[1]
    gens = ARGS[2:end]

    for gen in gens

        # read idx1 coords
        coords_raw = read(`rg $gen $inpath.idx2`, String)
        coords_i1 = parse.(Int64, split(coords_raw, '\t')[2:3])
        l1 = coords_i1[2] - coords_i1[1]

        # read fasta coords
        io_1 = open("$(inpath).idx1", "r")
        seek(io_1, coords_i1[1])
        raw = read(io_1, l1) |> String
        records = split(raw, "\n")[[1, end]]
        coords = [parse(Int64, split(records[i], '\t')[i+1]) for i in 1:2]
        l = coords[2] - coords[1]
        close(io_1)

        # extract records
        io = open(inpath, "r")
        seek(io, coords[1] - 1)
        raw = read(io, l) |> String
        close(io)

        # print to stdout
        println(stdout, raw)

    end

    return 0

end


function read_idx2(f::AbstractString)
    f2 = "$f.idx2"
    lns = split.(readlines(f2), '\t')
    return Dict(x[1] => parse.(Int, x[2:end]) for x in lns)
end

function read_idx1(rio, v)
    seek(rio, v[1])
    d = read(rio, v[2] - v[1]) |> StringViews.StringView
    lnsr = split(d, '\n')
    lns = split.(lnsr, '\t')
    return Dict(x[1] => parse.(Int, x[2:end]) for x in lns)
end

function read_record(io::IOStream, hidx::Vector{Int}, record_delim)
    seek(io, hidx[1] - 1)
    if record_delim != 0x0a && peek(io) != record_delim
        seek(io, hidx[1] - 2)
    end
    d = read(io, diff(hidx)[1] + 1)
    o = StringViews.StringView(d) |> chomp
    println(stdout, o[1] == '\n' ? o[2:end] : o)
end

function get_headers(p::String)::Vector{String}
    tn = tempname()
    println("copying file to $tn")
    cp(p, tn)
    io = open(tn, "r")
    d = mmap(io)
    chunks = create_chunk_indexes(length(d), nthreads() * 4)
    fix_indexes!(d, chunks)
    vs = [String[] for _ in eachindex(chunks)]
    gds = Vector{Dict}(undef, length(vs))
    println("extracting headers")
    @threads for i in 1:length(chunks)
        v = String[]
        get_chunk_headers!(v, d, chunks[i])
        sort!(v, lt=NaturalSort.natural)
        vs[i] = v
        gds[i] = make_gendict(v)
    end
    close(io)
    rm(tn)
    println("creating output vector")
    dksv = [collect(keys(x)) for x in gds]
    dks = sort(reduce(vcat, dksv) |> unique, lt=NaturalSort.natural)
    o = Vector{String}(undef, length.(vs) |> sum)
    ptr = 0
    for k in dks
        Is = [i for i in eachindex(gds) if k in keys(gds[i])]
        v = String[]
        for i in Is
            vi = gds[i][k]
            append!(v, vs[i][vi])
        end
        (length(Is) > 1) && sort!(v, lt=NaturalSort.natural)
        l = length(v)
        ii = ptr+1:ptr+l
        o[ii] = v
        ptr = ptr + l
    end
    return o
end

function create_chunk_indexes(len, nchunks)
    chunk_len = ceil(len / nchunks) |> Int
    [(i-1)*chunk_len+1:(i * chunk_len < len ? i * chunk_len : len) for i in 1:nchunks]
end

function find_newline(x::Vector{UInt8}, inds)
    @inbounds for i in inds
        x[i] == 0x0a && return i
    end
    return inds[1]
end

function fix_indexes!(d, chunks)
    @inbounds for i in 2:length(chunks)
        i1 = chunks[i][1]
        i2 = max(1, i1 - 50_000)
        i3 = find_newline(d, i1:-1:i2)
        chunks[i-1] = chunks[i-1].start:i3
        chunks[i] = i3+1:chunks[i].stop
    end
end

function get_chunk_headers!(hdrs::Vector{String}, d::Vector{UInt8}, i)
    rs = i[d[i].==0x3e]
    @inbounds for i1 in rs
        i2 = find_newline(d, i1:min(i1 + 500, i.stop))
        Base.push!(hdrs, String(d[i1+1:i2-1]))
    end
end

function make_gendict(v)
    gens = hdrgen.(v, "#")
    o = Dict()
    i1 = 1
    cv = gens[1]
    for (i, x) in enumerate(gens)
        if x != cv
            o[cv] = i1:i-1
            i1 = i
            cv = x
        end
    end
    o[gens[end]] = i1:length(v)
    return o
end

