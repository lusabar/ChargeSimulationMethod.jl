using Combinatorics


fields = Array{Any}(undef, 15)
c = 1
for i in 1:5
    for j in (i+1):6
        fields[c] = :($(Symbol("v" * string(i) * string(j)))::Int)
        global c += 1
    end
end
@eval struct AngRep
    $(fields...)
end
AngRep(str::String) = AngRep(conf_to_ang(str)...)

mutable struct Configuration
    name::String
    angrep::AngRep
    equiv::Vector{String}
    inv::Vector{String}
end

function Configuration(str::String, symb::Symbol)
    if symb == :hex
        angrep = AngRep(str)
    end
    Configuration(str, angrep, [], [])
end


phases = Dict(
    'a' => 0,
    'b' => -60,
    'c' => -120,
    'd' => -180,
    'e' => -240,
    'f' => -300,
)
phasestri = Dict(
    'a' => 0,
    'b' => -120,
    'c' => 120,
)

function fixtuple(t::NTuple{15,Int64})
    a = []
    for element in t
        if -360 <= element <= -180
            push!(a, element + 360)
        elseif 180 < element < 360
            push!(a, element - 360)
        else
            push!(a, element)
        end
    end
    return Tuple(i for i in a)
end

function rm_reversed(arr::Vector{NTuple{15,Int64}})
    seen = Set{Tuple}()
    result = []

    for el in arr
        reversed = reverse(el)
        if el ∉ seen && reversed ∉ seen
            push!(seen, el)
            push!(seen, reversed)
            push!(result, el)
        end
    end
    return result
end

function rm_equiv(arr::Vector{Configuration})
    seen = []
    result = []

    for (id, conf) in enumerate(arr)
        if conf.angrep ∉ seen
            push!(seen, conf.angrep)
            push!(result, conf)
        else
            for (i, c) in enumerate(result)
                if c.angrep == conf.angrep
                    push!(result[i].equiv, conf.name)
                end
            end
        end
    end
    return result
end

global huh

function rm_reversed(arr::Vector{Configuration})
    seen = []
    result = []

    has_reverse = []

    for confs in combinations(arr, 2) |> collect
        if !isdisjoint([confs[1].inv...], [confs[2].name confs[2].equiv...])
            push!(result, confs[1])
            append!(has_reverse, [confs[1].name confs[2].name confs[1].equiv... confs[2].equiv...])
        end
    end

    for conf in arr
        for name in [conf.name conf.equiv... conf.inv...]
            if name ∉ has_reverse
                push!(result, conf)
                break
            end
        end
    end

    global huh = has_reverse

    return result
end

# ╔═╡ 8fc1e54f-4a30-405b-958c-86349842aea5
function phase_diff(c1::Char, c2::Char)
    return phases[c1] - phases[c2]
end

# ╔═╡ 6a887cfb-0204-4388-96b1-01f8b3168767
function conf_to_ang(str::String)
    m = []
    for i in 1:5
        for j in (i+1):6
            append!(m, phase_diff(str[i], str[j]))
        end
    end

    return fixtuple(Tuple(i for i in m))
end





hex_strs = String.(collect(Combinatorics.permutations("abcdef")))
hex_confs = Array{Configuration}(undef, 720)
for (id, str) in enumerate(hex_strs)
    hex_confs[id] = Configuration(str, :hex)
end

new_arr = rm_equiv(hex_confs)
#println("The length of the new array is $(length(new_arr))")

# Names the duplicates for each config
for (id, conf) in enumerate(new_arr)
    for str in [conf.name conf.equiv...]
        rev = reverse(str)
        push!(new_arr[id].inv, rev)
    end
end

#for i in new_arr
#    println(i)
#end

newer_arr::Vector{Configuration} = new_arr

newernewer_arr = rm_reversed(newer_arr)

for i in newernewer_arr
    println(i)
end

#println("\n\nAgora vai:")
#for str in hex_strs
#    if str ∉ huh
#        println(str)
#    end
#end

# exit()
# show(conf1)
#
# exit()
# angles_rep = []
#
# for conf in config_phases_hex
#     push!(angles_rep, conf_to_ang(conf))
#
# end
#
#
# y = fixtuple.(angles_rep)
#
# z = rm_reversed(y)
#
# length(z)
#
# conf_to_ang("ecadfb")
#
# conf_to_ang("dfbeca")
