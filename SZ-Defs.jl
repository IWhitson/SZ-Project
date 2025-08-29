# For the first cell with Expected rel and non rel SZ
function load_numeric_matrix(filename)
    lines = readlines(filename)
    data = []
    for line in lines
        if isempty(line) || startswith(strip(line), "#")
            continue
        end
        push!(data, [parse(Float64, x) for x in split(strip(line), r"[,\s]+")])
    end
    return reduce(vcat, [reshape(row, 1, :) for row in data])
end

function polynomial(coeffs, Te)
    sum(coeffs .* [Te^i for i in 0:(length(coeffs)-1)])
end


# Second Cell Conpute Chi squared rel and non using Te_10 files
# Helper: parse only numeric lines
function load_numeric_matrix(filename)
    lines = readlines(filename)
    data = []
    for line in lines
        if isempty(line) || startswith(strip(line), "#")
            continue
        end
        push!(data, [parse(Float64, x) for x in split(strip(line), r"[,\s]+")])
    end
    return reduce(vcat, [reshape(row, 1, :) for row in data])
end

function expected_nrSZ(y, yconv, Tconv, nband)
    [1.0 / yconv[i] * Tconv[i] * y for i in 1:nband]
end

function expected_relSZ(y, poly_pars, Te, Tconv, nband)
    [polynomial(poly_pars[i, :], Te) * Tconv[i] * y for i in 1:nband]
end

function load_vector(filename)
    lines = readlines(filename)
    data = Float64[]
    for line in lines
        if isempty(line) || startswith(strip(line), "#")
            continue
        end
        for x in split(strip(line), r"[,\s]+")
            push!(data, parse(Float64, x))
        end
    end
    return data
end

function chi2(observed, expected, errors)
    sum(((observed .- expected) ./ errors).^2)
end
