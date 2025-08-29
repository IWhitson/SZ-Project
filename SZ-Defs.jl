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
# Helper to load expected and error values from a two-column text file
function load_expected_and_error(filename)
	expected = Float64[]
	errors = Float64[]
	for line in eachline(filename)
	    if isempty(line) || startswith(strip(line), "#")
	        continue
	    end
	    vals = split(strip(line))
	    push!(expected, parse(Float64, vals[1]))
	    push!(errors, parse(Float64, vals[2]))
	end
	return expected, errors
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

# 3rd Cell, chi squared mimization w/w/o Noise
# Function to generate expected SZ signal for a given y (non-relativistic)
function expected_nrSZ(y)
    [1.0 / yconv[i] * Tconv[i] * y for i in 1:length(band_inds)]
end

# Function to generate expected SZ signal for a given y and Te (relativistic)
function expected_relSZ_2d(y, Te)
    [polynomial(poly_pars[i, :], Te) * Tconv[i] * y for i in 1:length(band_inds)]
end

# Function to add Gaussian noise to model data
function add_noise_to_data(model_data, errorbars)
    return model_data .+ randn(length(model_data)) .* errorbars
end

# --- Minimize WITHOUT noise ---
function chi2_model(y, expected_func, observed, errors)
    model = expected_func(y)
    sum(((model .- observed) ./ errors).^2)
end

function chi2_rel_2d(params, observed, errors)
    y, Te = params
    model = expected_relSZ_2d(y, Te)
    sum(((model .- observed) ./ errors).^2)
end
	
