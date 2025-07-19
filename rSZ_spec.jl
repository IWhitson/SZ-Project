### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 037e9950-b7c3-44ac-bd7f-7e72235b002a
begin
	using DelimitedFiles, Printf, PyPlot
	
	# Constants
	h_const = 6.62607015e-34  # Js
	k_B = 1.38064852e-23       # J/K
	mec2 = 511.0               # keV
	c_light = 2.99792458e8     # m/s
	Tcmb0 = 2.726              # K
	Dn_DI_conversion = 13.33914078 * Tcmb0^3
	
end

# ╔═╡ 3f4875c3-eabb-4731-9346-5b5c7d8a4cb2
begin
	function coth(x)
	    return 1.0 / tanh(x)
	end
	
	function polynomial(params, x)
	    y = zeros(length(x))
	    for i in 1:length(params)
	        y .+= params[i] .* x.^(i - 1)
	    end
	    return y
	end
	
	function mod_BB(nu, beta, T, A)
	    nu0 = 857e9
	    return A .* (nu ./ nu0).^(beta + 3) .* (exp.(h_const * nu0 / k_B / T) .- 1) ./ (exp.(h_const .* nu / k_B / T) .- 1)
	end
	
end

# ╔═╡ 13962d83-0c94-43c5-88dc-1b0e34356c7f
begin
	poly_pars = readdlm("YrSZ2KCMB_polyfits.txt")
	Tconv = readdlm("KCMB2MJysr.txt")[:]
	yconv = readdlm("KCMB2YSZ.txt")[:]
	bands = readdlm("new_bands.txt", header=true)[1]
	
	band_inds = bands[:, 1]
	nband = length(band_inds)
	band_freq = hcat(bands[:, 2], bands[:, 3])
	band_x = h_const .* band_freq .* 1e9 ./ k_B ./ Tcmb0
	
end

# ╔═╡ b80e8381-6ed5-4ee4-8665-fff0b9c0dc40
begin
	Te = 10.0  # keV
	y = 1e-4
	nu = 10.0 .^ range(log10(25), stop=log10(2000), length=100) .* 1e9
	x = h_const .* nu ./ k_B ./ Tcmb0
	
	Y0 = (2 * h_const / c_light^2) * (k_B * Tcmb0 / h_const)^3 .* x.^4 .* exp.(x) ./ (exp.(x) .- 1).^2 .* (x .* coth.(x / 2) .- 4) .* y * 1e26 * 1e-6
	
	rel_path = "Te_$(Te)_Yrel.npy"
	if isfile(rel_path)
	    Yrel = npyread(rel_path)
	else
	    Yrel = Y0  # Placeholder
	    Yrel .*= Dn_DI_conversion .* x.^3
	    npysave(rel_path, Yrel)
	end
	
end

# ╔═╡ 810f722f-dece-4dcc-8408-17d97e3cf930
begin
	fig, ax = subplots(2, 1, sharex=true)
	fig.subplots_adjust(top=0.95, right=0.95)
	
	ax[1].plot(nu .* 1e-9, Y0, "k", label="Te = 0")
	ax[1].plot(nu .* 1e-9, Yrel, "k--", label=@sprintf("Te = %.1f keV", Te))
	ax[2].plot(nu .* 1e-9, Yrel .- Y0, "k--")
	ax[2].axhline(y=0.0, color="k")
	ax[1].legend(loc="best")
	
	SZsig = zeros(nband)
	nrSZsig = zeros(nband)
	for i in 1:nband
	    SZsig[i] = polynomial(poly_pars[i, :], [Te])[1] * Tconv[i] * y
	    nrSZsig[i] = 1.0 / yconv[i] * Tconv[i] * y
	    ax[1].plot(band_freq[i, :], fill(SZsig[i], 2), "k")
	    ax[2].plot(band_freq[i, :], fill(SZsig[i] - nrSZsig[i], 2), "k")
	end
	
end

# ╔═╡ d49a2555-e544-4a93-a28c-4a85a707a39a
instr_sim_results = readdlm("sensitivity_calculations.txt", header=true)[1]
sens = instr_sim_results[:, 2] .* 1e-3
beam_sizes = instr_sim_results[:, 3]
ref = findfirst(b -> b == 8, band_inds)
refSNR = 50
err_ref = SZsig[ref] / refSNR
err = err_ref .* sens ./ sens[ref] ./ (beam_sizes ./ beam_sizes[ref])
nu_eff = readdlm("nu_eff.txt")[:]

ax[1].errorbar(nu_eff, SZsig, yerr=err, fmt="k.", markersize=10)
ax[2].errorbar(nu_eff, SZsig .- nrSZsig, yerr=err, fmt="k.", markersize=10)
ax[1].set_xscale("log")
ax[2].set_xlabel("Freq / GHz")
ax[1].set_ylabel("I / (MJy sr⁻¹)")
ax[2].set_ylabel("ΔI / (MJy sr⁻¹)")
savefig("Te_$(Te)_sig_errs.png")
close("all")


# ╔═╡ 5efc4fff-010b-4753-946e-997d7c06b491
Adust_857 = 0.2 * 0.52
Tdust = 18.44
beta = 1.5
dust_sig = mod_BB.(mean(band_freq, dims=2) .* 1e9, beta, Tdust, Adust_857)[:]

fig2, ax2 = subplots()
fig2.subplots_adjust(top=0.95, right=0.95)
ax2.plot(nu .* 1e-9, Y0, "k", label="Te = 0")
ax2.plot(nu .* 1e-9, Yrel, "k--", label="Te = $(Te) keV")
ax2.plot(nu .* 1e-9, mod_BB.(nu, beta, Tdust, Adust_857), "r:", label="Dust")
ax2.plot(nu .* 1e-9, Y0 .+ mod_BB.(nu, beta, Tdust, Adust_857), "r")
ax2.plot(nu .* 1e-9, Yrel .+ mod_BB.(nu, beta, Tdust, Adust_857), "r--")

for i in 1:nband
    ax2.plot(band_freq[i, :], fill(SZsig[i] + dust_sig[i], 2), "r")
end
ax2.errorbar(nu_eff, SZsig .+ dust_sig, yerr=err, fmt="r.", markersize=10)
ax2.set_xscale("log")
ax2.set_xlabel("Freq / GHz")
ax2.set_ylabel("I / (MJy sr⁻¹)")
savefig("Te_$(Te)_dust_sig_errs.png")
close("all")


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"

[compat]
DelimitedFiles = "~1.9.1"
PyPlot = "~2.11.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "e44b9d8de66de3f0d94323a2a6f8ca4344759b18"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "b19db3927f0db4151cb86d073689f2428e524576"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "44f6c1f38f77cafef9450ff93946c53bd9ca16ff"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "9816a3826b0ebf49ab4926e2b18842ad8b5c8f04"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.4"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "d2c2b8627bbada1ba00af2951946fb8ce6012c05"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.6"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"
"""

# ╔═╡ Cell order:
# ╠═037e9950-b7c3-44ac-bd7f-7e72235b002a
# ╠═3f4875c3-eabb-4731-9346-5b5c7d8a4cb2
# ╠═13962d83-0c94-43c5-88dc-1b0e34356c7f
# ╠═b80e8381-6ed5-4ee4-8665-fff0b9c0dc40
# ╠═810f722f-dece-4dcc-8408-17d97e3cf930
# ╠═d49a2555-e544-4a93-a28c-4a85a707a39a
# ╠═5efc4fff-010b-4753-946e-997d7c06b491
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
