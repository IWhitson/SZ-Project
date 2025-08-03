### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ e06835d0-7017-11f0-0d27-5de33db62d77
begin
    using Pkg
    Pkg.add("Downloads")
    Pkg.add("DelimitedFiles")
end


# ╔═╡ a1f82f16-87a4-4318-ac60-23afa2b31dc9
begin
	using DelimitedFiles
	
	# Example: read conversion factors
	kcmb2ysz = readdlm("SZ-Project")
	
end

# ╔═╡ Cell order:
# ╠═e06835d0-7017-11f0-0d27-5de33db62d77
# ╠═a1f82f16-87a4-4318-ac60-23afa2b31dc9
