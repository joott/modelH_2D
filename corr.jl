cd(@__DIR__)

using Plots
using DelimitedFiles
using LaTeXStrings
using Measures

function scatter_style(xl,yl)
    scatter!(
        	ylabel=yl, xlabel=xl,
            grid = :off,
            box = :on,
            foreground_color_legend = nothing,
            fontfamily = "serif-roman",
            font="CMU Serif",
            xtickfontsize = 12,
            ytickfontsize = 12,
            xguidefontsize = 12,
            yguidefontsize = 12,
            thickness_scaling=1,
            legendfontsize=10,
            yguidefontrotation=0,
            #legend_font_pointsize=14,
            #legendtitlefontsize=14,
            markersize=1,
            legend=:topright,
            margin=5mm
        )
end

function autocor_loc_2(x, beg, max, n=2)
	C = zeros(Complex{Float64},max+1)
	N = zeros(Int64,max+1)
	Threads.@threads for tau in 0:max
		for i in beg:length(x)-max
			j = i + tau
			@inbounds C[tau+1] = C[tau+1] +  (x[i]*conj(x[j]))^n
			@inbounds N[tau+1] = N[tau+1] + 1
		end
	end
	(collect(0:max),  C ./ N)
end

p = plot()

z=3.7

for N in (8 , 12,16) 
  df=readdlm("output_"*string(N)*".dat",' ')
  c = Any[] 
  (t,tmp) = autocor_loc_2((df[:,5].+df[:,6].*1.0im), 1, N^3  , 1)
  
  plot!(p, 10*t/N^z, real.(tmp)/real.(tmp)[1],xlim=(0,0.3), label=L"L=%$N")  
end 

scatter_style(L"t/L^{%$z}",L"C")

savefig("test.pdf")

#df_16=readdlm("output_16.dat",' ')
# df_16=readdlm("output_12.dat",' ')
#df_8=readdlm("output_8.dat",' ')


#
#
# (t_8,tmp) = autocor_loc_2((df_8[:,5].+df_8[:,6].*1.0im), 1, 8^2 , 1)
#
# t_16.=t_16.-1.0
# t_8.=t_8.-1.0
#
# for i in 1:5
# 	local tmp 
# 	(tmp_t,tmp) = autocor_loc_2(df_16[:,3+2*i].+df_16[:,4+2*i].*1.0im, 1, 16^3  , 1)
# 	push!(c_16,real.(tmp))
#
#   (tmp_t,tmp) = autocor_loc_2(df_8[:,3+2*i].+df_8[:,4+2*i].*1.0im, 1, 8^2  , 1)
# 	push!(c_8,real.(tmp))
# end
#
#
# z = 3.6
#
# scatter(10*t_16/12.0^z,c_16[1]/c_16[1][1],label=L"L=16, C(t,k=2\pi/L)")
# scatter!(10*t_8/8.0^z,c_8[1]/c_8[1][1],label=L"L=8, C(t,k=2\pi/L)")
#
#
#
# scatter!(8*10*t_16/16.0^z,c_16[2]/c_16[2][1],label=L"L=16, C(t,k=3\pi/L)")
# scatter!(8*10*t_8/8.0^z,c_8[2]/c_8[2][1],label=L"L=8, C(t,k=3\pi/L)")
#
# scatter!(xlim=(0,0.3),ylabel = L"C(t)",xlabel=L"t/L^{%$z}")
#
# savefig("c36.pdf")
#
#
#
# z = 4
#
# scatter(10*t_16/16.0^z,c_16[1]/c_16[1][1],label=L"L=16, C(t,k=2\pi/L)")
# scatter!(10*t_8/8.0^z,c_8[1]/c_8[1][1],label=L"L=8, C(t,k=2\pi/L)")
#
#
#
# scatter!(8*10*t_16/16.0^z,c_16[2]/c_16[2][1],label=L"L=16, C(t,k=3\pi/L)")
# scatter!(8*10*t_8/8.0^z,c_8[2]/c_8[2][1],label=L"L=8, C(t,k=3\pi/L)")
#
# scatter!(xlim=(0,0.3),ylabel = L"C(t)",xlabel=L"t/L^{%$z}")
#
# savefig("c4.pdf")
