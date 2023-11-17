using LaTeXStrings

df = DataFrame(CSV.File("data/M_spin_hist.csv"))
hist_mm = Matrix(df)
histogram(hist_mm,show =true, xlab = "M", 
ylab = "Counts",linecolor = nothing
,label = ["1Spin" "2Spins" "3Spins" "4Spins" "5Spins"],
xscale = :log, yscale = :lin)
