using PlotlyJS
using ORCA

p1 = Plot(scatter(;y=randn(3)))
p2 = Plot(histogram(;x=randn(50), nbinsx=4))
p3 = Plot(scatter(;y=cumsum(randn(12)), name="Random Walk"))
p4 = Plot([scatter(;x=1:4, y=[0, 2, 3, 5], fill="tozeroy"),
           scatter(;x=1:4, y=[3, 5, 1, 7], fill="tonexty")])
plot([p1 p2
      p3 p4])
#p = [p1 p2]
#savefig(p, "G://Project-SFEM//PlotlyJS//plot1.eps")
