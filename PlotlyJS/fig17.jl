using PlotlyJS
using ORCA

X = [0.5,1,1.5,2]
Y = [17.5075,19.5532,19.4680,20.7903]
p = scatter(;x = X, y = Y, mode = "lines+markers")
layout = Layout(;xaxis_range=[0, 2.4],
                 yaxis_range=[10, 24])
plot(p,layout)
