using Plots


plotlyjs()

x = 1:10
y = rand(10)
plt = plot(x, y)
gui()
savefig(plt,"F://Project-SFEM//Plots//myplot.svg")
