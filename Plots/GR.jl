using Plots

gr()
x = 1:10
y = rand(10)
plt = plot(x, y,
           xlabel = "my label",
           xlims = (0,10),
           xticks = 0:0.5:10,
           xtickfont = font(10, "Times New Roman"),
           ylabel = "my label",
           ylims = (0,1),
           yticks = 0:0.05:1,
           ytickfont = font(10, "Times New Roman"),
           #size = (1280, 920)
          )
display(plt)
savefig(plt,"F://Project-SFEM//Plots//myplot_gr.pdf")
