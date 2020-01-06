using PlotlyJS
using ORCA

function bar9()
    trace1 = bar(;x=["0.56334","1.10847","1.57959","2.0775"],
                  y=[3.11,5.77,8.09,10.03],
                  name="Multi-core Parallel",
                  marker_color="rgb(55, 83, 109)")
    trace2 = bar(;x=["0.56334","1.10847","1.57959","2.0775"],
                  y=[54.53,112.83,157.67,208.67],
                  name="Single-core ",
                  marker_color="rgb(26, 118, 255)")

    data = [trace1, trace2]

    layout = Layout(;xaxis=attr(tickfont_size= 14,
                                tickfont_color="rgb(107, 107, 107)"),
                     yaxis=attr(title="USD (millions)",
                                titlefont=attr(size=16,
                                               color="rgb(107, 107, 107)"),
                                tickfont=attr(size=14,
                                              color="rgb(107, 107, 107)")),
                     legend=attr( bgcolor="rgba(255, 255, 255, 0)",
                                 bordercolor="rgba(255, 255, 255, 0)"),
                     barmode="group",
                     #bargap=0.15,
                     #bargroupgap=0.1)
                     )
    plot(data, layout)
end
bar9()
