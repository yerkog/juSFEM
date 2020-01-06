include("LinearTriangleAssemble.jl")
include("LinearTriangleElementStiffness.jl")
include("LinearTriangleElementStresses.jl")

using LinearAlgebra
using Gadfly

#E = 2.06e11
E = 21e10
NU = 0.3
t = 0.01

stifflike1 = LinearTriangleElementStiffness(E,NU,t,0.4,0.08,0.36,0.06,0.4,0.06)
stifflike5 = LinearTriangleElementStiffness(E,NU,t,0.4,0.08,0.36,0.08,0.36,0.06)

KE = zeros(110,110)
for i = 1:49
    if i%5 != 0
    j=i
    global KE=LinearTriangleAssemble(KE,stifflike5,j,j+5,j+6)
    global KE=LinearTriangleAssemble(KE,stifflike1,j,j+6,j+1)
    end
end

TEMPQ = zeros(110,1)
TEMPQ[2,1] = -4000
#=TEMPQ[2,1] = -200
for i = (12:10:92)
   global TEMPQ[i,1] = -400
end
=#

TEMPMatrix = zeros(110,110)
TEMPMatrix[1:100,1:100] = KE[1:100,1:100]
TEMPMatrix[101:110,101:110] = Matrix{Float16}(I,10,10)

d=TEMPMatrix\TEMPQ

y = zeros(80,3)
q = 0

for i = 1:49
   if i%5 == 1
      j=2*i
      u=[d[j-1] d[j] d[j+11] d[j+12] d[j+1] d[j+2]]
      u=u'
      xl = 0.40
      yl = 0.08
      xm = 0.36
      ym = 0.06
      xn = 0.40
      yn = 0.06
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 2
      j=2*i
      u=[d[j-1] d[j] d[j+11] d[j+12] d[j+1] d[j+2]]
      u=u'
      xl=0.40
      yl=0.06
      xm=0.36
      ym=0.04
      xn=0.40
      yn=0.04
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 3
      j=2*i
      u=[d[j-1] d[j] d[j+11] d[j+12] d[j+1] d[j+2]]
      u=u'
      xl=0.40
      yl=0.04
      xm=0.36
      ym=0.02
      xn=0.40
      yn=0.02
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 4
      j=2*i
      u=[d[j-1] d[j] d[j+11] d[j+12] d[j+1] d[j+2]]
      u=u'
      xl=0.40
      yl=0.02
      xm=0.36
      ym=0.00
      xn=0.40
      yn=0.00
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   else
      global q=q+3
   end
end

global q=4

for i=1:49
   if i%5 == 1
      j=2*i
      u=[d[j-1] d[j] d[j+9] d[j+10] d[j+11] d[j+12]]
      u=u'
      xl=0.40
      yl=0.08
      xm=0.36
      ym=0.08
      xn=0.36
      yn=0.06
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 2
      j=2*i
      u=[d[j-1] d[j] d[j+9] d[j+10] d[j+11] d[j+12]]
      u=u'
      xl=0.40
      yl=0.06
      xm=0.36
      ym=0.06
      xn=0.36
      yn=0.04
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 3
      j=2*i
      u=[d[j-1] d[j] d[j+9] d[j+10] d[j+11] d[j+12]]
      u=u'
      xl=0.40
      yl=0.04
      xm=0.36
      ym=0.04
      xn=0.36
      yn=0.02
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
   elseif i%5 == 4
      j=2*i
      u=[d[j-1] d[j] d[j+9] d[j+10] d[j+11] d[j+12]]
      u=u'
      xl=0.40
      yl=0.02
      xm=0.36
      ym=0.02
      xn=0.36
      yn=0.00
      y[i+q,:]=LinearTriangleElementStresses(E,NU,xl,yl,xm,ym,xn,yn,u)'
      xl=xl-0.04
      xm=xm-0.04
      xn=xn-0.04
	else
	    global q=q+3
   end
end

N=y[73:80,1]

XX = collect(0:0.04:0.4)
YY = d[106:-10:6]
plot(x=XX,y=YY,Geom.point,Geom.smooth(method=:loess,smoothing=0.5)
,Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"))
