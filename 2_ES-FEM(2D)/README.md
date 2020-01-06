### Ⅰ. $S-FEM$ 分类
| $S-FEM$类型 |           对应名称           |
| :---------: | :--------------------------: |
|  $NS-FEM$   |      光滑节点域有限元法      |
|  $CS-FEM$   |     光滑子单元域有限元法     |
|  $ES-FEM$   | 光滑边域有限元法（ $2D$ 优） |
|  $FS-FEM$   | 光滑面域有限元法（ $3D$ 优） |

---
<br>

### Ⅱ. $ES-FEM$ 公式推导

​        对于光滑有限元模型，在假设位移场和计算相容应变场后，需要修正相容应变场，或者直接使用假设位移场来构建一个新的应变场，而不计算原来的相容应变场。用修正或直接构建的应变场来计算应变势能，并且用一个合适的能量弱形式来构造离散模型。

#### 2.1 位移场表达

##### 2.1.1 光滑域划分

​        对于 $ES-FEM$，光滑域与单元的边 $k$ 相关，基于边的光滑域是通过连接边的两个顶点及邻接单元中心点创建的。

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/10.png width=900 height=450 />

##### 2.1.2 形函数的计算

**方法一：**与 $FEM$ 构建方法相同

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/11.png width=600 height=170 />

**方法二：** 线性 $PIM（Point Interpolation Method）$ —— 更通用的方法

​        主要针对多边形构建形函数，$Dai$ 在设计一个简洁但很重要的模型，称作线性 $PIM$ ，计算普通多边形单元形函数的值。（超出范围暂不介绍）


##### 2.1.3 高斯点处形函数计算

​        当使用沿着光滑域边界的相容（连续）位移场时，光滑应变矩阵计算时只需要光滑域边界线段的中点（ $Gauss Points$ ）的形函数值。在每个高斯点处的形函数值可以通过对包含高斯点的线段的两个端点作简单的线性 $PIM$ 来计算。例如，如图所示线段上的点`#g1` 的形函数值可以用线段两个端点 `#1` 和 `#A` 的形函数值的平均值来表示。因此，为了便于计算 ES-FEM 中高斯点的形函数值，我们需要先计算线段端点例如域节点（( `#1` ,  `#6`  等）和形心（`#A`, `#B` 等）处的形函数值。

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/1.png width=600 height=370 />

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/2.png width=980 height=445 />



#### 2.2 应变场矩阵

​        当相容应变场可以方便获得时，如 $T3$ 或 $T4$ 单元，通过适当修正相容应变场可以得到光滑应变场。另外，光滑应变场也可以在光滑操作下直接用位移场构造。这种光滑操作是基于单元网格上建立的一组光滑域。

##### 2.2.1 基于修正的方法

​        当使用线性三角形单元网格时，光滑应变矩阵 $B_I$ 可以通过下面简单的式来组装：
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/3.png width=250 height=80 />
</div>

​        其中 $n_k^e$ 是边 $k$ 周围单元的数目，$A_j^e$ 为边 $k$ 周围的第 $j$ 个单元的面积；$B_j^e$为边 $k$ 周围的第 $j$ 个单元的相容应变矩阵。注意到在这个公式中，$ES-FEM$ 的系统刚度矩阵的计算只需要用面积和标准 $FEM$ 三角形单元的相容应变矩阵 $B_j^e$ 来计算。公式虽然很简单，但是只适用于线性插值的三角形单元$（T3）$。

##### 2.2.2 基于重构的方法
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/4.png width=600 height=80 />
</div>

这个光滑应变矩阵 $B_I$ 通过以下公式计算：
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/5.png width=400 height=120 />
</div>

上面的方程可以进一步简化为一个求和公式：
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/6.png width=500 height=100 />
</div>

其中 $X_p^G$ 是高斯点，它的长度和单位外法向量分别记作 $l_p$ 和 $n_{h,p}$ 。



#### 2.3 建立离散线性代数方程组

​        我们考虑定义的固体力学问题，并使用 $S-FEM$ 静态分析的一般公式，$ES-FEM$ 的线性系统方程的形式为：
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/7.png width=180 height=70 />
</div>

其中 $K^{ES-FEM} $ 是光滑刚度矩阵，其具体形式如下：
<div  align="center">    
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/8.png width=430 height=90 />
</div>

##### 2.3.1 组装方法

​        除了 $S-FEM$ 的求和是在光滑域而非单元上实施的之外**（ 按边循环 ）**，其他的方式和标准的 $FEM$ 中相似。如同在 $FEM$ 中一样，只要将节点适当编号， $K$ 可成为带状矩阵。 $K$ 的带宽取决于构成光滑域的单元节点编号的最大差值。

##### 2.3.2 本质边界条件施加

        求解过程和标准 $FEM$ 中相同（使用解线性方程组的常规方法）。当网格加密后，$S-FEM$ 的解趋于带有平方可积外力的原始问题的精确解。获得节点位移后，我们可以通过计算位移场、再得到光滑应变场，然后根据本构方程求出应力场。最后可以使用应力、应变在整个问题域上积分得到固体的应变能的解。

### Ⅲ. 算例实践

#### 3.1 算例模型
<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/9.png width=900 height=320 />
        如下图所示的悬臂梁，受载荷作用如图。$E＝2.1×10^5(N/mm^2)$ ，$μ＝0.3$ ，厚度$h＝10(mm)$ 。现用有限元法分析其位移及应力。梁可视为平面应力状态，先按图示尺寸划分为均匀的三角形网格，共有 $8×10＝80$ 个单元， $5×11=55$ 个节点，坐标轴以及单元与节点的编号如图。将均布载荷分配到各相应节点上，把有约束的节点 $51、52、53、54、55$ 视作固定铰链，建立如图所示的离散化计算模型。

#### 3.2 程序流程图

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/program.jpg width=300 height=640 />

#### 3.3 程序说明

 ① $Vector3.jl$ :

通过顺时针输入三角形三节点坐标计算出每条边的外法向单位向量，$Vector4.jl$ 的计算原理相同，需要逆时针输入坐标点即可。

```julia
function vector3(x1,y1,x2,y2,x3,y3)
    A = y2-y1
    B = x1-x2
    v1 = [A/sqrt(A^2+B^2) B/sqrt(A^2+B^2)]
    D = y3-y2
    E = x2-x3
    v2 = [D/sqrt(D^2+E^2) E/sqrt(D^2+E^2)]
    F = y1-y3
    G = x3-x1
    v3 = [F/sqrt(F^2+G^2) G/sqrt(F^2+G^2)]
    return [v1 v2 v3]
end
```

<br>

② $BInner.jl$ :

当为内部边时，光滑域 $\Omega_k^s$ 是四边形，支持 $\Omega_k^s$ 的是一个四边形。按照内部边的端点（两个中的任意一个均可）为起点，逆时针输入四个点的坐标，按照公式计算每个点的光滑应变矩阵，最后组装成“光滑域单元刚度矩阵”（此处的单元和三角形单元不同），例如 $B_Q = [B_{I1},B_{I2},B_{I3},B_{I4}](3×6)$ 。高斯点形函数值的计算见理论推导。$BOuter.jl$ 的计算类似。

```julia
include("Vector4.jl")
function BQ(x1,y1,x2,y2,x3,y3,x4,y4)
    xi1 = (x1+x2+x3)/3
    yi1 = (y1+y2+y3)/3
    xi2 = (x1+x3+x4)/3
    yi2 = (y1+y3+y4)/3
     lp = [sqrt((x1-xi1)^2+(y1-yi1)^2) sqrt((xi1-x3)^2+(yi1-y3)^2) sqrt((x3-xi2)^2+(y3-yi2)^2) sqrt((xi2-x1)^2+(yi2-y1)^2)]
      n = vector4(x1,y1,xi1,yi1,x3,y3,xi2,yi2)
    Gp1 = [4/6 1/6 1/6 4/6]
    b11 = n[1]*Gp1[1]*lp[1]
    b12 = n[3]*Gp1[2]*lp[2]
    b13 = n[5]*Gp1[3]*lp[3]
    b14 = n[7]*Gp1[4]*lp[4]
    bx1 = (3/0.0008)*(b11+b12+b13+b14)
    b15 = n[2]*Gp1[1]*lp[1]
    b16 = n[4]*Gp1[2]*lp[2]
    b17 = n[6]*Gp1[3]*lp[3]
    b18 = n[8]*Gp1[4]*lp[4]
    by1 = (3/0.0008)*(b15+b16+b17+b18)
     B1 = [bx1 0; 0 by1; by1 bx1]
    Gp2 = [1/6 1/6 0 0]
    b21 = n[1]*Gp2[1]*lp[1]
    b22 = n[3]*Gp2[2]*lp[2]
    b23 = n[5]*Gp2[3]*lp[3]
    b24 = n[7]*Gp2[4]*lp[4]
    bx2 = (3/0.0008)*(b21+b22+b23+b24)
    b25 = n[2]*Gp2[1]*lp[1]
    b26 = n[4]*Gp2[2]*lp[2]
    b27 = n[6]*Gp2[3]*lp[3]
    b28 = n[8]*Gp2[4]*lp[4]
    by2 = (3/0.0008)*(b25+b26+b27+b28)
     B2 = [bx2 0; 0 by2; by2 bx2]
    Gp3 = [1/6 4/6 4/6 1/6]
    b31 = n[1]*Gp3[1]*lp[1]
    b32 = n[3]*Gp3[2]*lp[2]
    b33 = n[5]*Gp3[3]*lp[3]
    b34 = n[7]*Gp3[4]*lp[4]
    bx3 = (3/0.0008)*(b31+b32+b33+b34)
    b35 = n[2]*Gp3[1]*lp[1]
    b36 = n[4]*Gp3[2]*lp[2]
    b37 = n[6]*Gp3[3]*lp[3]
    b38 = n[8]*Gp3[4]*lp[4]
    by3 = (3/0.0008)*(b35+b36+b37+b38)
     B3 = [bx3 0; 0 by3; by3 bx3]
    Gp4 = [0 0 1/6 1/6]
    b41 = n[1]*Gp4[1]*lp[1]
    b42 = n[3]*Gp4[2]*lp[2]
    b43 = n[5]*Gp4[3]*lp[3]
    b44 = n[7]*Gp4[4]*lp[4]
    bx4 = (3/0.0008)*(b41+b42+b43+b44)
    b45 = n[2]*Gp4[1]*lp[1]
    b46 = n[4]*Gp4[2]*lp[2]
    b47 = n[6]*Gp4[3]*lp[3]
    b48 = n[8]*Gp4[4]*lp[4]
    by4 = (3/0.0008)*(b45+b46+b47+b48)
     B4 = [bx4 0; 0 by4; by4 bx4]
      B = [B1 B2 B3 B4]
    return B
end
```

<br>

③ $Stiffness.jl$ :

负责将上一步计算得到的光滑应变矩阵进行刚度矩阵计算，内部边的刚度矩阵得到的是 $(8\times8)$ 的矩阵，外部边的刚度矩阵是 $(6\times6)$ 。（$D$ 为弹性系数矩阵，此问题为平面应力问题）

```julia
function Stiffness(B,A)
    D = (E/(1-NU*NU))*[1 NU 0; NU 1 0; 0 0 (1-NU)/2]
    ke = t*A*B'*D*B
    return ke
end
#test pass
```

<br>

④ $AssembleQuad.jl$ :

$AssembleQuad.jl$ 是负责将内部边构成的光滑刚度矩阵组装到整体刚度矩阵的函数，编写和 $FEM$ 中用到的组装函数思路完全一致，$AssembleTriangle.jl$ 负责组装内部边组成的光滑刚度矩阵。

```julia
function AssembleTriangle(KE,ke,i,j,m)
    KE[2*i-1,2*i-1] = KE[2*i-1,2*i-1] + ke[1,1]
    KE[2*i-1,2*i] = KE[2*i-1,2*i] + ke[1,2]
    KE[2*i-1,2*j-1] = KE[2*i-1,2*j-1] + ke[1,3]
    KE[2*i-1,2*j] = KE[2*i-1,2*j] + ke[1,4]
    KE[2*i-1,2*m-1] = KE[2*i-1,2*m-1] + ke[1,5]
    KE[2*i-1,2*m] = KE[2*i-1,2*m] + ke[1,6]
    KE[2*i,2*i-1] = KE[2*i,2*i-1] + ke[2,1]
    KE[2*i,2*i] = KE[2*i,2*i] + ke[2,2]
    KE[2*i,2*j-1] = KE[2*i,2*j-1] + ke[2,3]
    KE[2*i,2*j] = KE[2*i,2*j] + ke[2,4]
    KE[2*i,2*m-1] = KE[2*i,2*m-1] + ke[2,5]
    KE[2*i,2*m] = KE[2*i,2*m] + ke[2,6]
    KE[2*j-1,2*i-1] = KE[2*j-1,2*i-1] + ke[3,1]
    KE[2*j-1,2*i] = KE[2*j-1,2*i] + ke[3,2]
    KE[2*j-1,2*j-1] = KE[2*j-1,2*j-1] + ke[3,3]
    KE[2*j-1,2*j] = KE[2*j-1,2*j] + ke[3,4]
    KE[2*j-1,2*m-1] = KE[2*j-1,2*m-1] + ke[3,5]
    KE[2*j-1,2*m] = KE[2*j-1,2*m] + ke[3,6]
    KE[2*j,2*i-1] = KE[2*j,2*i-1] + ke[4,1]
    KE[2*j,2*i] = KE[2*j,2*i] + ke[4,2]
    KE[2*j,2*j-1] = KE[2*j,2*j-1] + ke[4,3]
    KE[2*j,2*j] = KE[2*j,2*j] + ke[4,4]
    KE[2*j,2*m-1] = KE[2*j,2*m-1] + ke[4,5]
    KE[2*j,2*m] = KE[2*j,2*m] + ke[4,6]
    KE[2*m-1,2*i-1] = KE[2*m-1,2*i-1] + ke[5,1]
    KE[2*m-1,2*i] = KE[2*m-1,2*i] + ke[5,2]
    KE[2*m-1,2*j-1] = KE[2*m-1,2*j-1] + ke[5,3]
    KE[2*m-1,2*j] = KE[2*m-1,2*j] + ke[5,4]
    KE[2*m-1,2*m-1] = KE[2*m-1,2*m-1] + ke[5,5]
    KE[2*m-1,2*m] = KE[2*m-1,2*m] + ke[5,6]
    KE[2*m,2*i-1] = KE[2*m,2*i-1] + ke[6,1]
    KE[2*m,2*i] = KE[2*m,2*i] + ke[6,2]
    KE[2*m,2*j-1] = KE[2*m,2*j-1] + ke[6,3]
    KE[2*m,2*j] = KE[2*m,2*j] + ke[6,4]
    KE[2*m,2*m-1] = KE[2*m,2*m-1] + ke[6,5]
    KE[2*m,2*m] = KE[2*m,2*m] + ke[6,6]
    return KE
end
#test pass
```

<br>

⑤ $beam-main.jl$ :

主函数主要负责调用以上所有函数，按照按边循环的思想，根据算例输入特定条件，最后注释掉的部分为调用 $Gadfly.jl$ 进行作图的代码。

```julia
include("AssembleQuad.jl")
include("AssembleTriangle.jl")
include("BInner.jl")
include("BOuter.jl")
include("Stiffness.jl")

using LinearAlgebra
using Gadfly
using Cairo
using Fontconfig

E = 21e10
NU = 0.3
t = 0.01

KE = zeros(110,110)

edge1 = BT(0.4,0.08,0.36,0.08,0.36,0.06)
s1 = Stiffness(edge1,0.0008/6)
for j = 1:5:46
    global KE = AssembleTriangle(KE,s1,j,j+5,j+6)
end
edge2 = BT(0.36,0,0.4,0,0.4,0.02)
s2 = Stiffness(edge2,0.0008/6)
for j = 10:5:55
    global KE = AssembleTriangle(KE,s2,j,j-5,j-6)
end
edge3 = BT(0.4,0.06,0.4,0.08,0.36,0.06)
s3 = Stiffness(edge3,0.0008/6)
for j = 2:1:5
    global KE = AssembleTriangle(KE,s3,j,j-1,j+5)
end
edge4 = BT(0,0.08,0,0.06,0.04,0.08)
s4 = Stiffness(edge4,0.0008/6)
for j = 51:1:54
    global KE = AssembleTriangle(KE,s4,j,j+1,j-5)
end

edge5 = BQ(0.4,0.08,0.36,0.08,0.36,0.06,0.4,0.06)
s5 = Stiffness(edge5,0.0008/3)
for j = 1:5:46
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 2:5:47
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 3:5:48
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
for j = 4:5:49
    global KE = AssembleQuad(KE,s5,j,j+5,j+6,j+1)
end
edge6 = BQ(0.36,0.06,0.36,0.04,0.4,0.06,0.4,0.08)
s6 = Stiffness(edge6,0.0008/3)
for j = 7:5:52
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
for j = 8:5:53
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
for j = 9:5:54
    global KE = AssembleQuad(KE,s6,j,j+1,j-5,j-6)
end
edge7 = BQ(0.36,0.08,0.32,0.06,0.36,0.06,0.4,0.08)
s7 = Stiffness(edge7,0.0008/3)
for j = 6:5:46
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 7:5:47
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 8:5:48
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end
for j = 9:5:49
    global KE = AssembleQuad(KE,s7,j,j+6,j+1,j-5)
end


TEMPQ = zeros(110,1)
TEMPQ[2,1] = -4000

TEMPMatrix = zeros(110,110)
TEMPMatrix[1:100,1:100] = KE[1:100,1:100]
TEMPMatrix[101:110,101:110] = Matrix{Float16}(I,10,10)

d=TEMPMatrix\TEMPQ
#=
XX = collect(0:0.04:0.4)
Y1 = d[106:-10:6]
Y2 = [0.0 -1.0822437722258944e-5 -4.035503570515363e-5 -8.588051757808588e-5 -0.00014535720223828857 -0.00021680303958405118 -0.0002982248351603411 -0.00038762920996542776 -0.0004830328435555055 -0.0005824293893356926 -0.0006833792103817733]
Y3 = map(i -> (((-4000*i^2)/(6*E*(t*0.08^3/12)))*(3*0.4-i)),0:0.04:0.4;)
line1 = layer(x=XX,y=Y1,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="green"))
line2 = layer(x=XX,y=Y2,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="blue"))
line3 = layer(x=XX,y=Y3,Geom.point,Geom.smooth(method=:loess,smoothing=0.8),Theme(default_color="red"))

plot(x=XX,y=YY,Geom.point,Geom.smooth(method=:loess,smoothing=0.9),
Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"))
set_default_plot_size(20cm,8cm)
Gadfly.with_theme(:default) do
f = plot(line1,line2,line3,Guide.manual_color_key("Methods", ["ES-FEM", "FEM", "analytical solution"],["green", "blue", "red"]),Guide.ylabel("Beam Deflection (m)"),Guide.xlabel("Beam Axial Distance (m)"),Guide.title("RESULT"))
    #,Coord.cartesian(xmin=0, xmax=0.18)
    draw(PDF("foo.pdf",20cm,8cm,dpi=1000),f)
end
=#
```

<br>

#### 3.4 结果展示

<img src = https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/2_ES-FEM(2D)/compare.jpg />

<br>

### Ⅳ. 书中结果对比方法总结

① 例 7.7.1 矩形悬臂梁(静态分析)
在实例 5.8.1 中描述的承受端部载荷的矩形悬臂梁问题被用来验证 $ES-FEM$ 模型。
矩形悬臂梁的几何形状和边界条件可参看图 5.6。图 5.7 和 6.6 分别展示了用四边形、多
边形和三角形单元离散的区域。问题的精确应变能为 $4.4746Nm$。结果有：

```
自由边受抛物型切向力的矩形悬臂梁,沿横坐标轴的位移分布:
自由边受抛物型切向力的矩形悬臂梁,沿横坐标轴的相对位移误差
使用 ES-FEM-T3 的自由边受抛物型切向力的矩形悬臂梁,正应力剪应力
自由边受抛物型切向力的矩形悬臂梁,不同方法的 CPU 计算时间。使用相同节点时 FEM-T3是最快的
应力云图
```

② 带圆孔的无限平板(静态分析)
```
对比内容类似
```

③ 四面体悬臂梁
```
对比内容类似
```
