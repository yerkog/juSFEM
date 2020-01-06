> 如下图所示的悬臂梁，受载荷作用如图。$E＝2.1×10^5 N/mm^2$，$μ＝0.3mm$，厚度$h＝10mm$。现用有限元法分析其位移及应力。梁可视为平面应力状态，先按图示尺寸划分为均匀的三角形网格，共有 $8×10＝80$个单元，$5×11=55$个节点，坐标轴以及单元与节点的编号如图。将均布载荷分配到各相应节点上，把有约束的节点 $51、52、53、54、55$ 视作固定铰链，建立如图所示的离散化计算模型。

![fig.1](https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/1_FEM(2D)/beam.png)

### 程序计算框图：

<div  align="center">    
<img src="https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/1_FEM(2D)/program.jpg" width=380 height=750 />
</div>



## 程序中的函数功能介绍及源代码

### 一、函数介绍

1. $LinearTriangleElementStiffness$

   该函数用于计算平面应力情况下弹性模量为 $E$、泊松比为 $NU$ 、厚度为 $t$ 、第一个节点坐标为 $(x_i,y_i)$ 、第二个节点坐标为 $(x_j,y_j)$ 、第三个节点坐标为 $(x_m，y_m)$ 时的线性三角形元的单元刚度矩阵.该函数返回 $6×6$的单位刚度矩阵 $ke$
2. $LinearTriangleAssemble$

   该函数将连接节点 $i,j,m$ 的线性三角形元的单元刚度矩阵$k$集成到整体刚度矩阵$K$。每集成一个单元，该函数都将返回 $2N×2N$ 的整体刚度矩阵 $KE$
3. $LinearTriangleElementStresses$

   该函数计算在平面应力情况下弹性模量为 $E$、泊松比为 $NU$、第一个节点坐标为 $(x_i,y_i)$ ,第二个节点坐标为 $(x_j,y_j)$ ,第三个节点坐标为 $(x_m，y_m)$ 以及单元位移矢量为 $u$ 时的单元应力。该函数返回单元应力矢量。

### 二、悬臂梁轴线位移图
![fig.3](https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/1_FEM(2D)/disp.png)
### 三、悬臂梁应力云图
![fig.3](https://raw.githubusercontent.com/Elliothuo/Storage-Sac/master/images/1_FEM(2D)/stress.png)
