# `JuliaLang` FEM & S-FEM Notes 📔


  ```shell
  Microsoft Windows
  (c) 2018 Microsoft Corporation 
                 _
     _       _ _(_)_     |  Documentation: https://docs.julialang.org
    (_)     | (_) (_)    |
     _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
    | | | | | | |/ _` |  |
    | | |_| | | | (_| |  |  Version 1.1.0 (2019-01-21)
   _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
  |__/                   |
  ```

### Ⅰ. Julia ENV
<br>

Windows : `cmd`
```shell
set http_proxy=http://127.0.0.1:1080
set https_proxy=http://127.0.0.1:1080
```
<br>

① **VScode** ⭐⭐⭐

| Extensions | Function | Note |
| :------: | :------: | :------: |
|<a href="https://github.com/JuliaEditorSupport/julia-vscode" target="_blank">JuliaEditorSupport</a>| julia | - |
| Markdown Julia | Julia syntax highlighting in Markdown fenced code block | - |
| Markdown Preview Enhanced | markdown | - |
| Python | run python | - |
| vscode-pdf | view pdf | - |
> Copy extension's name and search in `VScode->EXTENSIONS`

<br>

② **Jupyter** ⭐⭐⭐
> Locate`julia` and `jupyter` path first
  ```julia
  julia> ENV["JUPYTER"]="jupyter"
  julia> using Pkg
  julia> Pkg.add("IJulia")
  julia> using IJulia
  ```
<br>

③ **Atom** ⭐⭐⭐⭐⭐
  ```julia
  julia> using Pkg
  julia> Pkg.add("Atom")
  ```
| Atom pkg name | Function |
| :------: | :------: |
| `Atom`中安装`JunoLab`的内容相关,`ink` | 运行`Julia` |
| Activate Power Mode | 打字效果 |
| Script | run langs |
| minimap | full code |
| tool-bar | tool bar （julia） |
| PDF View | 查看PDF |

| Atom Themes | Name |
| :------: | :------: |
| Monokai | - |
| One Dark | - |
<br>

### Ⅱ.Julia PKGS
<br>

| Pkgname | Function | Note |
| :------: | :------: | :------: |
| <a href="https://github.com/JuliaPy/PyPlot.jl" target="_blank">PyPlot</a> | 绘图包 | - |
| <a href="https://github.com/JuliaData/CSV.jl" target="_blank">CSV</a> | 读写`CSV`文件 | `df=DataFrame(var)`<br>`CSV.write("name.csv",df)` |
| <a href="https://github.com/JuliaData/DataFrames.jl" target="_blank">DataFrames</a> | 数据管理包 | - |
| <a href="https://github.com/JuliaDiffEq/DifferentialEquations.jl" target="_blank">DifferentialEquations</a> | 求解微分方程 | - |
| <a href="https://github.com/JuliaDiffEq/DifferentialEquations.jl" target="_blank">PyCall</a> | 调用`Python` | `@pyimport` |
| <a href="https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html" target="_blank">LinearAlgebra</a> | 调用线性代数 | - |
| <a href="https://github.com/JuliaCI/BenchmarkTools.jl" target="_blank">BenchmarkTools</a> | 基准测试 | `Pkg.add("BenchmarkTools")` |
| <a href="https://github.com/sunoru/RandomNumbers.jl" target="_blank">RandomNumbers</a> | Random Number Generators for the Julia Language | `Pkg.add("RandomNumbers")` |
| <a href="https://github.com/JuliaGizmos/Interact.jl" target="_blank">Interact</a> | Interactive widgets to play with your Julia code | - |
| <a href="https://github.com/JuliaMath/IterativeSolvers.jl" target="_blank">IterativeSolvers</a> | Iterative algorithms for solving linear systems, eigensystems, and singular value problems | - |
| <a href="https://github.com/symengine/SymEngine.jl" target="_blank">SymEngine</a> | IJulia wrappers of SymEngine | - |
| <a href="https://github.com/konsim83/TriangleMesh.jl" target="_blank">TriangleMesh</a> | Generate and refine unstructured 2D triangular meshes from polygons with Julia | - |
| <a href="https://github.com/giordano/Cuba.jl" target="_blank">Cuba</a> | Library for multidimensional numerical integration with four independent algorithms: Vegas, Suave, Divonne, and Cuhre | - |
