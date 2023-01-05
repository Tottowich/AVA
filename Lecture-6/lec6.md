# Lecture 6
* **Divergence**
<!-- Definition Box -->
<div class="DefinitionBox" style="border: 2px solid Blue; padding: 10px 5px">

The **divergence** of a continously differentiable vector field $\bar{A}$, then the divergence of $\bar{A}$ is a continously differentiable scalar field, and its value at the point $p$ is defined as:
$$(div\bar{A})(p)=\lim_{\Delta V\rightarrow0} \frac{1}{\Delta V}\oiint_S d\bar{S}\cdot\bar{A}$$
where $\Delta V$ is the infinitesimal volume around the given point $p$, and $d\bar{S}$ is the perpendicular infinitesimal surface element directed outward from the surface $S$.
</div>
<br>
<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Orange; padding: 10px 
5px">

## Theorem of Divergence
The function is more commonly.
$$div\bar{A}=\nabla\cdot\bar{A}=\frac{\partial}{\partial x_1}\bar{A}_1+\frac{\partial}{\partial x_2}\bar{A}_2+\frac{\partial}{\partial x_3}\bar{A}_3$$
</div>
<br>
<!-- Proof box -->
<div class="ProofBox" style="border: 2px solid Magenta; padding: 10px 5px">

Let:
$$ p = \{x_0,y_0,z_0\} $$

Then the divergence can be written as.

$$ \oiint_S d\bar{S}\cdot\bar{A} = F_x+F_y+F_z $$
$$ F_x = \iint_{S_1} dS\bar{A}_x + \it_{S_2} dS(-\bar{A}_x)$$
$S_x$ is the projection of $S_1$ and $S_2$ on the $y-z$ plane.
Then:
$$F_x = \iint_{S_x}dydz\left (A_x(x_0+ \frac{\Delta x}{2},y,z)-A_x(x_0-\frac{\Delta x}{2},y,z)\right )$$
Using the mean value theorem we recieve:
$$F_x = \left (A_x(x_0+ \frac{\Delta x}{2},y,z)-A_x(x_0-\frac{\Delta x}{2},y,z)\right )\Delta y\Delta z$$
Once again using the mean value theorem we recieve:
$$F_x = \frac{\partial A_x}{\partial x}(x_0+\theta \frac{\Delta x}{2},y,z)\Delta x\Delta y\Delta z= \frac{\partial A_x}{\partial x}(p')\Delta V$$
where $|\theta|<1$.
When the box $p'\rightarrow p$ we get:
$$\lim_{\Delta V \rightarrow 0} = \frac{1}{\Delta V} F_x= \lim_{p'\rightarrow p}\frac{\partial A_x}{\partial x}(p')=\frac{\partial A_x}{\partial x}(p)$$
This could also be done for $F_y$ and $F_z$.
Which gives us:
$$(div\bar{A})_p = \frac{\partial A_x}{\partial x}(p)+\frac{\partial A_y}{\partial y}(p)+\frac{\partial A_z}{\partial z}(p)$$
</div>
<br>

* **Nabla operator** $\nabla = \frac{\partial}{\partial x_1}\hat{x}+\frac{\partial}{\partial x_2}\hat{y}+\frac{\partial}{\partial x_3}\hat{z}$
Using the nabla operator we can express the divergence as:
$$div\bar{A}=\nabla\cdot\bar{A}=\frac{\partial}{\partial x_1}\bar{A}_1+\frac{\partial}{\partial x_2}\bar{A}_2+\frac{\partial}{\partial x_3}\bar{A}_3+...+\frac{\partial}{\partial x_n}\bar{A}_n$$
* **Curl operator** $\nabla\times\bar{A}$ Will be discussed in the later chapters.
* **Gradient** could also be expressed as 
$$\nabla f = \frac{\partial f}{\partial x_1}\hat{x}+\frac{\partial f}{\partial x_2}\hat{y}+\frac{\partial f}{\partial x_3}\hat{z}$$
<!-- Proof box -->
<div class="ProofBox" style="border: 2px solid Magenta; padding: 10px 5px">

## Gauss Theorem
$$ \oiint_S d\bar{S}\cdot\bar{A} = \iiint_V d\bar{V}\;\nabla\cdot\bar{A}$$
Consider first volume with a $S$ that intersects lines parallel to the $x,y$ or $z$ axis twice at most, i.e. a simple body. The body can not be folded onto itself even without contact. This senario will be called a **simple body**.

$$ \iiint_{dV}dV\;\nabla\cdot\bar{A} = \iiint dxdydz\left(\frac{\partial A_x}{\partial x}+\frac{\partial A_y}{\partial y}+\frac{\partial A_z}{\partial z}\right)$$

* Focus on the $\frac{\partial A_z}{\partial z}$ term for now, the other terms are similar.
$$\iiint dxdydz\frac{\partial A_z}{\partial z} =\iint_{S_p} dxdy\int_{f(x,y)}^{g(x,y)}dz$$
The first term could be thought of as the integral of the $S_p$ which is the projection of the volume onto the $x-y$ plane.

Using $S_p$: projection onto x-y plane
$$ \int_{S_p}dxdy\int_{g(x,y)}^{f(x,y)}dz\frac{\partial A_z}{\partial z}$$
$$ = \int_{S_p}dxdy\left[A_z(x,y,f(x,y)-A_z(x,y,g(x,y)\right]$$
$$ \bar{r}_2 = x\bar{e}_x+y\bar{e}_y+f(x,y)\bar{e}_z$$
$$ d\bar{S}_2 = \left(\frac{\partial \bar{r}}{\partial x}\times\frac{\partial \bar{r}}{\partial y}\right)dxdy=\left(-\frac{\partial f}{\partial x}\bar{e}_x-\frac{\partial f}{\partial y}\bar{e}_y+\bar{e}_z\right)dxdy$$
$$\bar{e}_z\cdot d\bar{S}_2 = dxdy$$
FILL THIS IN
<br>
<br>
$$\iiint dV\;\nabla\cdot\bar{A} = \oiint_S d\bar{S}\cdot\left(A_z\bar{e}_z+\bar{A}_x\bar{e}_x+\bar{A}_y\bar{e}_y\right)$$
</div>
<br>

Any volume can be divided into subvolumes for which the previous assumtions hold.
$$\iiint dV\;\nabla\cdot\bar{A} = \sum_{i=1}^n\iiint_{V_i}dV\;\nabla\cdot\bar{A}=\sum_{i=1}^n\oiint_{S_i}d\bar{S}\cdot\left(A_z\bar{e}_z+\bar{A}_x\bar{e}_x+\bar{A}_y\bar{e}_y\right)$$
$$ = \oiint_S d\bar{S}\cdot\bar{A}$$

<!-- Example box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## Applications fluid of model.
The fluid flows at a velocity $\bar{v}$ through the surface $S$.
$$\oiint d\bar{S}\cdot h\bar{v}$$
Where $h(\bar{r})$ is the density of the fluid at $\bar{r}$ and is a scalar filed. And $\bar{v}$ is a vector field.
Using Gauss Theorem we get:
$$\oiint_S d\bar{S}\cdot h\bar{v}=\iiint dV\;\nabla\cdot h\bar{v}$$
Assume that $h$ is constant in the volume $V$, incompressible fluid and that there are no sources or sinks in the volume.
$$\int_V dV\;\nabla\cdot h\bar{v} = h\int_V dV\;\nabla\cdot\bar{v}=\oiint d\bar{S}\cdot h\bar{v}=0$$
But if h is not constant in the volume $V$ we get:
$$\frac{d}{dt}\int dVh = \int dV\frac{\partial u}{\partial t}