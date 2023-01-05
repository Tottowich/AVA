# Lecture 3
## Surface integrals:
Flux of a fluid.
* How much fluid is flowing through a surface?
$f$ is the voluime of fluid passing $S$ per unit time.

<div style="text-align:center"><img src="LectureImages\Lecture5\Lecture5FluxInitEx.png" /></div>

$$x = |\bar{v}|\cdot\Delta{t}$$
$$F = \Delta{V}/\Delta{t},\; \Delta{V} = S|\vec{v}|\Delta{t}$$
$$F = \frac{S|\vec{v}|\Delta{t}}{\Delta{t}} = S|\vec{v}|$$


* At an angle the perpendicular component of the velocity does not contribute to the flux.

<!-- Image of anlge with a border-->
<div style="text-align:center"><img src="LectureImages\Lecture5\Lecture5FluxAngleEx.png" /></div>

* Note! To only obtain the parallel component of the velocity we use the dot product.
$$F = S\cdot\bar{v}\cdot\cos{\theta}=S\bar{v}\cdot\hat{n}$$
We could rewrite this as:
$$F = \bar{v}\cdot\bar{S}$$
Where $\bar{S}$ is the area vector of the surface which is perpendicular to the surface.
$$\bar{S} = S\cdot\hat{n}$$
The decomposion of $\bar{v}$ can be as:
$$\bar{v} = (v\cdot\hat{n})\hat{n}+v^{\perp},\; v^{\perp}\cdot\hat{n} = 0$$

For a general surface $S$ the flux is:

$$\bar{F} \approx \sum_{i=1}^{n}F_i=\sum_{i=1}^{n} \Delta \bar{S_i}\cdot\bar{v}_i$$

$$\bar{F} = \lim_{\Delta{S}\rightarrow0}\sum_{i=1}^{n} \Delta \bar{S_i}\cdot\bar{v}_i=\int_{S}\bar{v}\cdot{d\bar{S}}$$
Which can be written as a dubble integral:
<!-- FormulaBox -->
<div class="FormulaBox" style="border: 2px solid Magenta; padding: 10px 5px">

## Flux of a fluid
$$\bar{F} = \int\int\bar{v}\cdot d\bar{S}= \int\int\bar{v}\cdot\hat{n}dS$$

</div>

1. Parametize the surface $S: \bar{r}(u,v)$
2. Express the field as a function of the parameters: $\bar{A} = \bar{A}(u,v)$
3. Express the surface elements $d\bar{S}$ as a function of the parameters: $d\bar{S} = d\bar{S}(u,v)$
4. Performe the double integral: $\bar{F} = \int\int\bar{v}\cdot{d\bar{S}}$

In order to express the surface elements we need to use the cross product between:
<!-- Image of surface element -->
<div style="text-align:center"><img src="LectureImages\Lecture5\Lecture5SurfaceElement.png" /></div>

$$\Delta{r_1} = \bar{r}(u+\Delta{u},v)-\bar{r}(u,v)$$
$$\Delta{r_2} = \bar{r}(u,v+\Delta{v})-\bar{r}(u,v)$$
Then using the cross product we obtain the surface element:
$$\Delta{S} = \Delta{r_1}\times\Delta{r_2}$$
As $\Delta{u}$ and $\Delta{v}$ are infinitesimal we can write:
$$\Delta{r_1} = \frac{\bar{r}(u+\Delta{u},v)-\bar{r}(u,v)}{\Delta{u}}\cdot\Delta{u}=d\bar{r_1} = \frac{\partial\bar{r}}{\partial{u}}\cdot d{u}$$
$$\Delta{r_2} = \frac{\bar{r}(u,v+\Delta{v})-\bar{r}(u,v)}{\Delta{v}}\cdot\Delta{v}=d\bar{r_2} = 
\frac{\partial\bar{r}}{\partial{v}}\cdot d{v}$$
Then to obatin the surface element we use the cross product since the cross product of two vectors generates a vector which size is the area of the parallelogram formed by the two vectors, remember Linear Algebra!

<div class="FormulaBox" style="border: 2px solid Magenta; padding: 10px 5px">

## Surface element
The surface element from a two dimensional surface is:
$$d\bar{S} = d\bar{r_u}\times d\bar{r_v} = \frac{\partial\bar{r}}{\partial{u}}\times\frac{\partial\bar{r}}{\partial{v}}dudv$$

</div>
<br>
<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid green; padding: 10px 5px">

## Ex. Compute the flux of the vector field $\bar{A}$
<hr>

$$\bar{A} = yz^2\bar{e}_x$$
Through the surface $S$:
<!-- Large math bracket -->
$$\bar{r}(x,y) = \begin{cases}x=y^2+z^2\\0<y<1\\0<z<1\end{cases}$$

$$\bar{r} = (u^2+v^2)\bar{e}_x+u\bar{e}_y+v\bar{e}_z$$
With the following limits:
$$u = 0\rightarrow1,\;\; v = 0\rightarrow1$$
Derivative with respecpt to the paramters $u$ and $v$:
$$\frac{\partial\bar{r}}{\partial{v}} = 2v\bar{e}_x+\bar{e}_z$$
$$\frac{\partial\bar{r}}{\partial{u}} = 2u\bar{e}_x+\bar{e}_y$$
Which gives the following surface element using the formula provided above:
$$d\bar{S} = (2u\bar{e}_x+\bar{e}_y)\times(2v\bar{e}_x+\bar{e}_z)du dv$$
Resulting in:
$$d\bar{S} = \left(\bar{e}_x-2u\bar{e}_y-2v\bar{e}_z\right)du dv$$
The vector field is given in terms of $u$ and $v$:
$$\bar{A} = yz^2\bar{e}_x = uv^2\bar{e}_x$$
Then the flux is:
$$\bar{F} = \int\int{d\bar{S}\cdot\bar{A}} = \int{d\bar{S}\cdot uv^2\bar{e}_x}$$
$$\bar{F} = \int\int{\left(\bar{e}_x-2u\bar{e}_y-2v\bar{e}_z\right)\cdot{uv^2\bar{e}_x}du dv}$$
Since $\bar{e}_x\cdot\bar{e}_x = 1$ and $\bar{e}_x\cdot\bar{e}_y = \bar{e}_x\cdot\bar{e}_z = 0$ we can the integral as:
$$\bar{F} = \int_0^1\int_0^1{uv^2du dv}= \frac{1}{2}\frac{1}{3}$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>
<!-- Proof Box -->

<div class="ProofBox" style="border: 2px solid Orange; padding: 10px 5px">

Independent of choice of parametrization:
Proof: 
Consider two parametrizations $\bar{r}(u,v)$ and $\bar{r}(s,t)$.

The surface element in terms of $u$ and $v$:
$$d\bar{S} = \frac{\partial\bar{r}}{\partial{u}}\times\frac{\partial\bar{r}}{\partial{v}}dudv$$
We known that the flux in terms of $u$ and $v$ is:
$$\bar{F} = \int\int{\frac{\partial\bar{r}}{\partial{u}}\times\frac{\partial\bar{r}}{\partial{v}}\cdot\bar{A}}\:dudv$$

If we let $s=s(u,v)$ and $t=t(u,v)$ then via the chain rule:
$$d\bar{S} = \frac{\partial\bar{r}}{\partial{u}}\times\frac{\partial\bar{r}}{\partial{v}}dudv = \left ( \frac{\partial\bar{r}}{\partial s}\frac{\partial s}{\partial u}+\frac{\partial\bar{r}}{\partial t}\frac{\partial t}{\partial u}\right ) \times \left ( \frac{\partial\bar{r}}{\partial s}\frac{\partial s}{\partial v}+\frac{\partial\bar{r}}{\partial t}\frac{\partial t}{\partial v}\right )dudv$$
The terms with the same derivative are parallell and does therefore not affect the cross product.
Then the flux integral can be written as:
$$\bar{F} = \int\int{d\bar{S}\cdot\bar{A}} = \int\int{\frac{\partial\bar{r}}{\partial{s}}\times\frac{\partial{\bar{r}}}{\partial{t}}\cdot{det(J)}}dudv\cdot \bar{A}$$
Where $det(J)$ is the determinant of the Jacobian matrix:
$$J = \begin{bmatrix}\frac{\partial{s}}{\partial{u}} & \frac{\partial{s}}{\partial{v}}\\\frac{\partial{t}}{\partial{u}} & \frac{\partial{t}}{\partial{v}}\end{bmatrix}$$
Which give the following result:
$$ det(J) = \frac{\partial{s}}{\partial{u}}\frac{\partial{t}}{\partial{v}}-\frac{\partial{s}}{\partial{v}}\frac{\partial{t}}{\partial{u}} = \frac{\partial(s,t)}{\partial(u,v)}$$
Using knowledge from multivariable calculus we can rewrite:
$$dudv\frac{\partial(s,t)}{\partial(u,v)}=dsdt$$
This means the at function can be written as:
$$\bar{F} = \int\int{dsdt\frac{\partial\bar{r}}{\partial{s}}\times\frac{\partial{\bar{r}}}{\partial{t}}}$$
Thereby the flux integral is independent of the parametrization.
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

## Post break
There are differnt types of surface integrals:
* $\int_S dS\Phi(\bar{r})$ where $\Phi$ is a scalar function
* $\int_S d\bar{S}\Phi(\bar{r})$ where $\Phi$ is a scalar function
* $\int_S d\bar{S}\cdot\bar{A}(\bar{r})$ where $\bar{A}$ is a vector field
* $\int_S d\bar{S}\times\bar{A}(\bar{r})$ where $\bar{A}$ is a vector field
Note that:
$$dS = |d\bar{S}| = \left|\frac{\partial\bar{r}}{\partial{u}}\times\frac{\partial\bar{r}}{\partial{v}}\right|dudv$$
<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## The surface area of a rotation paraboloid?
* Note surface area can be obtained by the surface integral of the normal vector field.

The surface $S$ is given by:
$$\bar{r}(u,v) = \begin{cases}x^2+y^2\leq1\\z=x^2+y^2\end{cases}$$

The surface is given by:

$$\int_S dS\Phi(\bar{r}),\; \Phi(\bar{r}) = 1$$
It is easier to work with cylindrical coordinates:
$$\begin{cases} x = \rho\cos\phi\\y = \rho\sin\phi\\z = \rho^2\end{cases}$$
If we call $u = \rho$ and $v = \phi$ we get:

$$\bar{r}(u,v) = \rho\bar{e}_\rho+z\bar{e}_z= \rho\bar{e}_\rho+\rho^2\bar{e}_z$$
Returning to cartesian coordinates we gete:
$$\bar{r}(u,v) = \rho\cos\phi\bar{e}_x+\rho\sin\phi\bar{e}_y+\rho^2\bar{e}_z$$
Finding the surface element by finding the partial derivatives:
$$\frac{\partial\bar{r}}{\partial{\rho}} = \cos\phi\bar{e}_x+\sin\phi\bar{e}_y+2\rho\bar{e}_z= \bar{e}_\rho+2\rho \bar{e}_z$$
$$\frac{\partial\bar{r}}{\partial{\phi}} = -\rho\sin\phi\bar{e}_x+\rho\cos\phi\bar{e}_y$$

Then the absolute surface element is:
$$\left|d\bar{S}\right| = \left|\frac{\partial\bar{r}}{\partial{\rho}}\times\frac{\partial\bar{r}}{\partial{\phi}}\right|= \left|\rho \bar{e}_z-2\rho^2\bar{e}_\rho\right|=\sqrt{4\rho^4+\rho^2}=(\sqrt{4\rho^2+1})\rho$$
The area is then:
$$ \int_S dS = \int_0^{2\pi}d\phi \int_0^1 (\sqrt{4\rho^2+1})\rho d\rho = \frac{\pi}{6}\left({5^{\frac{3}{2}}-1}\right)$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## The force on the paraboloid if it is filled wih a fluid of constant density.

$$p = \rho\cdot{g}\cdot{z}$$
Choose units such that $g\cdot{\rho} = 1$.
Then the pressure is gien by:
$$p = 1-z$$

Let the surface $S$ be the paraboloid:

$$S  = \begin{cases}\rho\leq1\\z=\rho^2\end{cases}$$
The parametization $\bar{r}$ is given by:
$$\bar{r} = \rho\bar{e}_\rho+\rho^2\bar{e}_z$$
Then the surface element $d\bar{S}$ is given by:
$$\frac{\partial\bar{r}}{\partial{\rho}} \times \frac{\partial\bar{r}}{\partial{\phi}} d\rho d\phi= (\rho\bar{e}_z-2\rho^2\bar{e}_\phi) d\rho d\phi$$

Then the force is given by:
$$\bar{F} = -\int_S{p d\bar{S}} = - \int_0^1 d\rho \int_0^{2\pi} \left( \rho\bar{e}_z-2\rho^2\bar{e}_\rho\right )\left(1-\rho^2\right)  d\phi$$
Since $\bar{e}_\rho$ depends on $\phi$ we can integrate over $\phi$ first using 
$$\bar{e}_\rho = \cos\phi\bar{e}_x+\sin\phi\bar{e}_y$$
And by using symmetry we get:
$$\int_0^{2\pi} \bar{e}_\rho d\phi = 0 $$
Since when we integrate over $\phi$ we get:
$$\int_0^{2\pi} \bar{e}_\rho d\phi = \int_0^{2\pi} \cos(\phi)\bar{e}_x + \sin(\phi)\bar{e}_y d\phi = 0\bar{e}_x+0\bar{e}_y= 0$$
Then the force is given by:
Then we get:
$$\bar{F} = -2\pi\bar{e}_z\int_0^1\rho(1-\rho^2)d\rho = -\frac{\pi}{2}\bar{e}_z$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

Compute the out of a sphere of radius $R$ when $\bar{A}$ is given by:
$$\bar{A} = -\left(\frac{1}{r^2}+\frac{\lambda}{r}\right)e^{-\lambda R}\bar{e_r}$$
A visualization of the sphere:

<div style="text-align:center"><img src="LectureImages\Lecture5\Lecture5SphericalFlux.png" /></div>

Through the surface $S$:
$$d\bar{S} = r^2\sin(\theta)d\theta\ d\phi\bar{e_r}$$
$$\frac{\partial\bar{r}}{\partial{\theta}} = r\frac{\partial}{\partial{\theta}}\bar{e_r}= r\bar{e_\theta}$$
$$\frac{\partial\bar{r}}{\partial{\phi}} = r\frac{\partial}{\partial{\phi}}\bar{e_r}= r\sin\theta\bar{e}_\phi$$
Using the cross product to obtain the surface element gives:
$$d\bar{S} = r^2\sin(\theta)d\theta\ d\phi\bar{e_r} = r^2\sin(\theta)\bar{e_\theta}\times r\sin\theta\bar{e}_\phi = r^2\sin(\theta)\bar{e}_\phid\theta\ d\phi$$

Then the flux is given by:
$$F = \int_S{d\bar{S}\cdot\bar{A}} = \int_S{d\bar{S}\cdot\bar{A}} = \int_S{r^2\sin(\theta)d\theta d\phi\bar{e}_\phi\cdot\left(-\left(\frac{1}{r^2}+\frac{\lambda}{r}\right)e^{-\lambda R}\bar{e_r}\right)}$$
But since we are on the surface of the sphere we can use the fact that $r=R$ and that $\bar{e}_r$ is a unit vector in the direction of $\bar{r}$.
Then we get:
$$F = -\int_0^(\pi)sin(\theta)d\theta\int_0^{2\pi}d\phi\left(\frac{1}{R^2}+\frac{\lambda}{R}\right)e^{-\lambda R} = -4\pi\left(\frac{1}{R^2}+\frac{\lambda}{R}\right)e^{-\lambda R}$$

