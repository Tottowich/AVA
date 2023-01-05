<!-- FormulaBox
<div class="FormulaBox" style="border: 2px solid Magenta; padding: 10px 5px"> -->
<!-- Set the margin for the document -->
<div class="Document" style="margin: 10px 85px">

# Lecture 3

* **Gradient** $\Phi(x,y,z)$ is a continously differentiable function at the point $\bar{r}=(x,y,z)$
$$\bar{r}=x\bar{e}_x+y\bar{e}_y+z\bar{e}_z$$
Moving an infinitesimal distance results in:
$$\bar{r}'=(x+\Delta x)\bar{e}_x+(y+\Delta y)\bar{e}_y+(z+\Delta z)\bar{e}_z$$
The difference in the function value is:
$$\Delta\bar{r}=\bar{r}'-\bar{r}=\Delta x\bar{e}_x+\Delta y\bar{e}_y+\Delta z\bar{e}_z$$
As $\Delta x,\Delta y,\Delta z$ are infinitesimal, we can rewrite:
$$\Delta\bar{r}=\bar{r}'-\bar{r}=d x\bar{e}_x+d y\bar{e}_y+d z\bar{e}_z$$
The change in the function value is:
$$\Delta\Phi=\Phi(x+\Delta x,y+\Delta y,z+\Delta z)-\Phi(x,y,z\overset{\Delta x\rightarrow0}{\overset{\Delta y\rightarrow0}{\overset{\Delta z\rightarrow0}{\rightarrow}}} d\Phi(x,y,z)$$
Which is the same as:
$$d\Phi=\frac{\partial\Phi}{\partial x}dx+\frac{\partial\Phi}{\partial y}dy+\frac{\partial\Phi}{\partial z}dz$$
Which can be written with a dot product as:
$$d\Phi=\left(\frac{\partial \Phi}{\partial x}\bar{e}_x+\frac{\partial \Phi}{\partial y}\bar{e}_y+\frac{\partial \Phi}{\partial z}\bar{e}_z\right)\cdot \left (dx\bar{e}_x+dy\bar{e}_y+dz\bar{e}_z\right)=\nabla\Phi\cdot d\bar{r}=grad\:\Phi\cdot d\bar{r}$$
* **Directional Derivative**
Write $d\bar{r}$ as:
$$d\bar{r}=ds\bar{e}$$
Where $\bar{e}$ is a unit vector in the direction of $d\bar{r}$, i.e. 
$$\bar{e}_s=\frac{d\bar{r}}{|d\bar{r}|}$$
The directional derivative is:
$$\frac{d\Phi}{ds}=\nabla\Phi\cdot\bar{e}$$
This tells us how much the function $\Phi$ changes in the direction of $\bar{e}$.
<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Magenta; padding: 10px 5px">

## **Theorem**:
* The vector $\nabla\Phi$ is the gradient of $\Phi$ and is the direction of maximum rate of change of $\Phi$.
## **Proof**
Let $\bar{e}$ be a unit vector in the direction of most rapid change of $\Phi$.
$$\frac{d\Phi}{ds}=\nabla\Phi\cdot\bar{e}=|\nabla\Phi|\cos\theta$$
Where $\theta$ is the angle between $\nabla\Phi$ and $\bar{e}$.
The derivative is maximum when $\cos\theta=1$, since $|\cos\theta|\leq1$ when $0\leq\theta\leq\pi$.
Therefore, $\nabla\Phi$ is the direction of maximum rate of change of $\Phi$.

<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>


</div>
<br>

## Level Surfaces:

A Level surface is a surface on which a scalar field $\Phi$ is constant.
$$\Phi(\bar{r}_c)=C$$
For all points $\bar{r}_c$ on the surface.
One surface for each value of C in the scalar field.

<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Magenta; padding: 10px 5px">

## **Theorem**: 
If $\Phi$ is a scalar field and has a max,min or saddle point at $\bar{r}$, then the gradient of $\Phi$ is $0$ at $\bar{r}$.
$$\nabla\Phi(\bar{r})=\bar{0}$$

</div>
<br>
<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Magenta; padding: 10px 5px">


## **Theorem**: 
If $\Phi$ is a scalar field, then the gradient of $\Phi$ is perpendicular to the level surfaces of $\Phi$.

## **Proof**:
Choosing two infinitesimal vectors $d\bar{r}$ and $d\bar{r}'$ such that then:
$$\Phi(\bar{r_p}+d\bar{r})=\Phi(\bar{r_p})\iff\Phi(\bar{r_p}+d\bar{r}')-\Phi(\bar{r_p})=0$$
But as showed earlier:
$$\Phi(\bar{r_p}+d\bar{r}')-\Phi(\bar{r_p})=d\Phi$$
$$d\Phi=\nabla\Phi\cdot d\bar{r}'=0$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Magenta; padding: 10px 5px">

## **Theorem**:
The distance $\Delta s$ between two iso-surfaces:
$$\Phi(\bar{p})=C,\;\Phi(\bar{p})=C+h$$
Can be approximated by:
$$\Delta s\approx \frac{h}{|\nabla\Phi(\bar{p})|}$$

## **Proof**:
Let
$$\Phi=C,\;\Phi=C+h$$
Be the two iso-surfaces. Choose a vector $d\bar{r}$ such that $d\bar{r} \perp$ to the level surfaces.
Then: 
$$h = \Delta\Phi \approx d\Phi = \nabla\Phi\cdot d\bar{r}=|\nabla\Phi|\Delta s\implies \Delta s\approx \frac{h}{|\nabla\Phi|}$$
Since $d\bar{r}$ is perpendicular to the level surfaces, $|\nabla\Phi|$ is the distance between the two level surfaces.
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>
<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## **Example**:
Find the normal to the surface $z=x^2+y^2$ at the point $(1,2,5)$.
This surface can be seen as a level surface $\Phi=0$ for the scalar field $\Phi=x^2+y^2-z$.
Then:
$$\nabla\Phi=\left(\frac{\partial \Phi}{\partial x}\bar{e}_x+\frac{\partial \Phi}{\partial y}\bar{e}_y+\frac{\partial \Phi}{\partial z}\bar{e}_z\right)=\left(2x\bar{e}_x+2y\bar{e}_y-1\bar{e}_z\right)$$
Evaluating at $(1,2,5)$:
$$\bar{n}=\left(2\bar{e}_x+4\bar{e}_y-1\bar{e}_z\right)$$
Normalizing:
$$\bar{n}=\frac{1}{\sqrt{21}}\left(2\bar{e}_x+4\bar{e}_y-1\bar{e}_z\right)$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

## Scalar Potential:
* If for a given vector filed $\bar{A}$ there is a scalar field $\Phi$ such that:
$$\bar{A}=\nabla\Phi$$
Then $\Phi$ is called the scalar potential of $\bar{A}$.
* If $\bar{A}$ has a potential $\Phi$, then $\Phi+C$ is also a potential for $\bar{A}$.
* The potential is often defined with a minus sign e.g. the electrostatic potential:
$$\bar{E}=-\nabla V$$

<!-- Theorem Box -->
<div class="TheoremBox" style="border: 2px solid Magenta; padding: 10px 5px">

## Conditions for existence of a scalar potential:
* A continously differentiable vector field $\bar{A}$ has a scalar potential if and only if $\bar{A}$ satisfies:
$$\begin{equation}
\frac{\partial A_x}{\partial y}=\frac{\partial A_y}{\partial x}
\end{equation}$$
$$\begin{equation}
\frac{\partial A_y}{\partial z}=\frac{\partial A_z}{\partial y}
\end{equation}$$
$$\begin{equation}
\frac{\partial A_z}{\partial x}=\frac{\partial A_x}{\partial z}
\end{equation}$$

## **Proof**:
Assume that $\bar{A} = \nabla\Phi$.
Then:
$$ A_x=\frac{\partial \Phi}{\partial x}$$
$$ A_y=\frac{\partial \Phi}{\partial y}$$
Using this we known that since $\bar{A}$ is continously differentiable:
$$\frac{\partial A_x}{\partial y}=\frac{\partial}{\partial y}\frac{\partial \Phi}{\partial x}=\frac{\partial}{\partial x}\frac{\partial \Phi}{\partial y}=\frac{\partial A_y}{\partial x}$$
The mixed partial derivatives of $\Phi$ are the same as the mixed partial derivatives of $\bar{A}$.

<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## **Example**: Check if there exists a scalar potential
Check if the vector field $\bar{A}=y\bar{e}_x-x\bar{e}_y$ has a scalar potential.
$$\frac{\partial A_x}{\partial y}=\frac{\partial}{\partial y}y=1\neq \frac{\partial A_y}{\partial x}=\frac{\partial}{\partial x}(-x)=-1$$
The vector field does not have a scalar potential since it does not satisfy equation 1.
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

<!-- Example Box -->
<div class="ExampleBox" style="border: 2px solid Green; padding: 10px 5px">

## **Example**: Find the scalar potential
Find the scalar potential of the vector field $\bar{A}=2y^2\bar{e}_x+(4xy+y^2z^2)\bar{e}_y+(\frac{2}{3}y^3z+z)\bar{e}_z$.
* Known that $\bar{A}=\nabla\Phi$.
We known that each term of $\bar{A}$ corresponds to the partial derivative of $\Phi$ with respect to $x$, $y$ or $z$.
$$A_x = \frac{\partial \Phi}{\partial x} = 2y^2\overset{\int dx}{\implies} \Phi=2y^2x+F(y,z)$$
Using this we can find the other terms succequently:
$$A_y = \frac{\partial \Phi}{\partial y} = 4xy+\frac{\partial}{\partial y}F(y,z)=4xy+y^2z^2\implies$$
Canceling $4xy$ and integrating with respect to $y$ yields:
$$F(y,z)=\frac{1}{3}y^3z^2+G(z)$$
Replacing $F(y,z)$ in the equation for $\Phi$:
$$\Phi= 2y^2x+\frac{1}{3}y^3z^2+G(z)$$
Deriving with respect to $z$:
$$A_z = \frac{\partial \Phi}{\partial z} = \frac{2}{3}y^3z+\frac{\partial}{\partial z}G(z)=\frac{2}{3}y^3z+z\implies G(z)=\frac{1}{2}z^2+C$$
Replacing $G(z)$ in the equation for $\Phi$ concludes the calculation.
$$\Phi= 2y^2x+\frac{1}{3}y^3z^2+\frac{1}{2}z^2+C$$
<!-- Check -->
<div style="text-align:right"> <span style="font-size: 2.5em;margin-right: 20px;">✔</span></div>
</div>
<br>

</div>


