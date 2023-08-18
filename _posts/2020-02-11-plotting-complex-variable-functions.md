---
title: Plotting complex variable functions
date: "2020-02-11"
created: '2020-02-11T21:54:53.008Z'
modified: '2020-04-12T22:09:40.883Z'
author_profile: false
toc: true
toc_sticky: true
toc_label: "Table of Contents"
toc_icon: "fas fa-list-ul"
---

<div style="text-align: justify">
When we try to plot a complex variable function we encounter a first problem: we need a four-dimensional space. This is because this type of function goes from the complex plane to the complex plane, that is, two dimensions are associated to the input variable of the function and another two to the result of the function.

In this post I show a small compilation of some methods that can be used to solve this problem. These methods are:
</div>

  - [Plotting the function's module](#plotting-the-functions-module)
  - [Plotting the real and imaginary part](#plotting-the-real-and-imaginary-part)
  - [Domain coloring](#domain-coloring)
  - [Conformal mapping](#conformal-mapping)
  - [Plotting as a vector field](#plotting-as-a-vector-field)


<div style="text-align: justify">
All the examples will be included with their respective <b>Python 3</b> code.
</div>

## Plotting the function's module

<div style="text-align: justify">
The first method that could be used would be the representation of the module of the complex function. With this, what we are doing is turning our four-dimensional representation into a three-dimensional graph of a surface. In this way we can obtain information about the singularities of our function, since we are given the possibility of detecting them quickly simply by seeing which parts of the graph go to infinity.

For example, by representing the module of the function $f(z) = \dfrac{1}{z^2+4}$, where $z \in \mathbb{C}$, we can detect two singularities, in this case two simple poles one at $z=-2$ and the other at $z=2$.
</div>
<br>


```python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 500
lim = 4

x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y
f = abs(1/(z**2-4))
f[f>1.3] = 1.3 # Cut the function at the poles for decoration purposes.

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection="3d", xlim=(-lim,lim), ylim=(-lim,lim), zlim=(0,1))

ax.plot_surface(x, y, f, cmap="viridis", shade=True, alpha=1)
ax.set_xlabel("$Re(z)$", size=14)
ax.set_ylabel("$Im(z)$", size=14)
ax.set_title("$|f(z)|=|\dfrac{1}{z^2-4}|$", size=18, pad=30)

plt.show()
```

<img src="/assets/images/2020-02-11/output_2_0.png" style="width:50%" class="center">{: .align-center}

<div style="text-align: justify">
Or representing the absolute value of the gamma function, $\Gamma(z)= \int_{0}^{\infty} t^{z-1} e^{-t} dt $, we can see its single poles at the negative integers and at $z=0$.
</div>

<br>


```python
from scipy.special import gamma

f = abs(gamma(z))
f[f>6] = 6

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection="3d", xlim=(-lim,lim), ylim=(-lim,lim), zlim=(0,6))

ax.plot_surface(x, y, f, cmap="viridis", shade=True, alpha=1)
ax.set_xlabel("$Re(z)$", size=14)
ax.set_ylabel("$Im(z)$", size=14)
ax.set_title("$|f(z)|=|\Gamma (z)|$", size=18, pad=30)
ax.view_init(azim=-130, elev=35)

plt.show()
```

<img src="/assets/images/2020-02-11/output_4_0.png" style="width:50%" class="center">{: .align-center}


## Plotting the real and imaginary part

<div style="text-align: justify">
The next method would be to graph separately the real and imaginary part of the function. A similar thing we could do is to graph separately the module and the argument of the function.

As an example, we have the graphs of the complex sine $f(z)= \sin z = \dfrac{e^{iz}-e^{-iz}}{2i}$. We can see how its components diverge as we move away from the origin parallel to the imaginary axis of the complex plane.
</div>
<br>

```python
N = 50
lim = 7

x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y
f = np.sin(z)

fig = plt.figure(figsize=(12,5))
fig.suptitle("$f(z) = \sin(z)$", fontsize=20)
ax1 = fig.add_subplot(121, projection="3d", xlim=(-lim,lim), ylim=(-lim,lim))
ax2 = fig.add_subplot(122, projection="3d", xlim=(-lim,lim), ylim=(-lim,lim))

ax1.plot_surface(x, y, f.real, cmap="viridis", shade=True, alpha=0.9, label="Re f(z)")
ax2.plot_surface(x, y, f.imag, cmap="viridis", shade=True, alpha=0.9, label="Im f(z)")

ax1.set_zlim(f.real.min(), f.real.max())
ax1.set_xlabel("$Re(z)$", fontsize=14)
ax1.set_ylabel("$Im(z)$", fontsize=14)
ax1.set_title("$Re \,\, f(z)$") # \, adds space

ax2.set_zlim(f.imag.min(), f.imag.max())
ax2.set_xlabel("$Re(z)$", fontsize=14)
ax2.set_ylabel("$Im(z)$", fontsize=14)
ax2.set_title("$Im \,\, f(z)$")

plt.show()
```

<img src="/assets/images/2020-02-11/output_6_0.png">{: .align-center}


## Domain coloring

<div style="text-align: justify">
Another technique that can be used in order to represent complex functions would be to color the domain of our function expressed in polar form $f(z)=|f(z)| \, e^{\, i  Arg f(z)} = r e^{i\theta}$, where $r \geq 0$ and $\theta \in (-\pi,\pi]$. The most common way to do domain coloring is by using the Hue, Saturation and Brightness (HLS) system. Then, what we would do is color each point of the domain according to the argument $\theta$ and modulus $r$ of the function.

An example of criteria for choosing the colors would be:
</div>

$$
\begin{cases}
H &= Arg \, f(z) \\[0.5ex]
L &= 1-a^{|f(z)|},\; 0<a<1 \\[0.5ex]
S &= 1
\end{cases}
$$

<div style="text-align: justify">
but, precisely, for the following demonstrations I will be using 
</div>

$$
\begin{cases}
r &\rightarrow \log_2(1+r) \\[0.5ex]
H &= \frac{Arg \, f(z)}{2\pi} \\[0.5ex]
L &= 1-0.45^{\ln\left(1+|f(z)|\right)} \\[0.5ex]
S &= 1
\end{cases}
$$

<div style="text-align: justify">
With this, a higher module value translates into lighter colors, $\theta=-\pi, \pi$ corresponds to the cyan color and $\theta=0$ to the red color.

For example we choose the function $f(z) = (z-1)^5 + 1$. We can see how its five zeros are in the unit radius circle centered at $z=1$.
</div>
<br>

```python
from colorsys import hls_to_rgb


def colorize(fz):

    """
    The original colorize function can be found at:
    https://stackoverflow.com/questions/17044052/mathplotlib-imshow-complex-2d-array
    by the user nadapez.
    """
    
    r = np.log2(1. + np.abs(fz))
    
    h = np.angle(fz)/(2*np.pi)
    l = 1 - 0.45**(np.log(1+r)) 
    s = 1

    c = np.vectorize(hls_to_rgb)(h,l,s) # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (m,n,3)
    c = np.rot90(c.transpose(2,1,0), 1) # Change shape to (m,n,3) and rotate 90 degrees
    
    return c


N = 500
lim = 3
x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y
f = (z-1)**5+1


from matplotlib.lines import Line2D

legend_elements = [Line2D([0], [0], marker='o', color='cyan', label='$Arg$ $f(z) =$ $-\pi$,$\pi$', markersize=10, lw=0),
                   Line2D([0], [0], marker='o', color='red', label='$Arg$ $f(z)=0$', markersize=10, lw=0)]

# Create the figure
fig, ax = plt.subplots()
ax.legend(handles=legend_elements, loc='upper left')

img = colorize(f)
ax.imshow(img, extent=[-lim,lim, -lim,lim])
ax.set_xlabel("$Re$ $Z$", fontsize=14)
ax.set_ylabel("$Im$ $Z$", fontsize=14)
ax.set_title(r"$f(z)=(z-1)^5+1$", fontsize=14)

plt.tight_layout()
plt.show()
```

<img src="/assets/images/2020-02-11/output_8_0.png" style="width:50%" class="center">{: .align-center}


<div style="text-align: justify">
In this way, in just one graph, the variation of the two components of the polar function can be seen on the complex plane, which makes it possible to locate points of interest that the function may contain.

In this example we can clearly see where the zeros of the function would be located.
But, in the following examples, we will also identify the order of the poles and the types of singularities.

To begin with, let's choose the function $f(z) = \dfrac{1}{(z-1)(z+1)^2}$.

We see that in $z = 1$ it has a simple pole, but in $z=-1$ a pole of order 2. This can be detected by looking at the values of the colour spectrum around the singularity. The color spectrum turns twice around $z=-1$, an order two pole, but once around $z=1$, an order one pole. This multiplicity of turns in the colour can also be seen in the zeros of order greater than 1 in other functions (informal note: $e^{i 2θ}$ turns twice as fast as $e^{iθ}$).
</div>
<br>

```python
lim = 2
x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y

f = 1/((z-1)*(z+1)**2)

fig, ax = plt.subplots(figsize=(5,5))
ax.legend(handles=legend_elements, loc='upper left')

img = colorize(f)
ax.imshow(img, extent=[-lim,lim, -lim,lim])
ax.set_xlabel("$Re$ $Z$", fontsize=14)
ax.set_ylabel("$Im$ $Z$", fontsize=14)
ax.set_title(r"$f(z)=\frac{1}{(z-1)(z+1)^2}$", fontsize=14)

plt.tight_layout()
plt.show()
```

<img src="/assets/images/2020-02-11/output_10_0.png" style="width:50%" class="center">{: .align-center}


<div style="text-align: justify">
The following examples are the functions $\displaystyle f(z)=e^{\frac{1}{z}}$ and $f(z)=\frac{\sin z}{z}$. 
</div>
<br>


```python

#__exp(1/z)___________________________________________

lim = 1
x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y

f = np.exp(1/z)

fig, (ax, ax2) = plt.subplots(1,2,figsize=(10,4))
ax.legend(handles=legend_elements, loc='upper left')

img = colorize(f)
ax.imshow(img, extent=[-lim,lim, -lim,lim])
ax.set_xlabel("$Re$ $Z$", fontsize=14)
ax.set_ylabel("$Im$ $Z$", fontsize=14)
ax.set_title(r"$f(z)=e^{1/z}$", fontsize=14)


#__sin(z)/z___________________________________________

lim = 10
x, y = np.meshgrid(np.linspace(-lim,lim,N),
                   np.linspace(-lim,lim,N))
z = x + 1j*y

f = np.sin(z)/z

ax2.legend(handles=legend_elements, loc='upper left')

img = colorize(f)
ax2.imshow(img, extent=[-lim,lim, -lim,lim])
ax2.set_xlabel("$Re$ $Z$", fontsize=14)
ax2.set_ylabel("$Im$ $Z$", fontsize=14)
ax2.set_title(r"$f(z)=\frac{\sin(z)}{z}$", fontsize=14)

plt.tight_layout()
plt.show()
```


<img src="/assets/images/2020-02-11/output_12_0.png">{: .align-center}


<div style="text-align: justify">
$e^{\frac{1}{z}}$ is an example of function with an essential singularity. We can see how at $z=0$ the color spectrum takes infinitely many turns and the module becomes zero when approaching $z=0$ from the left but becomes infinite when approaching it from the right.
The function ends up taking all the values of $\mathbb{C}$ as verified by the Casorati-Weierstrass theorem (see references).

On the other hand, we can see the zeros of the function $\frac{sin(z)}{z}$ in black at the multiples of $\pi$ and we can see that the supposed singularity at $z=0$ does not appear, it is avoided. $\frac{sin(z)}{z}$ has an avoidable singularity at $z=0$.
</div>

## Conformal mapping

<div style="text-align: justify">
Another way to represent a complex function is by drawing how the curves and points on the complex plane are transformed when the function is applied. This method would be the transformation graph or conformal mapping <i>(remark: a conformal map is a function that preserves the angles)</i>.

A common way to do this is by taking the horizontal and vertical lines in the complex plane. 

If we have a function f(z) we can represent what the contour lines would be by equating its real and imaginary components to real constants $C_1$ and $C_2$ respectively. 
</div>
$$f(z) = \mbox{Re} f(z) + i \, \mbox{Im} f(z)$$

$$
\begin{align}
\mbox{Re} f(z) & = C_1 \\
\mbox{Im} f(z) & = C_2
\end{align}
$$

<div style="text-align: justify">
We would be representing the contour lines as if we were in the real plane $\mathbb{R}^2$. For example we take the function $f(z)=z+\frac{1}{z}$, with $z = x+iy$ where $x=\mbox{Re} \, z$ and $y=\mbox{Im} \, z$ and $\overline{z} = x-iy $ is the complex conjugate of $z$:
</div>
<br>
<div style="text-align: center">
$$
\begin{split}
f(z) & =z+\frac{1}{z} = \\[2ex]
& =  z+\frac{\overline{z}}{|z|^2} = \\[2ex]
& =  x + iy + \frac{x}{x^2+y^2} - i \frac{y}{x^2+y^2} = \\[2ex]
& = \left( x + \frac{x}{x^2+y^2} \right) + i \left( y - \frac{y}{x^2+y^2} \right) = \\[2ex]
& = \mbox{Re} f(z) + i \, \mbox{Im} f(z)
\end{split}
$$
</div>
<br>

```python
def f(z):
    return z + 1/z


# The x and y coordinates
lim = 3
N = 300
Xv = np.linspace(-lim, lim, N)
Yv = np.linspace(-lim, lim, N)
X, Y = np.meshgrid(Xv, Yv)

# Values of f as a function of z=x+iy
Z = X+1j*Y
Fv = f(Z)
lv = np.linspace(-10,10,31)

# Contours of constant Re f(z) and Im f(z) as a function of x and y

fig, (ax, ax2) = plt.subplots(1,2,figsize=(10,5))
fig.suptitle(r"$f(z)=z+\frac{1}{z}$", fontsize=18)
ax.set_aspect("equal")
ax2.set_aspect("equal")

ax.contour(Xv, Yv, X, colors='blue', linestyles='solid', levels=lv, linewidths=1)
ax.contour(Xv, Yv, Y, colors='red', linestyles='solid', levels=lv, linewidths=1)
ax.set_xlabel("$Re$ $Z$", fontsize=14)
ax.set_ylabel("$Im$ $Z$", fontsize=14)
ax.set_title("Before")

ax2.contour(Xv, Yv, np.real(Fv), colors='blue', linestyles='solid', levels=lv, linewidths=1)
ax2.contour(Xv, Yv, np.imag(Fv), colors='red', linestyles='solid', levels=lv, linewidths=1)
ax2.set_xlabel("$Re$ $Z$", fontsize=14)
ax2.set_ylabel("$Im$ $Z$", fontsize=14)
ax2.set_title("After")

plt.subplots_adjust(wspace=0.5)
plt.show()
```


<img src="/assets/images/2020-02-11/output_15_0.png">{: .align-center}


<div style="text-align: justify">
This method is used, for example, in fluid dynamics to represent the flow of a fluid around objects. As a concrete example, I show how the function z+1/z could describe the flow around the section of a cylindrical body of radius the unity. This transform is called a Joukowski transform (see references, N.Hall). 
</div>
<br>

```python
from matplotlib.patches import Circle

fig, ax = plt.subplots()
ax.set_aspect("equal")

Circ = Circle((0,0), radius=1, facecolor="black", alpha=1, zorder=10)
ax.add_patch(Circ)

ax.contour(Xv, Yv, np.imag(Fv), colors='red', linestyles='solid', levels=lv, linewidths=1)
ax.set_xlabel("Re $Z$", fontsize=14)
ax.set_ylabel("Im $Z$", fontsize=14)

plt.tight_layout()
plt.show()
```


<img src="/assets/images/2020-02-11/output_17_0.png" style="width:50%" class="center">{: .align-center}


<div style="text-align: justify">
The procedure of conformal mapping, apart from fluid mechanics, is also applied in cartography, general relativity, electrostatics, scattering and diffraction problems, medical physics for brain surface mapping, etc (see S.Ganguli).
</div>

## Plotting as a vector field

<div style="text-align: justify">
Another form of representation complementary to the one we have just seen, also used in the areas mentioned before, would be the representation of the vector field generated by the complex function.

In order to represent the complex function as a vector field, what is done is to separate the real part of the function from the imaginary one. Then, if we take any point in space , say $z=x+iy$, the two components of our vector at that point would be the value of the real component of the function at that point and the value of the imaginary component respectively: $\vec{v} = \left(\mbox{Re} f(z), \mbox{Im} f(z) \right)$.
</div>
<br>

```python
def f(z):
    return np.cos(z)


# The x and y coordinates
lim = 3
N = 30
Xv = np.linspace(-lim-1, lim+1, N)
Yv = np.linspace(-lim, lim, N)
X, Y = np.meshgrid(Xv, Yv)

# Values of f as a function of z=x+iy
Z = X+1j*Y
Fv = f(Z)
lv = np.linspace(-10,10,31)

# Plotting part
fig, ax = plt.subplots(figsize=(6,5))
ax.set_aspect("equal")
ax.set_xlabel("Re $Z$", fontsize=14)
ax.set_ylabel("Im $Z$", fontsize=14)
ax.set_title(r"$f(z)=\cos \, z$", fontsize=14)
ax.quiver(X, Y, np.real(Fv), np.imag(Fv), color='blue', pivot="middle", norm=True, headwidth=6, headlength=7 )


plt.tight_layout()
plt.show()
```


<img src="/assets/images/2020-02-11/output_20_0.png" style="width:70%" class="center">{: .align-center}


<div style="text-align: justify">
In this figure we can see that the real part of $\cos z$ has zeros at the points where $\mbox{Re} z = \pm \frac{\pi}{2} \approx \pm 1.57 $. We can say so because at those points the vectors $\vec{v} = \left(\mbox{Re} f(z), \mbox{Im} f(z) \right)$ only have vertical (or imaginary) component. In addition, we can determine if the function is positive or negative at those points just by looking if the vectors are positive or negative in the vertical direction. The modulus of the vectors can give information about the modulus of the function. 
</div>

$$|\vec{v}| = \sqrt{(\mbox{Re} f(z))^2 + (\mbox{Im} f(z))^2} = |f(z)|$$

<br>

***

<br>

# References

- C. Fernández González (n.d.). _¿Cómo representar funciones en variable compleja?_ [online]. Portal.uned.es. Available at: [http://portal.uned.es/pls/portal/docs/PAGE/UNED_MAIN/LAUNIVERSIDAD/<br/>
UBICACIONES/01/DOCENTE/CARLOS_FERNANDEZ_GONZALEZ/<br/>
VARIABLE%20COMPLEJA%20GEOGEBRA/REPRESENTARVC.PDF](http://portal.uned.es/pls/portal/docs/PAGE/UNED_MAIN/LAUNIVERSIDAD/UBICACIONES/01/DOCENTE/CARLOS_FERNANDEZ_GONZALEZ/VARIABLE%20COMPLEJA%20GEOGEBRA/REPRESENTARVC.PDF) 
<br>

- C. Maggi, et al. (n.d.). _Aplicaciones gráficas de funciones complejas_.
[online] academia.edu. Available at: [https://www.academia.edu/2211067/Aplicaciones_gr%C3%A1ficas_de_funciones_complejas](https://www.academia.edu/2211067/Aplicaciones_gr%C3%A1ficas_de_funciones_complejas)
<br>

- G. Viswanathan (2014). _Domain coloring for visualizing complex
functions_. [online] Gandhi Viswanathan’s Blog. Available at: [https://gandhiviswanathan.wordpress.com/2014/10/07/domain-coloring-for-visualizing-complex-functions/](https://gandhiviswanathan.wordpress.com/2014/10/07/domain-coloring-for-visualizing-complex-functions/)
<br>

- En.wikipedia.org. (2019). _Domain coloring_. [online] Available at: [https://en.wikipedia.org/wiki/Domain_coloring](https://en.wikipedia.org/wiki/Domain_coloring)
<br>

- Wikidot.com. (n.d.). _The Casorati-Weierstrass theorem_. [online] Available at: [http://mathonline.wikidot.com/the-casorati-weierstrass-theorem](http://mathonline.wikidot.com/the-casorati-weierstrass-theorem)
<br>

- N. Hall. (2018). _Conformal mapping, Joukowsky transformation_. [online] grc.nasa.gov. Available at: [https://www.grc.nasa.gov/www/k-12/airplane/map.html](https://www.grc.nasa.gov/www/k-12/airplane/map.html)
<br>

- S. Ganguli. (2008). _Conformal Mapping and its Applications_. [online] Iiserpune.ac.in. Available at: [http://www.iiserpune.ac.in/~p.subramanian/conformal_mapping1.pdf](http://www.iiserpune.ac.in/~p.subramanian/conformal_mapping1.pdf)
<br>
