# Two-phase Comptonization model of soft thermal seed photons for GRB prompt emission - grbjet

## Description

<p> 
This model computes the time-averaged flux over the signal duration  determined by the radiation curvature effect for a single pulse emitted by a relativistic top-hat jet (<a href="https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.5723F/abstract">Farinelli et al., 2021</a>).
</p>

## Parameters in grbjet

<table>
<tr><td>1</td>  <td>thobs</td>     <td>observer viewing angle (degrees)</td></tr>
<tr><td>2</td>  <td>thjet</td>     <td>jet half opening angle (degrees)</td></tr>
<tr><td>3</td>  <td>gamma</td>     <td>jet gamma Lorentz factor</td></tr>
<tr><td>4</td>  <td>r12</td>       <td>jet radius (10^12 cm)</td></tr> 
<tr><td>5</td>  <td>p1</td>        <td>low-energy index of the comoving frame broken powerlaw spectrum</td></tr>
<tr><td>6</td>  <td>p2</td>        <td>high-energy index of the comoving frame broken powerlaw spectrum</td></tr>
<tr><td>7</td>  <td>E0</td>        <td>break energy (keV)</td></tr>
<tr><td>8</td>  <td>delta</td>     <td>smoothness of the transition between the two powerlaws</td></tr>
<tr><td>9</td>  <td>index_pl</td>  <td>energy index of the comoving-frame cutoff powerlaw spectrum</td></tr>
<tr><td>10</td> <td>ecut</td>      <td>cut-off energy (keV)</td></tr>
<tr><td>11</td> <td>ktbb</td>      <td>comoving frame blackbody temperature (keV)</td></tr>
<tr><td>12</td> <td>model</td>     <td>flag to choose the comoving frame emissivity law. 1: broken powerlaw 2: cutoff powerlaw; 3: blackbody</td></tr>
<tr><td>13</td> <td>redshift</td>  <td>source redshift</td></tr>
<tr><td>14</td> <td>norm</td>      <td>comoving frame emissivity normalization (10^20 erg cm^-2 s^-1 keV^-1 sr^-1)</td></tr>
</table>


## Installing the model


First set a folder where source code files are placed and go there. 
We label it as `/model_folder/`

Then from the XSPEC command prompt line, type:

```
XSPEC12> initpackage grbjet grbjet.dat /model_folder/
XSPEC12> lmod grbjet /model_folder/ 
```

If everything ran correctly, you should find the library file in the folder.
It is named `libgrbjet.so` on Linux systems and `libgrbjet.dylib` on OS X systems.

In order to load the model library every time that XSPEC is launched from the 
command line, edit the file `$HOME/.xspec/xspec.rc`  and add the following line:

```
lmod grbjet /model_folder/
```

This can be checked after running XSPEC, when you should see the message

```
"Model package grbjet successfully loaded."
```

Now the model is ready to be used.


