# grbjet

## Parameters in grbjet

<table>
<tr><td>1</td>  <td>thobs</td>     <td> </td></tr>
<tr><td>2</td>  <td>thjet</td>     <td> </td></tr>
<tr><td>3</td>  <td>gamma</td>     <td> </td></tr>
<tr><td>4</td>  <td>r12</td>       <td> </td></tr>
<tr><td>5</td>  <td>p1</td>        <td> </td></tr>
<tr><td>6</td>  <td>p2</td>        <td> </td></tr>
<tr><td>7</td>  <td>E0</td>        <td> </td></tr>
<tr><td>8</td>  <td>delta</td>     <td> </td></tr>
<tr><td>9</td>  <td>index_pl</td>  <td> </td></tr>
<tr><td>10</td> <td>ecut</td>      <td> </td></tr>
<tr><td>11</td> <td>ktbb</td>      <td> </td></tr>
<tr><td>12</td> <td>model</td>     <td>switch </td></tr>
<tr><td>13</td> <td>redshift</td>  <td> </td></tr>
</table>


## Installing the model


First set a folder where source code files are placed and go there. 
We label it as /model_folder/

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



initpackage grbjet grbjet.dat /Users/klrutkow/Documents/Xspec/localmodels/grbjet
lmod grbjet /Users/klrutkow/Documents/Xspec/localmodels/grbjet
