# Multiplicative model zxipab: power-law distribution of ionised absorber

<p>This is an extension of the pwab model in which the complex absorption is modelled as a power-law distribution of covering fraction and column of ionised absorbers. Unlike pwab, which is based on neutral absorbers, zxipab utilizes pre-calculated grid of XSTAR photo-ionisation model also used by zxipcf model, with ionising parameter log (xi). See Islam & Mukai (2021) for details.</p>


## Parameters in zxipab

<table>
  <tr>
    <td>1</td>
    <td>nHmin</td>
    <td>minimum equivalent hydrogen column (10^{22} cm^{-2})</td>
  </tr>
  <tr>
    <td>2</td>
    <td>nHmax</td>
    <td>maximum equivalent hydrogen column (10^{22} cm^{-2})</td>
  </tr>
  <tr>
    <td>3</td>
    <td>beta</td>
    <td>power law index for covering fraction</td>
  </tr>
  <tr>
    <td>4</td>
    <td>logxi</td>
    <td>ionising parameter log (xi)</td>
  </tr>
</table>
