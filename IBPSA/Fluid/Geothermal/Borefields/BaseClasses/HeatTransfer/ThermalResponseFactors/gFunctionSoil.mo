within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors;
function gFunctionSoil "Evaluate the g-function of a bore field"
  extends Modelica.Icons.Function;

  input Modelica.SIunits.Height hBor "Borehole length";
  input Modelica.SIunits.Height dBor "Borehole buried depth";
  input Modelica.SIunits.Radius rBor "Borehole radius";
  input Modelica.SIunits.Distance r
    "Distance from the borehole wall at which the G-function is evaluated";
  input Modelica.SIunits.ThermalDiffusivity aSoi "Ground thermal diffusivity used in g-function evaluation";
  input Integer nTimSho "Number of time steps in short time region";
  input Integer nTimLon "Number of time steps in long time region";
  input Real ttsMax "Maximum adimensional time for gfunc calculation";

  output Modelica.SIunits.Time tGFun[nTimSho+nTimLon] "Time of g-function evaluation";
  output Real g[nTimSho+nTimLon] "g-function";

protected
  Modelica.SIunits.Time ts = hBor^2/(9*aSoi) "Characteristic time";
  Modelica.SIunits.Time tSho_min = 1 "Minimum time for short time calculations";
  Modelica.SIunits.Time tLon_max = ts*ttsMax "Maximum time for long time calculations";
  Modelica.SIunits.Time tVec[nTimSho+nTimLon-1] "Time vector for short time calculations";

  Real FLS "Finite line source solution";
  Real ILS "Infinite line source solution";
  Real CHS "Cylindrical heat source solution";

algorithm
  // Generate geometrically expanding time vectors
  tVec :=
    IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.timeGeometric(
    tSho_min,
    tLon_max,
    nTimSho + nTimLon-1);
  // Concatenate the short- and long-term parts
  tGFun :=  cat(1,{0},tVec);

  // -----------------------
  // Short time calculations
  // -----------------------

  g[1] := 0.;
  for k in 2:nTimSho+nTimLon loop
    // Finite line source solution
    FLS :=
      IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.finiteLineSource(
      tGFun[k],
      aSoi,
      r,
      hBor,
      dBor,
      hBor,
      dBor);
    // Infinite line source solution
    ILS := 0.5*
      IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.infiniteLineSource(
      tGFun[k],
      aSoi,
      r);
    // Cylindrical heat source solution
    CHS := 2*Modelica.Constants.pi*
      IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.cylindricalHeatSource(
      tGFun[k],
      aSoi,
      r,
      rBor);
    // Correct finite line source solution for cylindrical geometry
    g[k] := FLS + (CHS - ILS);
  end for;

annotation (
Documentation(info="<html>
<p>
This function implements the <i>g</i>-function evaluation method introduced by
Cimmino and Bernier (see: Cimmino and Bernier (2014), and Cimmino (2018)) based
on the <i>g</i>-function function concept first introduced by Eskilson (1987).
The <i>g</i>-function gives the relation between the variation of the borehole
wall temperature at a time <i>t</i> and the heat extraction and injection rates
at all times preceding time <i>t</i> as
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/GFunction_01.png\" />
</p>
<p>
where <i>T<sub>b</sub></i> is the borehole wall temperature,
<i>T<sub>g</sub></i> is the undisturbed ground temperature, <i>Q</i> is the
heat injection rate into the ground through the borehole wall per unit borehole
length, <i>k<sub>s</sub></i> is the soil thermal conductivity and <i>g</i> is
the <i>g</i>-function.
</p>
<p>
The <i>g</i>-function is constructed from the combination of the combination of
the finite line source (FLS) solution (see
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.finiteLineSource\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.finiteLineSource</a>),
the cylindrical heat source (CHS) solution (see
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.cylindricalHeatSource\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.cylindricalHeatSource</a>),
and the infinite line source (ILS) solution (see
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.infiniteLineSource\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.infiniteLineSource</a>).
To obtain the <i>g</i>-function of a bore field, each borehole is divided into a
series of <code>nSeg</code> segments of equal length, each modeled as a line
source of finite length. The finite line source solution is superimposed in
space to obtain a system of equations that gives the relation between the heat
injection rate at each of the segments and the borehole wall temperature at each
of the segments. The system is solved to obtain the uniform borehole wall
temperature required at any time to maintain a constant total heat injection
rate (<i>Q<sub>tot</sub> = 2&pi;k<sub>s</sub>H<sub>tot</sub>)</i> into the bore
field. The uniform borehole wall temperature is then equal to the finite line
source based <i>g</i>-function.
</p>
<p>
Since this <i>g</i>-function is based on line sources of heat, rather than
cylinders, the <i>g</i>-function is corrected to consider the cylindrical
geometry. The correction factor is then the difference between the cylindrical
heat source solution and the infinite line source solution, as proposed by
Li et al. (2014) as
</p>
<p align=\"center\">
<i>g(t) = g<sub>FLS</sub> + (g<sub>CHS</sub> - g<sub>ILS</sub>)</i>
</p>
<h4>Implementation</h4>
<p>
The calculation of the <i>g</i>-function is separated into two regions: the
short-time region and the long-time region. In the short-time region,
corresponding to times <i>t</i> &lt; 1 hour, heat interaction between boreholes
and axial variations of heat injection rate are not considered. The
<i>g</i>-function is calculated using only one borehole and one segment. In the
long-time region, corresponding to times <i>t</i> &gt; 1 hour, all boreholes
are represented as series of <code>nSeg</code> line segments and the
<i>g</i>-function is evaluated as described above.
</p>
<h4>References</h4>
<p>
Cimmino, M. and Bernier, M. 2014. <i>A semi-analytical method to generate
g-functions for geothermal bore fields</i>. International Journal of Heat and
Mass Transfer 70: 641-650.
</p>
<p>
Cimmino, M. 2018. <i>Fast calculation of the g-functions of geothermal borehole
fields using similarities in the evaluation of the finite line source
solution</i>. Journal of Building Performance Simulation. DOI:
10.1080/19401493.2017.1423390.
</p>
<p>
Eskilson, P. 1987. <i>Thermal analysis of heat extraction boreholes</i>. Ph.D.
Thesis. Department of Mathematical Physics. University of Lund. Sweden.
</p>
<p>
Li, M., Li, P., Chan, V. and Lai, A.C.K. 2014. <i>Full-scale temperature
response function (G-function) for heat transfer by borehole heat exchangers
(GHEs) from sub-hour to decades</i>. Applied Energy 136: 197-205.
</p>
</html>", revisions="<html>
<ul>
<li>
October 18, 2018 by Iago Cupeiro:<br/>
Adapted model for evaluation at different distances.
</li>
<li>
March 22, 2018 by Massimo Cimmino:<br/>
First implementation.
</li>
</ul>
</html>"));
end gFunctionSoil;
