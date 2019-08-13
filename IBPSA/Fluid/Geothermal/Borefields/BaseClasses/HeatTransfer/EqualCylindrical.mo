within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer;
model EqualCylindrical
  "Heat conduction in a cylinder using equally radial discretization"

  parameter IBPSA.Fluid.Geothermal.Borefields.Data.Soil.Template soiDat
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));
  parameter Modelica.SIunits.Height h "Height of the cylinder";
  parameter Modelica.SIunits.Radius r_a "Internal radius";
  parameter Modelica.SIunits.Radius r_b "External radius";
  parameter Integer nSta(min=1) = 10 "Number of state variables";
  parameter Modelica.SIunits.Temperature TInt_start
    "Initial temperature at port_a, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Modelica.SIunits.Temperature TExt_start
    "Initial temperature at port_b, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Boolean steadyStateInitial=false
    "true initializes dT(0)/dt=0, false initializes T(0) at fixed temperature using T_a_start and T_b_start"
    annotation (Dialog(group="Initialization"), Evaluate=true);

  parameter Real gridFac(min=1) = 2 "Grid factor for spacing";

  parameter Modelica.SIunits.Radius r[nSta + 1](each fixed=false)
   "Radius to the boundary of the i-th domain";

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a(T(start=TInt_start))
    "Heat port at surface a" annotation (Placement(transformation(extent={{-110,
            -10},{-90,10}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b(T(start=TExt_start))
    "Heat port at surface b" annotation (Placement(transformation(extent={{90,-10},
            {110,10}},rotation=0)));

  Modelica.SIunits.Temperature T[nSta](
    start={TInt_start +
      (TExt_start - TInt_start)/Modelica.Math.log(r_b/r_a)*
      Modelica.Math.log((r_a + (r_b - r_a)/(nSta)*(i - 0.5))/r_a) for i in 1:nSta})
    "Temperature of the states";

  Modelica.SIunits.TemperatureDifference dT "port_a.T - port_b.T";

  Modelica.SIunits.HeatFlowRate Q_flow[nSta + 1]
    "Heat flow rate from state i to i+1";

protected
  parameter Modelica.SIunits.Radius rC[nSta](each fixed=false)
    "Radius to the center of the i-th domain";

  final parameter Modelica.SIunits.SpecificHeatCapacity c=soiDat.cSoi
    "Specific heat capacity";
  final parameter Modelica.SIunits.ThermalConductivity k=soiDat.kSoi
    "Thermal conductivity of the material";
  final parameter Modelica.SIunits.Density d=soiDat.dSoi
    "Density of the material";

  parameter Modelica.SIunits.ThermalConductance G[nSta + 1](each fixed=false)
    "Heat conductance between the temperature nodes";
  parameter Modelica.SIunits.HeatCapacity C[nSta](each fixed=false)
    "Heat capacity of each state";

  parameter Real gridFac_sum(fixed=false);
  parameter Real gridFac_sum_old(fixed=false);

initial algorithm
  for i in 0:nSta - 3 - 1 loop
    if i == 0 then
      gridFac_sum := gridFac^i;
      gridFac_sum_old := gridFac_sum;
    else
      gridFac_sum := gridFac_sum_old + gridFac^i;
      gridFac_sum_old := gridFac_sum;
    end if;
  end for;

initial equation
  assert(r_a < r_b, "Error: Model requires r_a < r_b.");
  assert(0 < r_a, "Error: Model requires 0 < r_a.");

  for i in 1:nSta + 1 loop
    r[i] = r_a + (r_b - r_a)/(nSta)*(i-1);
  end for;

  assert(abs(r[nSta + 1] - r_b) < 1E-3,
    "Error: Wrong computation of radius. r[nSta+1]=" + String(r[nSta + 1]));

  // Radius at middle of resistance
  for i in 1:nSta loop
    rC[i] = (r[i] + r[i + 1])/2;
  end for;

// Conductance between nodes (which are in the center of the domain)
//      G[1] = 2*Modelica.Constants.pi*k*h/Modelica.Math.log(rC[1]/r_a);
//      G[nSta + 1] = 2*Modelica.Constants.pi*k*h/Modelica.Math.log(r_b/rC[nSta]);
//      for i in 2:nSta loop
//        G[i] = 2*Modelica.Constants.pi*k*h/Modelica.Math.log(rC[i]/rC[i - 1]);
//      end for;

     G = {
507.5042165205714,
19662.136771088783,
19132.327557531484,
18561.176087961536,
17988.262463979874,
17436.161804781827,
16916.072477254864,
16457.164555152893,
16036.637217269954,
15637.279108753712,
15268.919965386684,
14930.443377861999,
14643.989974389779,
14377.897065240086,
14124.93526360193,
13890.148479918877,
13690.959541830836,
13512.688913937176,
13339.941751555785,
13182.298141024874,
13050.198038387733,
12929.219989331192,
12813.591250358555,
12716.289349018345,
12631.683480958702,
12551.351650062124,
12484.829232399174,
12430.370006117695,
12378.378682730918,
12338.076465234557,
12307.101736043423,
12279.804677860748,
12262.204702785551,
12250.599304472833,
12244.527527150649,
12245.682585002702,
12251.753167435387,
12263.47399540132,
12280.626934378804,
12302.860723809283,
12329.721401652743,
12361.382144743926,
12397.500829947787,
12437.990867554387,
12483.032739577815,
12531.403930791346,
12584.738015402496,
12641.455628066822,
12702.159103661417,
12766.813584628433,
12834.385827702534,
12906.403089906677,
12980.737847212537,
13059.25385886371,
13140.743241718252,
13225.03751354971,
13313.191470407153,
13403.510502050307,
13497.75585806065,
13594.32322730159,
13693.708800553935,
13796.733306735952,
13901.64859278736,
14010.654703645934,
14120.9613342,
14234.899768946092,
14351.71474572582,
14470.595977741204,
14593.093797728765,
14716.93381066978,
14844.70621022523,
14974.670854378845,
30139.897374826243};

  // Heat capacity of each segment
  for i in 1:nSta loop
    C[i] = 100*(d*Modelica.Constants.pi*c*h*((r[i + 1])^2 - (r[i])^2));
  end for;
  // The initialization is only done for materials that store energy.
  if not soiDat.steadyState then
    if steadyStateInitial then
      der(T) = zeros(nSta);
    else
      for i in 1:nSta loop
        T[i] = TInt_start + (TExt_start - TInt_start)/Modelica.Math.log(r_b/r_a)
          *Modelica.Math.log(rC[i]/r_a);
      end for;
    end if;
  end if;
equation
  dT = port_a.T - port_b.T;
  port_a.Q_flow = +Q_flow[1];
  port_b.Q_flow = -Q_flow[nSta + 1];
  Q_flow[1] = G[1]*(port_a.T - T[1]);
  Q_flow[nSta + 1] = G[nSta + 1]*(T[nSta] - port_b.T);
  for i in 2:nSta loop
    Q_flow[i] = G[i]*(T[i - 1] - T[i]);
    // Q_flow[i] represents the heat flowing between two nodes
  end for;
  if soiDat.steadyState then
    for i in 2:nSta + 1 loop
      Q_flow[i] = Q_flow[1];
    end for;
  else
    for i in 1:nSta loop
      der(T[i]) = (Q_flow[i] - Q_flow[i + 1])/C[i];
    end for;
  end if;
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                    graphics={
        Rectangle(
          extent={{-94,4},{92,-4}},
          lineColor={0,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{12,8},{14,8},{16,4},{18,-2},{18,-6},{16,-12},{10,-16},{6,-20},
              {-2,-22},{-6,-18},{-12,-12},{-14,-2},{-12,4},{-10,8},{-8,10},{-6,
              12},{-2,14},{2,14},{8,12},{12,8}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-6,-16},{2,-18},{8,-16},{14,-14},{10,-16},{6,-20},{-2,-22},{
              -8,-20},{-12,-12},{-14,-2},{-12,4},{-10,8},{-8,10},{-10,0},{-10,-8},
              {-6,-16}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-110,-74},{-26,-86}},
          lineColor={0,0,255},
          textString="%r_a"),
        Text(
          extent={{-22,-62},{20,-76}},
          lineColor={0,0,255},
          textString="%nSta"),
        Text(
          extent={{16,-76},{102,-88}},
          lineColor={0,0,255},
          textString="%r_b"),
        Polygon(
          points={{-50,60},{-38,34},{-32,0},{-36,-30},{-50,-60},{-62,-60},{-48,
              -30},{-44,0},{-50,34},{-62,60},{-50,60}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillPattern=FillPattern.Backward,
          fillColor={175,175,175}),
        Polygon(
          points={{52,60},{64,34},{70,0},{66,-30},{52,-60},{40,-60},{54,-30},{
              58,0},{52,34},{40,60},{52,60}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillPattern=FillPattern.Backward,
          fillColor={175,175,175}),
        Text(
          extent={{-100,100},{100,60}},
          lineColor={0,0,255},
          textString="%name")}),
    defaultComponentName="lay",
    Documentation(info="<html>
<p>
Model for radial heat transfer in a hollow cylinder.
</p>
<p>
If the heat capacity of the material is non-zero, then this model computes transient heat conduction, i.e., it
computes a numerical approximation to the solution of the heat equation
</p>
<p align=\"center\" style=\"font-style:italic;\">
   &rho; c ( &part; T(r,t) &frasl; &part;t ) = 
    k ( &part;&sup2; T(r,t) &frasl; &part;r&sup2; + 1 &frasl; r &nbsp;  &part; T(r,t) &frasl; &part;r ),
</p>
<p>
where 
<i>&rho;</i>
is the mass density,
<i>c</i>
is the specific heat capacity per unit mass,
<i>T</i>
is the temperature at location <i>r</i> and time <i>t</i> and
<i>k</i> is the heat conductivity. 
At the locations <i>r=r<sub>a</sub></i> and <i>r=r<sub>b</sub></i>, 
the temperature and heat flow rate are equal to the 
temperature and heat flow rate of the heat ports.
</p>
<p>
If the heat capacity of the material is set to zero, then steady-state heat flow is computed using
</p>
<p align=\"center\" style=\"font-style:italic;\">
   Q = 2 &pi; k (T<sub>a</sub>-T<sub>b</sub>)&frasl; ln(r<sub>a</sub> &frasl; r<sub>b</sub>),
</p>
<p>
where
<i>r<sub>a</sub></i> is the internal radius,
<i>r<sub>b</sub></i> is the external radius,
<i>T<sub>a</sub></i> is the temperature at port a and
<i>T<sub>b</sub></i> is the temperature at port b.
</p>
<h4>Implementation</h4>
<p>
To spatially discretize the heat equation, the construction is 
divided into compartments with <code>nSta &ge; 1</code> state variables. 
The state variables are connected to each other through thermal conductors. 
There is also a thermal conductor
between the surfaces and the outermost state variables. Thus, to obtain
the surface temperature, use <code>port_a.T</code> (or <code>port_b.T</code>)
and not the variable <code>T[1]</code>.
</p>
</html>", revisions="<html>
<ul>
<li>
January, 2014, by Damien Picard:<br/>
Modify the discretization of the cilindrical layer so that the first three layers have an equal thickness the following an exponentionally growing thickness.
This follows the guidelines of Eskilson (P. Eskilson. Thermal analysis of heat extraction
boreholes. PhD thesis, Dep. of Mathematical
Physics, University of Lund, Sweden, 1987).
</li>
<li>
March 9, 2012, by Michael Wetter:<br/>
Removed protected variable <code>der_T</code> as it is not required.
</li>
<li>
April 14 2011, by Pierre Vigouroux:<br/>
First implementation.
</li>
</ul>
</html>"));
end EqualCylindrical;
