within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer;
model GroundTemperatureResponse_ContinuousCSV
  "Model calculating discrete load aggregation"
  parameter Modelica.SIunits.Time tLoaAgg(final min = Modelica.Constants.eps)=3600
    "Time resolution of load aggregation";
  parameter Integer nCel(min=1)=5 "Number of cells per aggregation level";
  parameter Boolean forceGFunCalc = false
    "Set to true to force the thermal response to be calculated at the start instead of checking whether it has been pre-computed";
  parameter IBPSA.Fluid.Geothermal.Borefields.Data.Borefield.Template borFieDat
    "Record containing all the parameters of the borefield model" annotation (
     choicesAllMatching=true, Placement(transformation(extent={{-80,-80},{-60,-60}})));

  Modelica.Blocks.Interfaces.RealOutput delTBor(unit="K")
    "Temperature difference current borehole wall temperature minus initial borehole wall temperature"
    annotation (Placement(transformation(extent={{100,-14},{126,12}}),
        iconTransformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput QBor_flow(unit="W")
    "Heat flow from all boreholes combined (positive if heat from fluid into soil)"
    annotation (Placement(transformation(extent={{-120,-10},{-100,10}}),
        iconTransformation(extent={{-120,-10},{-100,10}})));

  Modelica.SIunits.HeatFlowRate[i] QAgg_flow
    "Vector of aggregated loads";
   Modelica.SIunits.HeatFlowRate[i] QFace
     "Vector of cell face values of aggregated loads";

protected
  constant Integer nSegMax = 1500 "Max total number of segments in g-function calculation";
  final parameter Integer nSeg = integer(if 12*borFieDat.conDat.nBor<nSegMax then 12 else floor(nSegMax/borFieDat.conDat.nBor))
    "Number of segments per borehole for g-function calculation";
  constant Integer nTimSho = 26 "Number of time steps in short time region";
  constant Integer nTimLon = 50 "Number of time steps in long time region";
  constant Real ttsMax = exp(5) "Maximum non-dimensional time for g-function calculation";
  constant Integer nTimTot = nTimSho+nTimLon
    "Total length of g-function vector";
  constant Real lvlBas = 2 "Base for exponential cell growth between levels";

  parameter String SHAgfun=
    IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.shaGFunction(
      nBor=borFieDat.conDat.nBor,
      cooBor=borFieDat.conDat.cooBor,
      hBor=borFieDat.conDat.hBor,
      dBor=borFieDat.conDat.dBor,
      rBor=borFieDat.conDat.rBor,
      aSoi=borFieDat.soiDat.aSoi,
      nSeg=nSeg,
      nTimSho=nTimSho,
      nTimLon=nTimLon,
      ttsMax=ttsMax) "String with encrypted g-function arguments";
  parameter Modelica.SIunits.Time timFin=
    (borFieDat.conDat.hBor^2/(9*borFieDat.soiDat.aSoi))*ttsMax
    "Final time for g-function calculation";
  parameter Integer i(min=1)=
    IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.countAggregationCells(
      lvlBas=lvlBas,
      nCel=nCel,
      timFin=timFin,
      tLoaAgg=tLoaAgg)
      "Number of aggregation cells";
  final parameter Real[nTimTot,2] timSer(each fixed=false)
    "g-function input from matrix, with the second column as temperature Tstep";
  final parameter Modelica.SIunits.Time[i] nu(each fixed=false)
    "Time vector for load aggregation";
  final parameter Real[i] kappa(each fixed=false)
    "Weight factor for each aggregation cell";
  final parameter Real[i] rCel(each fixed=false) "Cell widths";

initial equation
  QAgg_flow = zeros(i);
  delTBor = 0;

  (nu,rCel) = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationCellTimes(
    i=i,
    lvlBas=lvlBas,
    nCel=nCel,
    tLoaAgg=tLoaAgg,
    timFin=timFin);

  kappa = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationWeightingFactors(
    i=i,
    nTimTot=nTimTot,
    TStep=timSer,
    nu=nu);

  timSer = transpose({{
0.00000000000000E+00,
1.00000000000000E+00,
2.30986142530842E+00,
4.02559837881945E+00,
6.27297603019975E+00,
9.21672932384306E+00,
1.30726382088110E+01,
1.81233445167347E+01,
2.47390698800458E+01,
3.34047533338816E+01,
4.47555978139947E+01,
5.96236311431699E+01,
7.90986944712566E+01,
1.04608328680156E+02,
1.38022414504122E+02,
1.81790236586880E+02,
2.39120018402847E+02,
3.14214088124931E+02,
4.12576913323311E+02,
5.41418583735025E+02,
7.10183317779631E+02,
9.31241732857096E+02,
1.22079762350688E+03,
1.60007571513987E+03,
2.09687745683452E+03,
2.74761889430638E+03,
3.60000000000000E+03,
8.70005510809538E+03,
1.59252112485425E+04,
2.61609597154520E+04,
4.06617600637413E+04,
6.12047825325647E+04,
9.03077121652983E+04,
1.31537307980288E+05,
1.89946533186106E+05,
2.72693829678003E+05,
3.89920433054276E+05,
5.55993248985895E+05,
7.91265613765851E+05,
1.12457173204076E+06,
1.59676050178171E+06,
2.26570182042133E+06,
3.21337892851689E+06,
4.55593600517463E+06,
6.45791241537404E+06,
9.15240811157249E+06,
1.29696515946663E+07,
1.83774716293321E+07,
2.60386327936259E+07,
3.68920617185977E+07,
5.22679188378067E+07,
7.40506184606309E+07,
1.04909776371600E+08,
1.48627389131424E+08,
2.10561343095166E+08,
2.98302059283313E+08,
4.22602750334232E+08,
5.98697298760349E+08,
8.48167271318230E+08,
1.20158688461274E+09,
1.70227008019156E+09,
2.41157893828858E+09,
3.41644401188720E+09,
4.84001824842979E+09,
6.85677020866089E+09,
9.71386635769923E+09,
1.37614629713721E+10,
1.94956201338922E+10,
2.76190972248179E+10,
3.91274807895598E+10,
5.54312003410274E+10,
7.85283859464210E+10,
1.11249752454591E+11,
1.57605522566774E+11,
2.23276906497422E+11,
3.16312372342685E+11},
{0.00000000000000E+00,
4.67327338532371E-07,
8.54768300010403E-07,
1.22029565548485E-06,
1.59251031763735E-06,
1.98697309251237E-06,
2.41445796668285E-06,
2.88406029167041E-06,
3.40440409155404E-06,
3.98419491726270E-06,
4.63251421758110E-06,
5.35898913933792E-06,
6.17389484110219E-06,
7.08819158607969E-06,
8.11321867945619E-06,
9.25823070408750E-06,
1.05400294443549E-05,
1.19730758134009E-05,
1.35693615812436E-05,
1.53427286990428E-05,
1.73068512171728E-05,
1.94749222298288E-05,
2.18592994659176E-05,
2.44711118336247E-05,
2.73198635409048E-05,
3.04130000718121E-05,
3.37555601770054E-05,
4.64223376139180E-05,
5.65805631056040E-05,
6.57201179283028E-05,
7.43613769559089E-05,
8.27428198248109E-05,
9.09876989599182E-05,
9.91630816597448E-05,
1.07305720197607E-04,
1.15435054729307E-04,
1.23560361178599E-04,
1.31684898598120E-04,
1.39808379624034E-04,
1.47928448499937E-04,
1.56103303544163E-04,
1.64584136879976E-04,
1.73944658472606E-04,
1.85089273402192E-04,
1.99040192856377E-04,
2.16775190419390E-04,
2.39078926162446E-04,
2.66550182055831E-04,
2.99477731521920E-04,
3.37701737317360E-04,
3.80555853455346E-04,
4.26964673489830E-04,
4.75624618772444E-04,
5.25195558036758E-04,
5.74432794164806E-04,
6.22255973142679E-04,
6.67767032100423E-04,
7.10238796173440E-04,
7.49089519389505E-04,
7.83844644834656E-04,
8.14106658191637E-04,
8.39578588492022E-04,
8.60163998578101E-04,
8.76066786803901E-04,
8.87797652828870E-04,
8.96080076601125E-04,
9.01704948487089E-04,
9.05402615658912E-04,
9.07770729031246E-04,
9.09256929057555E-04,
9.10175443894487E-04,
9.10736668309760E-04,
9.11076720260699E-04,
9.11281508981634E-04,
9.11404296990298E-04,
9.11477686416417E-04}});

equation
  delTBor = QAgg_flow[:]*kappa[:];

//    // "Upwind" scheme
//    der(QAgg_flow[1]) = -1/(rCel[1]*tLoaAgg)*(QAgg_flow[1] - QBor_flow);
//    for j in 2:i-1 loop
//      der(QAgg_flow[j]) = -1/(rCel[j]*tLoaAgg)*(QAgg_flow[j] - QAgg_flow[j-1]);
//    end for;
//    der(QAgg_flow[i]) = 1/(rCel[i]*tLoaAgg)*(QAgg_flow[i-1]);

   // "QUICK" scheme
   QFace[1] = QBor_flow;
   QFace[2] = 0.5*(QAgg_flow[1] + QAgg_flow[2]) - 0.125*(0.5*rCel[1] + 0.5*rCel[2])^2/rCel[1]*((QAgg_flow[2] - QAgg_flow[1])/(0.5*rCel[1] + 0.5*rCel[2]) - (QAgg_flow[1] - QBor_flow)/(0.5*rCel[1]));
   for j in 3:i loop
     QFace[j] = 0.5*(QAgg_flow[j-1] + QAgg_flow[j]) - 0.125*(0.5*rCel[j-1] + 0.5*rCel[j])^2/rCel[j-1]*((QAgg_flow[j] - QAgg_flow[j-1])/(0.5*rCel[j-1] + 0.5*rCel[j]) - (QAgg_flow[j-1] - QAgg_flow[j-2])/(0.5*rCel[j-1] + 0.5*rCel[j-2]));
   end for;
   for j in 1:i-1 loop
     der(QAgg_flow[j])*rCel[j]*tLoaAgg = QFace[j] - QFace[j+1];
   end for;
   der(QAgg_flow[i]) = 1/(rCel[i]*tLoaAgg)*(QAgg_flow[i-1]);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-100,30},{100,-100}},
          lineColor={0,0,0},
          fillColor={127,127,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-52,30},{-94,-100}},
          lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-66,-4},{72,-4}},
          color={255,0,0},
          arrow={Arrow.None,Arrow.Filled}),
        Rectangle(
          extent={{-100,30},{-94,-100}},
          lineColor={0,0,0},
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-151,147},{149,107}},
          lineColor={0,0,255},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255},
            textString="%name")}),
  Diagram(
        coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>
This model calculates the ground temperature response to obtain the temperature
at the borehole wall in a geothermal system where heat is being injected into or
extracted from the ground.
</p>
<p>
A load-aggregation scheme based on that developed by Claesson and Javed (2012) is
used to calculate the borehole wall temperature response with the temporal superposition
of ground thermal loads. In its base form, the
load-aggregation scheme uses fixed-length aggregation cells to agglomerate
thermal load history together, with more distant cells (denoted with a higher cell and vector index)
representing more distant thermal history. The more distant the thermal load, the
less impactful it is on the borehole wall temperature change at the current time step.
Each cell has an <em>aggregation time</em> associated to it denoted by <code>nu</code>,
which corresponds to the simulation time (since the beginning of heat injection or
extraction) at which the cell will begin shifting its thermal load to more distant
cells. To determine <code>nu</code>, cells have a temporal size <i>r<sub>cel</sub></i>
(<code>rcel</code> in this model)
which follows the exponential growth
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_02.png\" />
</p>
<p>
where <i>n<sub>Cel</sub></i> is the number of consecutive cells which can have the same size.
Decreasing <i>r<sub>cel</sub></i> will generally decrease calculation times, at the cost of
precision in the temporal superposition. <code>rcel</code> is expressed in multiples
of the aggregation time resolution (via the parameter <code>tLoaAgg</code>).
Then, <code>nu</code> may be expressed as the sum of all <code>rcel</code> values
(multiplied by the aggregation time resolution) up to and including that cell in question.
</p>
<p>
To determine the weighting factors, the borefield's temperature
step response at the borefield wall is determined as
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_03.png\" />
</p>
<p>
where <i>g(&middot;)</i> is the borefield's thermal response factor known as the <em>g-function</em>,
<i>H</i> is the total length of all boreholes and <i>k<sub>s</sub></i> is the thermal
conductivity of the soil. The weighting factors <code>kappa</code> (<i>&kappa;</i> in the equation below)
for a given cell <i>i</i> are then expressed as follows.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_04.png\" />
</p>
<p>
where <i>&nu;</i> refers to the vector <code>nu</code> in this model and
<i>T<sub>step</sub>(&nu;<sub>0</sub>)</i>=0.
</p>
<p>
At every aggregation time step, a time event is generated to perform the load aggregation steps.
First, the thermal load is shifted. When shifting between cells of different size, total
energy is conserved. This operation is illustred in the figure below by Cimmino (2014).
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_01.png\" />
</p>
<p>
After the cell-shifting operation is performed, the first aggregation cell has its
value set to the average thermal load since the last aggregation step.
Temporal superposition is then applied by means
of a scalar product between the aggregated thermal loads <code>QAgg_flow</code> and the
weighting factors <i>&kappa;</i>.
</p>
<p>
Due to Modelica's variable time steps, the load aggregation scheme is modified by separating
the thermal response between the current aggregation time step and everything preceding it.
This is done according to
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_05.png\" />
<br/>
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_06.png\" />
</p>
<p>
where <i>T<sub>b</sub></i> is the borehole wall temperature,
<i>T<sub>g</sub></i>
is the undisturbed ground temperature,
<i>Q</i> is the ground thermal load per borehole length and <i>h = g/(2 &pi; k<sub>s</sub>)</i>
is a temperature response factor based on the g-function. <i>t<sub>k</sub></i>
is the last discrete aggregation time step, meaning that the current time <i>t</i>
satisfies <i>t<sub>k</sub>&le;t&le;t<sub>k+1</sub></i>.
<i>&Delta;t<sub>agg</sub>(=t<sub>k+1</sub>-t<sub>k</sub>)</i> is the
parameter <code>tLoaAgg</code> in the present model.
</p>
<p>
Thus,
<i>&Delta;T<sub>b</sub>*(t)</i>
is the borehole wall temperature change due to the thermal history prior to the current
aggregation step. At every aggregation time step, load aggregation and temporal superposition
are used to calculate its discrete value. Assuming no heat injection or extraction until
<i>t<sub>k+1</sub></i>, this term is assumed to have a linear
time derivative, which is given by the difference between <i>&Delta;T<sub>b</sub>*(t<sub>k+1</sub>)</i>
(the temperature change from load history at the next discrete aggregation time step, which
is constant over the duration of the ongoing aggregation time step) and the total
temperature change at the last aggregation time step, <i>&Delta;T<sub>b</sub>(t)</i>.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_09.png\" />
</p>
<p>
The second term <i>&Delta;T<sub>b,q</sub>(t)</i> concerns the ongoing aggregation time step.
To obtain the time derivative of this term, the thermal response factor <i>h</i> is assumed
to vary linearly over the course of an aggregation time step. Therefore, because
the ongoing aggregation time step always concerns the first aggregation cell, its derivative (denoted
by the parameter <code>dTStepdt</code> in this model) can be calculated as
<code>kappa[1]</code>, the first value in the <code>kappa</code> vector,
divided by the aggregation time step <i>&Delta;t</i>.
The derivative of the temperature change at the borehole wall is then expressed
as the multiplication of <code>dTStepdt</code> (which only needs to be
calculated once at the start of the simulation) and the heat flow <i>Q</i> at
the borehole wall.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_10.png\" />
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_11.png\" />
</p>
<p>
With the two terms in the expression of <i>&Delta;T<sub>b</sub>(t)</i> expressed
as time derivatives, <i>&Delta;T<sub>b</sub>(t)</i> can itself also be
expressed as its time derivative and implemented as such directly in the Modelica
equations block with the <code>der()</code> operator.
</p>
<p align=\"center\">
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_07.png\" />
<br/>
<img alt=\"image\" src=\"modelica://IBPSA/Resources/Images/Fluid/Geothermal/Borefields/LoadAggregation_08.png\" />
</p>
<p>
This load aggregation scheme is validated in
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Validation.Analytic_20Years\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Validation.Analytic_20Years</a>.
</p>
<h4>References</h4>
<p>
Cimmino, M. 2014. <i>D&eacute;veloppement et validation exp&eacute;rimentale de facteurs de r&eacute;ponse
thermique pour champs de puits g&eacute;othermiques</i>,
Ph.D. Thesis, &Eacute;cole Polytechnique de Montr&eacute;al.
</p>
<p>
Claesson, J. and Javed, S. 2012. <i>A load-aggregation method to calculate extraction temperatures of borehole heat exchangers</i>. ASHRAE Transactions 118(1): 530-539.
</p>
</html>", revisions="<html>
<ul>
<li>
August 30, 2018, by Michael Wetter:<br/>
Refactored model to compute the temperature difference relative to the initial temperature,
because the model is independent of the initial temperature.
</li>
<li>
April 5, 2018, by Alex Laferri&egrave;re:<br/>
First implementation.
</li>
</ul>
</html>"));
end GroundTemperatureResponse_ContinuousCSV;
