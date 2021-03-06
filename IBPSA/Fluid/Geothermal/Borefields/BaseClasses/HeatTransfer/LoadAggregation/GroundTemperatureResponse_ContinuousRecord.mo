within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation;
model GroundTemperatureResponse_ContinuousRecord
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

  parameter Modelica.SIunits.Time tStep = 604800
  "Time-step of the long-term predictions";
  parameter Real[16] intervals = {1,2,3,4,8,12,16,20,24,28,32,36,40,44,48,52}
  "Array with the long-term intervals to be evaluated";
  Modelica.SIunits.TemperatureDifference[16-1] delTBor_LT
  "Array with the long-term predictions of deltaT";
  parameter Real[i, 16-1] kappa_LT(each fixed=false)
    "Weight factor for each aggregation cell for long-term predictions";
//   Modelica.SIunits.Time[16] futTime;
  Modelica.Blocks.Interfaces.RealInput curTime(unit="s")
    annotation (Placement(transformation(extent={{-120,-50},{-100,-30}}),
        iconTransformation(extent={{-120,-50},{-100,-30}})));
  Real wT = IBPSA.Utilities.Math.Functions.spliceFunction(1,0, time-(curTime+6*86400), 86400/4)
  "weighting function to take only the last value of optimization";
  Real[16-1] QBor_LT(unit="W")
  "Long-term prediction of the ground loads";

  parameter Data.GFunctions.SquareConfig_9bor_3x3_B6 gFunc
    annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
  parameter Data.GFunctions.SquareConfig_9bor_3x3_B6 gFuncOriginal
    annotation (Placement(transformation(extent={{-8,-80},{12,-60}})));
  Modelica.Blocks.Interfaces.RealOutput delTBorOriginal(unit="K")
    "Temperature difference current borehole wall temperature minus initial borehole wall temperature"
    annotation (Placement(transformation(extent={{100,-66},{126,-40}}),
        iconTransformation(extent={{100,-76},{120,-56}})));


  Modelica.SIunits.HeatFlowRate[15] Qinj
  "Heat flow injected into the field";
  Modelica.SIunits.HeatFlowRate[15] Qgb(min=0)
  "Gas boiler heat flow";
  Modelica.SIunits.HeatFlowRate[15] Qcon(min=0)
  "Condenser heat flow";
  Real[15] COP
  "Heat pump COP";
  Modelica.SIunits.Temperature[15] TevaOutLT
  "Steady state prediction of the outlet evaporator temperature";
  Real[15] costLT
  "long-term cost";
  parameter Real electricityPrice;
  parameter Real gasPrice;


  Modelica.Blocks.Interfaces.RealInput Qext[15]
    annotation (Placement(transformation(extent={{-120,60},{-100,80}})));

 Modelica.Blocks.Interfaces.RealInput[15] Qbuih(unit="W")
    annotation (Placement(transformation(extent={{-120,-78},{-100,-58}}),
        iconTransformation(extent={{-120,-78},{-100,-58}})));
 Modelica.Blocks.Interfaces.RealInput[15] Qbuic(unit="W")
    annotation (Placement(transformation(extent={{-120,-100},{-100,-80}}),
        iconTransformation(extent={{-120,-100},{-100,-80}})));

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
  final parameter Real[76,2] timSer(each fixed=false)
    "g-function input from matrix, with the second column as temperature Tstep";
  final parameter Real[76,2] timSerOriginal(each fixed=false)
    "g-function input from matrix, with the second column as temperature Tstep";
  final parameter Modelica.SIunits.Time[i] nu(each fixed=false)
    "Time vector for load aggregation";
  final parameter Modelica.SIunits.Time[i,16] nu_LT(each fixed=false)
    "Time vector for load aggregation for long-term predictions, i.e. start point is increased";
  final parameter Real[i] kappa(each fixed=false)
    "Weight factor for each aggregation cell";
  final parameter Real[i] kappaOriginal(each fixed=false)
    "Weight factor for each aggregation cell";
  final parameter Real[i] rCel(each fixed=false) "Cell widths";

  parameter Real[16-1,16-1] deltaG(each fixed=false)
  "Evaluation of the g-function at the long-term intervals";

initial equation
  QAgg_flow = zeros(i);
  delTBor = 0;

  (nu,rCel) = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationCellTimes(
    i=i,
    lvlBas=lvlBas,
    nCel=nCel,
    tLoaAgg=tLoaAgg,
    timFin=timFin);

  for j in 1:16 loop
    nu_LT[:,j] = nu[:] + tStep*intervals[j]*ones(i);
  end for;

  kappa = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationWeightingFactors(
    i=i,
    nTimTot=nTimTot,
    TStep=timSer,
    nu=nu);


  kappaOriginal = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationWeightingFactors(
    i=i,
    nTimTot=nTimTot,
    TStep=timSerOriginal,
    nu=nu);

  timSer[:,1] = gFunc.timExp[:];
  timSer[:,2] = gFunc.gFunc[:];

  timSerOriginal[:,1] = gFuncOriginal.timExp[:];
  timSerOriginal[:,2] = gFuncOriginal.gFunc[:];


  // curTime = time;

  for k in 1:16-1 loop
    for j in 1:k loop
      deltaG[j,k] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[k+1]-intervals[j]))
                  - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[k+1]-intervals[j+1]));
    end for;
    for j in k+1:16-1 loop
      deltaG[j,k] = 0;
    end for;
  end for;

 for j in 2:16 loop
       kappa_LT[1,j-1] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[j]-1) + nu[1])
                     - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[j]-1));
      for k in 2:i loop
      kappa_LT[k,j-1] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[j]-1) + nu[k])
                    - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], tStep*(intervals[j]-1) + nu[k-1]);
      end for;
 end for;

//  for j in 1:15 loop
//   Qbuih[j] = 0;
//   Qbuic[j] = 0;
//  end for;

equation
  assert(size(gFunc.timExp,1) == 76, "The size of the time series and the g-function does not match", AssertionLevel.error);
  assert(size(gFunc.timExp,1) == 76, "The size of the time series and the g-function does not match", AssertionLevel.error);
  delTBor = QAgg_flow[:]*kappa[:];
  delTBorOriginal = QAgg_flow[:]*kappaOriginal[:];

  QBor_LT = Qinj + Qext;

  Qinj = -Qbuic;

  Qbuih = Qgb + Qcon;

  Qext = (-Qcon).*(COP-ones(15))./COP;

  COP = 5.15*ones(15) + 0.1*delTBor_LT;

  TevaOutLT = 4.46*ones(15) + (6/7)*delTBor_LT;

 //  curTime = time;
  //  der(curTime) = 0;
//   for j in 1:16 loop
//      futTime[j] = der(QAgg_flow[j]curTime + tStep*intervals[j];
//   end for;

   for j in 2:16 loop
      //Qbuih[j-1] = sum(IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolateVector(Qbui[:,1], Qbui[:,2], mod((time + (intervals[j-1]-1)*tStep)/86400 + k,365)) for k in 1:(intervals[j]-intervals[j-1])*tStep/86400)/((intervals[j] - intervals[j - 1])*tStep/86400);
      //Qbuic[j-1] = sum(IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolateVector(Qbui[:,1], Qbui[:,3], mod((time + (intervals[j-1]-1)*tStep)/86400 + k,365)) for k in 1:(intervals[j]-intervals[j-1])*tStep/86400)/((intervals[j] - intervals[j - 1])*tStep/86400);
      costLT[j-1] = (gasPrice/1000*Qgb[j-1] + electricityPrice/1000*((Qcon[j-1]/COP[j-1])))*((intervals[j] - intervals[j - 1])*tStep/3600);
   end for;

//         Qbuih[j-1] =sum(buiNeeds.Qbuih[integer(mod(time + intervals[j - 1]*
//       tStep + 3600*k, 31536000)/3600 + 1)] for k in 1:(intervals[j] - intervals[
//       j - 1])*tStep/3600)/((intervals[j] - intervals[j - 1])*tStep/3600);
//         Qbuic[j-1] =sum(buiNeeds.Qbuih[integer(mod(time + intervals[j - 1]*
//       tStep + 3600*k, 31536000)/3600 + 1)] for k in 1:(intervals[j] - intervals[
//       j - 1])*tStep/3600)/((intervals[j] - intervals[j - 1])*tStep/3600);

// for j in 1:16-1 loop
//       kappa_LT[1,j] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j+1] - time + nu[1])
//                     - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j+1] - time);
//       for k in 2:i loop
//       kappa_LT[k,j] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j+1] - time + nu[k])
//                     - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j+1] - time + nu[k-1]);
//     end for;
// end for;


 // for j in 1:16-1 loop
 //      kappa_LT[1,j] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j] - time + nu[1])
 //                    - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j] - time);
 //      for k in 2:i loop
 //      kappa_LT[k,j] = IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j] - time + nu[k])
 //                    - IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.interpolate(timSer[:,1], timSer[:,2], futTime[j] - time + nu[k-1]);
 //    end for;
 // end for;

  for i in 1:16-1 loop
    delTBor_LT[i] = QAgg_flow[:]*kappa_LT[:,i] + QBor_LT[:]*deltaG[:,i];
  end for;



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
</html>"),
    experiment(
      StopTime=315360000,
      Interval=300,
      Tolerance=1e-06,
      __Dymola_Algorithm="Euler"));
end GroundTemperatureResponse_ContinuousRecord;
