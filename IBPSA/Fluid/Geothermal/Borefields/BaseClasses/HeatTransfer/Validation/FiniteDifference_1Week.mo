within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Validation;
model FiniteDifference_1Week
  "Short term validation of ground temperature response model"
  extends Modelica.Icons.Example;

  parameter Modelica.SIunits.Temperature T_start = 283.15
    "Initial soil temperature";
  IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Cylindrical soi(
    final steadyStateInitial=false,
    final soiDat=borFieDat.soiDat,
    final r_a=borFieDat.conDat.rBor,
    final TInt_start=T_start,
    final TExt_start=T_start,
    gridFac=1.2,
    final nSta=50,
    final r_b=15,
    final h=borFieDat.conDat.hBor)
                   "Heat conduction in the soil"
    annotation (Placement(transformation(extent={{-12,50},{8,70}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo
    "Prescribed heat flow to soil for the discretized model"
    annotation (Placement(transformation(extent={{-50,50},{-30,70}})));

  parameter IBPSA.Fluid.Geothermal.Borefields.Data.Borefield.Template borFieDat(
    soiDat(
      kSoi=1,
      cSoi=1,
      dSoi=1e6),
    filDat(
      kFil=0,
      cFil=Modelica.Constants.small,
      dFil=Modelica.Constants.small,
      steadyState=true),
    conDat(
      borCon=Types.BoreholeConfiguration.SingleUTube,
      mBor_flow_nominal=0.3,
      dp_nominal=5e4,
      rBor=0.05,
      dBor=4,
      rTub=0.02,
      kTub=0.5,
      eTub=0.002,
      xC=0.05,
      nBor=4,
      cooBor=[0,0; 0,6; 6,0; 6,6],
      hBor=100))
                "Borefield parameters"
    annotation (Placement(transformation(extent={{-90,-88},{-70,-68}})));

  Modelica.Blocks.Sources.Constant
                               const(k=2*borFieDat.soiDat.kSoi*Modelica.Constants.pi
        *borFieDat.conDat.hBor*4)
                   "Heat flow signal"
    annotation (Placement(transformation(extent={{-92,-10},{-72,10}})));
  IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse
    groTemRes(
    borFieDat=borFieDat,
    nCel=5,
    r=soi.rC,
    tLoaAgg=300,
    forceGFunCalc=true)
              "Heat conduction in the soil"
    annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
  Modelica.Blocks.Sources.Constant groTem(k=T_start)
    "Ground temperature signal"
    annotation (Placement(transformation(extent={{-94,-46},{-74,-26}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature preTem
    "Ground temperature for discretized model"
    annotation (Placement(transformation(extent={{48,50},{28,70}})));
  Modelica.Blocks.Math.Add deltaT(
    k2=-1,
    y(unit="K"))
    "Temperature difference between borehole with discrete ground and borehole with analytical ground"
    annotation (Placement(transformation(extent={{30,-10},{50,10}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temSenDis
    "Borehole wall temperature sensor for the discretized model" annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={0,20})));
  Modelica.Blocks.Math.Add[51] TSoiAna(y(unit="K"))
    "Ground temperature from ground response model"
    annotation (Placement(transformation(extent={{-12,-40},{8,-20}})));
  Modelica.Blocks.Routing.Replicator                reaRep(nout=51)
    annotation (Placement(transformation(extent={{-48,-46},{-28,-26}})));
equation
  connect(preHeaFlo.port, soi.port_a) annotation (Line(
      points={{-30,60},{-12,60}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(preTem.port, soi.port_b)
    annotation (Line(points={{28,60},{18,60},{8,60}}, color={191,0,0}));
  connect(preTem.T, groTem.y)
    annotation (Line(points={{50,60},{60,60},{60,88},{-66,88},{-66,-36},{-73,
          -36}},                                             color={0,0,127}));
  connect(temSenDis.T, deltaT.u1)
    annotation (Line(points={{10,20},{20,20},{20,6},{28,6}}, color={0,0,127}));
  connect(soi.port_a, temSenDis.port) annotation (Line(points={{-12,60},{-20,60},
          {-20,20},{-10,20}}, color={191,0,0}));

  connect(groTemRes.QBor_flow, const.y)
    annotation (Line(points={{-51,0},{-71,0}}, color={0,0,127}));
  connect(groTemRes.delTBor, TSoiAna.u1) annotation (Line(points={{-29,0},{-20,
          0},{-20,-24},{-14,-24}}, color={0,0,127}));
  connect(groTem.y, reaRep.u)
    annotation (Line(points={{-73,-36},{-50,-36}}, color={0,0,127}));
  connect(reaRep.y, TSoiAna.u2)
    annotation (Line(points={{-27,-36},{-14,-36}}, color={0,0,127}));
  connect(TSoiAna[1].y, deltaT.u2) annotation (Line(points={{9,-30},{20,-30},{
          20,-6},{28,-6}}, color={0,0,127}));
  connect(const.y, preHeaFlo.Q_flow) annotation (Line(points={{-71,0},{-60,0},{
          -60,60},{-50,60}}, color={0,0,127}));
  annotation (
    __Dymola_Commands(file=
          "modelica://IBPSA/Resources/Scripts/Dymola/Fluid/Geothermal/Borefields/BaseClasses/HeatTransfer/Validation/FiniteDifference_1Week.mos"
        "Simulate and plot", file="Resources/Scripts/Dymola/Fluid/Geothermal/Borefields/BaseClasses/HeatTransfer/LoadAggregation/Validation/FiniteDifference_1Week.mos"),
    experiment(Tolerance=1e-6, StopTime=604800.0),
    Documentation(info="<html>
<p>
This example validates the implementation of
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse</a>
for the evaluation of the borehole wall temperature at a short time scale.
</p>
<p>
After a short delay, a sinusoidal heat flow rate is applied to borehole heat
exchanger. The temperature at the borehole wall evaluated with
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse</a>
is compared to the temperature obtained with
<a href=\"modelica://IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Cylindrical\">
IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.Cylindrical</a>
</p>
</html>", revisions="<html>
<ul>
<li>
June 13, 2018, by Massimo Cimmino:<br/>
First implementation.
</li>
</ul>
</html>"));
end FiniteDifference_1Week;
