within IBPSA.Fluid.Geothermal.Borefields;
model OneUTubeContinuous
  extends OneUTube(redeclare
      BaseClasses.HeatTransfer.GroundTemperatureResponse_ContinuousRecord
      groTemRes);
  Modelica.Blocks.Interfaces.RealOutput QBor_flow(
    final unit="W",
    displayUnit="W",
    start=0) "Heat flow of ground"
    annotation (Placement(transformation(extent={{100,80},{120,100}})));
  Modelica.Blocks.Interfaces.RealOutput deltaT(
    final quantity="ThermodynamicTemperature",
    final unit="K",
    displayUnit="degC",
    start=0)
    "temperature difference between the original g-function and the cyclic one"
    annotation (Placement(transformation(extent={{100,58},{120,78}})));
protected
  Modelica.Blocks.Math.Add deltaTgFunc(each final k1=1, each final k2=-1)
    annotation (Placement(transformation(extent={{80,52},{94,66}})));
equation
  connect(gaiQ_flow.y, QBor_flow) annotation (Line(points={{1,80},{6,80},{6,96},
          {92,96},{92,90},{110,90}}, color={0,0,127}));
  connect(deltaTgFunc.y, deltaT) annotation (Line(points={{94.7,59},{100.35,59},
          {100.35,68},{110,68}}, color={0,0,127}));
  connect(deltaTgFunc.u1, groTemRes.delTBor) annotation (Line(points={{78.6,
          63.2},{52,63.2},{52,80},{41,80}}, color={0,0,127}));
  connect(deltaTgFunc.u2, groTemRes.delTBorOriginal) annotation (Line(points={{
          78.6,54.8},{46,54.8},{46,73.4},{41,73.4}}, color={0,0,127}));
end OneUTubeContinuous;
