within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers;
model BorefieldTwoUTube_r
  extends BorefieldTwoUTube(redeclare
      GroundHeatTransfer.GroundTemperatureResponse_r groTemRes(r=r));
  Modelica.Blocks.Interfaces.RealOutput TSoi[size(groTemRes.TSoi, 1)](each unit="K",each displayUnit="degC")
    annotation (Placement(transformation(extent={{100,44},{120,64}})));
  parameter Modelica.SIunits.Distance r[:]
    "Radial distance from borehole wall at which the soil temperature is evaluated";
equation
  connect(groTemRes.TSoi, TSoi)
    annotation (Line(points={{-59,54},{110,54}}, color={0,0,127}));
end BorefieldTwoUTube_r;
