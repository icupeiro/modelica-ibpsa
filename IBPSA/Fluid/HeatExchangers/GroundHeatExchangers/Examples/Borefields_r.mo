within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Examples;
model Borefields_r
  "Borefield with a double U-Tube configuration in parallel"
extends Modelica.Icons.Example;
  package Medium = IBPSA.Media.Water;

  parameter Modelica.SIunits.Time tLoaAgg=60
    "Time resolution of load aggregation";

  Data.BorefieldData.ExampleBorefieldData borFieUTubDat(conDat=
        IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Data.ConfigurationData.ExampleConfigurationData(
        borHolCon=IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.Types.BoreHoleConfiguration.SingleUTube))
    annotation (Placement(transformation(extent={{80,-100},{100,-80}})));
  Modelica.Blocks.Sources.Constant TGro(k=283.15) "Ground temperature"
    annotation (Placement(transformation(extent={{-80,30},{-60,50}})));

  BorefieldOneUTube_r                                                borFieUTub(
    redeclare package Medium = Medium,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    dynFil=false,
    r={borFieUTubDat.conDat.rBor*2,borFieUTubDat.conDat.rBor*4,borFieUTubDat.conDat.rBor
        *8,borFieUTubDat.conDat.rBor*16})
                "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-22,-78},{20,-42}})));
  Sources.MassFlowSource_T sou(
    redeclare package Medium = Medium,
    nPorts=1,
    use_T_in=false,
    m_flow=borFieUTubDat.conDat.m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-100,-70},{
            -80,-50}}, rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn(redeclare package Medium = Medium,
      m_flow_nominal=borFieUTubDat.conDat.m_flow_nominal)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,-70},{-40,-50}})));
  Sources.Boundary_pT sin(
    redeclare package Medium = Medium,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{100,-70},{80,
            -50}}, rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut(redeclare package Medium = Medium,
      m_flow_nominal=borFieUTubDat.conDat.m_flow_nominal)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,-70},{60,-50}})));

equation
  connect(sou.ports[1], TUTubIn.port_a)
    annotation (Line(points={{-80,-60},{-60,-60}}, color={0,127,255}));
  connect(TUTubIn.port_b, borFieUTub.port_a)
    annotation (Line(points={{-40,-60},{-22,-60}}, color={0,127,255}));
  connect(borFieUTub.port_b, TUTubOut.port_a)
    annotation (Line(points={{20,-60},{30,-60},{40,-60}}, color={0,127,255}));
  connect(TUTubOut.port_b, sin.ports[1])
    annotation (Line(points={{60,-60},{70,-60},{80,-60}}, color={0,127,255}));
  connect(TGro.y, borFieUTub.TGro) annotation (Line(points={{-59,40},{-59,40},{-36,
          40},{-36,-49.2},{-26.2,-49.2}}, color={0,0,127}));
  annotation (__Dymola_Commands(file="Resources/Scripts/Dymola/Fluid/HeatExchangers/GroundHeatExchangers/Examples/Borefields_r.mos"
        "Simulate and Plot"));
end Borefields_r;
