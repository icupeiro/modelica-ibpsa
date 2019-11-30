within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.Boreholes;
model TwoUTubeLin "Double U-tube borehole heat exchanger"
  extends
    IBPSA.Fluid.Geothermal.Borefields.BaseClasses.Boreholes.BaseClasses.PartialBoreholeLin;

  IBPSA.Fluid.Geothermal.Borefields.BaseClasses.Boreholes.BaseClasses.InternalHEXTwoUTube
    intHex(
    redeclare each final package Medium = Medium,
    final borFieDat=borFieDat,
    final hSeg=borFieDat.conDat.hBor/nSeg,
    final dp1_nominal=dp_nominal,
    final dp3_nominal=dp_nominal,
    final dp2_nominal=dp_nominal,
    final dp4_nominal=dp_nominal,
    final show_T=show_T,
    final energyDynamics=energyDynamics,
    final m1_flow_nominal=m_flow_nominal/2,
    final m2_flow_nominal=m_flow_nominal/2,
    final m3_flow_nominal=m_flow_nominal/2,
    final m4_flow_nominal=m_flow_nominal/2,
    final m1_flow_small= borFieDat.conDat.mBor_flow_small/2,
    final m2_flow_small= borFieDat.conDat.mBor_flow_small/2,
    final m3_flow_small= borFieDat.conDat.mBor_flow_small/2,
    final m4_flow_small= borFieDat.conDat.mBor_flow_small/2,
    final dynFil=dynFil,
    final mSenFac=mSenFac,
    final allowFlowReversal1=allowFlowReversal,
    final allowFlowReversal2=allowFlowReversal,
    final allowFlowReversal3=allowFlowReversal,
    final allowFlowReversal4=allowFlowReversal,
    final from_dp1=from_dp,
    final linearizeFlowResistance1=linearizeFlowResistance,
    final deltaM1=deltaM,
    final from_dp2=from_dp,
    final linearizeFlowResistance2=linearizeFlowResistance,
    final deltaM2=deltaM,
    final from_dp3=from_dp,
    final linearizeFlowResistance3=linearizeFlowResistance,
    final deltaM3=deltaM,
    final from_dp4=from_dp,
    final linearizeFlowResistance4=linearizeFlowResistance,
    final deltaM4=deltaM,
    final p1_start=p_start,
    final p2_start=p_start,
    final p3_start=p_start,
    final p4_start=p_start,
    final TFlu_start=TFlu_start,
    final TGro_start=TGro_start) "Discretized borehole segments"
    annotation (Placement(transformation(extent={{-10,-30},{10,10}})));

equation
  // Couple borehole port_a and port_b to first borehole segment.
  connect(port_a, intHex.port_a1) annotation (Line(
      points={{-100,5.55112e-016},{-52,5.55112e-016},{-52,6},{-10,6}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(port_b, intHex.port_b4) annotation (Line(
      points={{100,5.55112e-016},{28,5.55112e-016},{28,-40},{-32,-40},{-32,-27},
          {-10,-27}},
      color={0,127,255},
      smooth=Smooth.None));
  if borFieDat.conDat.borCon == Types.BoreholeConfiguration.DoubleUTubeParallel then
    // 2U-tube in parallel: couple both U-tube to each other.
    connect(port_a, intHex.port_a3) annotation (Line(
        points={{-100,5.55112e-016},{-52,5.55112e-016},{-52,-16.4},{-10,-16.4}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(port_b, intHex.port_b2) annotation (Line(
        points={{100,5.55112e-016},{28,5.55112e-016},{28,-40},{-32,-40},{-32,-4},
            {-10,-4}},
        color={0,127,255},
        smooth=Smooth.None));
  elseif borFieDat.conDat.borCon == Types.BoreholeConfiguration.DoubleUTubeSeries then
    // 2U-tube in serie: couple both U-tube to each other.
    connect(intHex.port_b2, intHex.port_a3) annotation (Line(
        points={{-10,-4},{-24,-4},{-24,-16},{-18,-16},{-18,-16.4},{-10,-16.4}},
        color={0,127,255},
        smooth=Smooth.None));
  end if;

  // Couple each layer to the next one

    connect(intHex.port_b1, intHex.port_a1) annotation (Line(
        points={{10,6},{10,10},{-10,10},{-10,6}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(intHex.port_a2, intHex.port_b2) annotation (Line(
        points={{10,-4},{10,0},{-10,0},{-10,-4}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(intHex.port_b3, intHex.port_a3) annotation (Line(
        points={{10,-16.2},{10,-12},{-10,-12},{-10,-16.4}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(intHex.port_a4, intHex.port_b4) annotation (Line(
        points={{10,-26},{10,-22},{-10,-22},{-10,-27}},
        color={0,127,255},
        smooth=Smooth.None));


  // Close U-tube at bottom layer
  connect(intHex.port_b1, intHex.port_a2)
    annotation (Line(
      points={{10,6},{16,6},{16,-4},{10,-4}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(intHex.port_b3, intHex.port_a4)
    annotation (Line(
      points={{10,-16.2},{14,-16.2},{14,-16},{18,-16},{18,-26},{10,-26}},
      color={0,127,255},
      smooth=Smooth.None));

  connect(intHex.port_wall, port_wall)
    annotation (Line(points={{0,10},{0,10},{0,100}}, color={191,0,0}));
  annotation (
    defaultComponentName="borHol",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2},
        initialScale=0.5), graphics={
        Rectangle(
          extent={{58,88},{50,-92}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-52,-92},{-44,88}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-50,-84},{56,-92}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{22,88},{14,-92}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-10,88},{-18,-92}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2},
        initialScale=0.5), graphics={Text(
          extent={{60,72},{84,58}},
          lineColor={0,0,255},
          textString=""), Text(
          extent={{50,-32},{90,-38}},
          lineColor={0,0,255},
          textString="")}),
    Documentation(info="<html>
<p>
Model of a double U-tube borehole heat exchanger. 
The borehole heat exchanger is vertically discretized into
<i>n<sub>seg</sub></i> elements of height
<i>h=h<sub>Bor</sub>&frasl;n<sub>seg</sub></i>.
Each segment contains a model for the heat transfer in the borehole, 
with a uniform borehole wall boundary temperature given by the
<code>port_wall</code> port.
</p>
<p>
The heat transfer in the borehole is computed using a convective heat transfer
coefficient that depends on the fluid velocity, a heat resistance between each
pair of pipes, and a heat resistance between the pipes and the borehole wall.
The heat capacity of the fluid and the heat capacity of the grout are taken
into account. The vertical heat flow is assumed to be zero. 
</p>
</html>", revisions="<html>
<ul>
<li>
July 2018, by Alex Laferri&egrave;re:<br/>
Following major changes to the structure of the IBPSA.Fluid.HeatExchangers.Ground package,
the documentation has been changed to reflect the new role of this model.
Additionally, this model now extends a partial borehole model.
</li>
<li>
July 2014, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</html>"));
end TwoUTubeLin;
