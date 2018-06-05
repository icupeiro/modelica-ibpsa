within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors;
model finiteCylindricalHeatSource
  extends Modelica.Icons.Function;

  input Real t "Time";
  input Real r "Radial distance from line source";
  input Real alpha "Ground thermal diffusivity";
  input Real len "Length of borehole";
  input Real burDep "Buried depth of emitting borehole";

  output Real h_finCylSou "Thermal response factor of finite line source at a radial distance r from the source";
protected
  Real lowBou = 1.0/sqrt(4*alpha*t) "Lower bound of integration";
  // Upper bound is infinite
  Real h_finLinSou "Thermal response factor of finite line source at a radial distance r from the source";
  Real h_cylLinSou "Thermal response factor of cylindrical source at a radial distance r from the source";
  Real h_infLinSou "Thermal response factor of infinite line source at a radial distance r from the source";
algorithm
  h_finLinSou := IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.finiteLineSource(
      lowBou=lowBou,
      dis=r,
      len1=len,
      burDep1=burDep,
      len2=len,
      burDep2=burDep,
      includeRealSource=true,
      includeMirrorSource=true);
  h_cylLinSou :=   IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.finiteLineSource(
      t=t,
      alpha=alpha,
      dis=r,
      rBor=rBor);
  h_infLinSou :=   IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.finiteLineSource(
      t=t,
      alpha=alpha,
      dis=r);

   h_finCylSou:=  h_finLinSou + h_cylLinSou - h_infLinSou;
end finiteCylindricalHeatSource;
