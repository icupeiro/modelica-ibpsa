within IBPSA.Fluid.Geothermal.Borefields.Data.Weather;
record Template "Template for external g-functions data records"
  extends Modelica.Icons.Record;
  parameter Modelica.SIunits.Time[8760] tim
  "Exponential time series";
  parameter Modelica.SIunits.HeatFlowRate[8760] Qbuih
  "Building heating needs";
  parameter Modelica.SIunits.HeatFlowRate[8760] Qbuic
  "Building cooling needs";

  annotation (  defaultComponentName="gFunc",
  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Template;
