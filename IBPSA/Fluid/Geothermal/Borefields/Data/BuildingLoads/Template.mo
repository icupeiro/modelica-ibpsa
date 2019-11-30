within IBPSA.Fluid.Geothermal.Borefields.Data.BuildingLoads;
record Template "Template for external g-functions data records"
  extends Modelica.Icons.Record;
  parameter Modelica.SIunits.Time[365] hou
  "Exponential time series";
  parameter Modelica.SIunits.HeatFlowRate[365] Qbuih
  "Building heating needs";
  parameter Modelica.SIunits.HeatFlowRate[365] Qbuic
  "Building cooling needs";

  annotation (  defaultComponentName="gFunc",
  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Template;
