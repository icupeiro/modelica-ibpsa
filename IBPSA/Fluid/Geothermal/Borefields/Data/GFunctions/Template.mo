within IBPSA.Fluid.Geothermal.Borefields.Data.GFunctions;
record Template
  "Template for external g-functions data records"
  extends Modelica.Icons.Record;
  parameter Modelica.SIunits.Time[76] timExp
  "Exponential time series";
  parameter Real[76] gFunc
  "Adimensional values of the gFunction";

  annotation (  defaultComponentName="gFunc",
  Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Template;
