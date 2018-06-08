within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors;
function shaGFunction2
  "Return a pseudo sha code of the formatted arguments for the g-function generation"
  extends Modelica.Icons.Function;
  input Real hBor "Borehole length";
  input Real dBor "Borehole buried depth";
  input Real rBor "Borehole radius";
  input Real alpha "Ground thermal diffusivity used in g-function evaluation";
  input Real r "Distance from the borehole wall at which the soil temperature is evaluated";

  output String sha
  "Pseudo sha code of the g-function arguments";

protected
  String shaStr =  "";
  String formatStr =  "1.3e";

algorithm
  shaStr := shaStr + String(hBor, format=formatStr);
  shaStr := shaStr + String(dBor, format=formatStr);
  shaStr := shaStr + String(rBor, format=formatStr);
  shaStr := shaStr + String(alpha, format=formatStr);
  shaStr := shaStr + String(r, format=formatStr);

  sha := IBPSA.Utilities.Cryptographics.sha(shaStr);
end shaGFunction2;
