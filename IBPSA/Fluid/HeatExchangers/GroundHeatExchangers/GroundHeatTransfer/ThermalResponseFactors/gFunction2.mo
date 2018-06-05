within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors;
function gFunction2 "Evaluate the g-function of a bore field"
  extends Modelica.Icons.Function;

  input Modelica.SIunits.Height hBor "Borehole length";
  input Modelica.SIunits.Height dBor "Borehole buried depth";
  input Modelica.SIunits.Radius rBor "Borehole radius";
  input Modelica.SIunits.Distance r
    "Distance from the borehole wall at which the G-function is evaluated";
  input Modelica.SIunits.ThermalDiffusivity alpha=1e-6
    "Ground thermal diffusivity used in g-function evaluation";
  input Integer nbTim=75 "Number of time steps";
  input Real ttsMax=exp(5) "Maximum adimensional time for gfunc calculation";

  output Real lntts[nbTim] "Logarithmic dimensionless time";
  output Real g[nbTim] "g-Function";

protected
  Modelica.SIunits.Time ts=hBor^2/(9*alpha) "Characteristic time";
  Modelica.SIunits.Time tVec_min=1 "Minimum time for short time calculations";
  Modelica.SIunits.Time tLon_max=ts*ttsMax
    "Maximum time for long time calculations";
  Modelica.SIunits.Time tVec[nbTim] "Time vector for short time calculations";
  Real FLS "Finite line source solution";
  Real ILS "Infinite line source solution";
  Real CHS "Cylindrical heat source solution";

algorithm

  // Generate geometrically expanding time vectors
  tVec :=
    IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.timeGeometric(
    tVec_min,
    tLon_max,
    nbTim) "Time vector for short time calculations";

  lntts := log(tVec/ts);

  // -----------------------
  // Short time calculations
  // -----------------------
  Modelica.Utilities.Streams.print(("Evaluation of short time g-function."));
  for k in 1:nbTim loop
    // Finite line source solution
    FLS :=
      IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.finiteLineSource(
      tVec[k],
      alpha,
      r,
      hBor,
      dBor,
      hBor,
      dBor);
    // Infinite line source solution
    ILS := 0.5*
      IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.infiniteLineSource(
      tVec[k],
      alpha,
      r);
    // Cylindrical heat source solution
    CHS := 2*Modelica.Constants.pi*
      IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.ThermalResponseFactors.cylindricalHeatSource(
      tVec[k],
      alpha,
      r,
      rBor);
    // Correct finite line source solution for cylindrical geometry
    g[k] := FLS + (CHS - ILS);
  end for;

end gFunction2;
