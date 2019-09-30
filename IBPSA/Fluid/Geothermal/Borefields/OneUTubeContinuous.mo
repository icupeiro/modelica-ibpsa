within IBPSA.Fluid.Geothermal.Borefields;
model OneUTubeContinuous
  extends OneUTube(redeclare
      BaseClasses.HeatTransfer.GroundTemperatureResponse_ContinuousRecord
      groTemRes);
end OneUTubeContinuous;
