within IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.LoadAggregation.Validation;
model CountAggregationCells "This validation case verifies the counting of the required length of aggregation vectors"
  extends Modelica.Icons.Example;

  Integer i;

equation
  i = IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.LoadAggregation.countAggregationCells(
      lvlBas=2,
      nCel=2,
      timFin=120,
      tLoaAgg=10);

annotation (
__Dymola_Commands(file="modelica://IBPSA/Resources/Scripts/Dymola/Fluid/HeatExchangers/GroundHeatExchangers/GroundHeatTransfer/LoadAggregation/Validation/CountAggregationCells.mos"
        "Simulate and plot"),
Documentation(info="<html>
<p>
This validation case counts the required length of the aggregation vectors for the
same fictional case as in
<a href=\"modelica://IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.LoadAggregation.Validation.AggregationCellTimes\">
IBPSA.Fluid.HeatExchangers.GroundHeatExchangers.GroundHeatTransfer.LoadAggregation.Validation.AggregationCellTimes</a>.
</p>
</html>", revisions="<html>
<ul>
<li>
July 18, 2018, by Alex Laferri&egrave;re:<br/>
First implementation.
</li>
</ul>
</html>"));
end CountAggregationCells;