within IBPSA.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer;
function interpolate "Interpolate linearly in a vector"
  extends Modelica.Icons.Function;
  input Real x[76]
    "Abscissa table vector (strict monotonically increasing values required)";
  input Real y[76] "Ordinate table vector";
  input Real xi "Desired abscissa value";
  input Integer iLast=1 "Index used in last search";
  output Real yi "Ordinate value corresponding to xi";
  output Integer iNew=1 "xi is in the interval x[iNew] <= xi < x[iNew+1]";
protected
  Integer i;
  Integer nx=76;
  Real x1;
  Real x2;
  Real y1;
  Real y2;
algorithm
  //assert(nx > 0, "The table vectors must have at least 1 entry.");
  if nx == 1 then
    yi := y[1];
  else
    // Search interval
    i := min(max(iLast, 1), nx - 1);
    if xi >= x[i] then
      // search forward
      while i < nx and xi >= x[i] loop
        i := i + 1;
      end while;
      i := i - 1;
    else
      // search backward
      while i > 1 and xi < x[i] loop
        i := i - 1;
      end while;
    end if;

    // Get interpolation data
    x1 := x[i];
    x2 := x[i + 1];
    y1 := y[i];
    y2 := y[i + 1];

    //assert(x2 > x1, "Abscissa table vector values must be increasing");
    // Interpolate
    yi := y1 + (y2 - y1)*(xi - x1)/(x2 - x1);
    iNew := i;
  end if;

  annotation (smoothOrder( normallyConstant=x, normallyConstant=y)=100,
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
// Real    x[:], y[:], xi, yi;
// Integer iLast, iNew;
        yi = Vectors.<strong>interpolate</strong>(x,y,xi);
(yi, iNew) = Vectors.<strong>interpolate</strong>(x,y,xi,iLast=1);
</pre></blockquote>
<h4>Description</h4>
<p>
The function call \"<code>Vectors.interpolate(x,y,xi)</code>\" interpolates
<strong>linearly</strong> in vectors
(x,y) and returns the value yi that corresponds to xi. Vector x[:] must consist
of monotonically increasing values. If xi &lt; x[1] or &gt; x[end], then
extrapolation takes places through the first or last two x[:] values, respectively.
If the x and y vectors have length 1, then always y[1] is returned.
The search for the interval x[iNew] &le; xi &lt; x[iNew+1] starts at the optional
input argument \"iLast\". The index \"iNew\" is returned as output argument.
The usage of \"iLast\" and \"iNew\" is useful to increase the efficiency of the call,
if many interpolations take place.
If x has two or more identical values then interpolation utilizes the x-value
with the largest index.
</p>

<h4>Example</h4>

<blockquote><pre>
  Real x1[:] = { 0,  2,  4,  6,  8, 10};
  Real x2[:] = { 1,  2,  3,  3,  4,  5};
  Real y[:]  = {10, 20, 30, 40, 50, 60};
<strong>algorithm</strong>
  (yi, iNew) := Vectors.interpolate(x1,y,5);  // yi = 35, iNew=3
  (yi, iNew) := Vectors.interpolate(x2,y,4);  // yi = 50, iNew=5
  (yi, iNew) := Vectors.interpolate(x2,y,3);  // yi = 40, iNew=4
</pre></blockquote>
</html>"));
end interpolate;
