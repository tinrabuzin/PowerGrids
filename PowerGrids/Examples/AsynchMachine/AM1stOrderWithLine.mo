within PowerGrids.Examples.AsynchMachine;

model AM1stOrderWithLine
  extends Modelica.Icons.Example;
  inner PowerGrids.Electrical.System systemPowerGrids(initOpt = PowerGrids.Types.Choices.InitializationOption.localSteadyStateFixedPowerFlow) annotation(
    Placement(visible = true, transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  PowerGrids.Electrical.Machines.AsynchronousMachine1stOrder am1(H = 2, PStart(displayUnit = "MW") = 100000, QStart(displayUnit = "MVA") = 43568.36380106667, RrPu = 0.00489, RsPu = 0.028, SNom(displayUnit = "V.A") = 0.1e6, UNom(displayUnit = "V") = 20000, UPhaseStart =-0.04925915594457064, UStart(displayUnit = "kV") = 19424.17645210658, XmPu = 4.236, XrPu = 0.1822, XsPu = 0.01) annotation(
    Placement(visible = true, transformation(origin = {0, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
Electrical.Branches.LineConstantImpedance line(R = 20, SNom(displayUnit = "W") = 1e6, UNom (displayUnit = "V") = 20000, UNomA(displayUnit = "V") = 20e3, X = 200) annotation(
    Placement(visible = true, transformation(origin = {0, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
Electrical.Buses.Bus bus(SNom(displayUnit = "V.A") = 0.1e6, UNom(displayUnit = "V") = 20e3)  annotation(
    Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Electrical.PowerFlow.SlackBus slack(SNom(displayUnit = "MVA") = 1000000, UNom(displayUnit = "kV") = 20000) annotation(
    Placement(visible = true, transformation(origin = {0, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(bus.terminal, am1.terminal) annotation(
    Line(points = {{0, -30}, {0, -50}}));
  connect(line.terminalB, bus.terminal) annotation(
    Line(points = {{0, -20}, {0, -30}}));
  connect(slack.terminal, line.terminalA) annotation(
    Line(points = {{0, 20}, {0, 0}}));
  annotation(
    Documentation(info = "<html><head></head><body>Infinite bus connected to an asynchronous machine.</body></html>"));
end AM1stOrderWithLine;
