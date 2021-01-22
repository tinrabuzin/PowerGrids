within PowerGrids.Electrical.Machines;

model AsynchronousMachineBase
  import PowerGrids.Types.Choices.LocalInitializationOption;
  import PowerGrids.Types.Choices.InitializationOption;
  extends Icons.Machine;
  extends BaseClasses.OnePortACdqPU(portVariablesPhases=true,
  localInit = if initOpt == InitializationOption.localSteadyStateFixedPowerFlow
                then LocalInitializationOption.PQ
                else LocalInitializationOption.none);

  parameter Types.PerUnit RsPu "Stator resistance in p.u.";
  parameter Types.PerUnit LsPu "Stator leakage inductance  in p.u.";
  parameter Types.PerUnit RrPu "Rotor resistance  in p.u.";
  parameter Types.PerUnit LrPu "Rotor leakage inductance  in p.u.";
  parameter Types.PerUnit LmPu "Magnetization inductance  in p.u.";
  parameter Modelica.SIunits.Time H "Kinetic constant = kinetic energy / rated power";
  parameter Real ExpCm "Exponent of the torque-speed dependance - omegaPu^ExpCm";
  parameter Types.Choices.InitializationOption initOpt = systemPowerGrids.initOpt "Initialization option";
  parameter Boolean useLoadInput = false;
  
  parameter Types.PerUnit LssPu = LsPu + LmPu;
  parameter Types.PerUnit LrrPu = LrPu + LmPu;
  
  constant Types.PerUnit omegaNomPu = 1 "Nominal frequency in p.u.";

  final parameter Types.PerUnit udpPuStart(fixed = false) "";
  final parameter Types.PerUnit uqpPuStart(fixed = false) "";
  final parameter Types.PerUnit slipPuStart(fixed = false) "Start value of phase-to-phase voltage phasor, phase angle";
    
  final parameter SI.AngularVelocity omegaBase = systemPowerGrids.omegaNom "Base angular frequency value";
  Types.PerUnit omegaRefPu = systemPowerGrids.omegaRef/systemPowerGrids.omegaNom "Reference frequency in p.u.";
  
  // Mechanical variables
  Types.PerUnit slipPu(start=0.000499) "Machine slip in p.u.";
  Types.PerUnit omegaPu(start=1-0.000499) "Angular frequency in p.u.";
  Types.PerUnit CePu(start=port.PStart/port.SNom) "Electrical torque in p.u. (base SNom/omegaBase)";
  Types.PerUnit CmPuScaled(start=port.PStart/port.SNom) "Mechanical torque in p.u. (base PNom/omegaBase)";
  Types.PerUnit CmPu(start=port.PStart/port.SNom) "The value of the mechanical torque when omegaPu = 1";
  
  Types.PerUnit udpPu(start=udpPuStart);
  Types.PerUnit uqpPu(start=uqpPuStart);
  
initial equation
  // The equation to compute the starting slip from the equivalent circuit
  port.QStart/port.SNom = (port.UStart/port.UNom)^2*(1 / LmPu + (LsPu + LrPu) 
                         / CM.'abs'(Complex(RsPu + RrPu / slipPuStart, LsPu + LrPu))^2); 

  // Start equations of the d-q voltages
  udPuStart -  udpPuStart = RsPu * idPuStart - 1*(LssPu - LmPu^2 / LrrPu) * iqPuStart;
  uqPuStart -  uqpPuStart = RsPu * iqPuStart + 1*(LssPu - LmPu^2 / LrrPu) * idPuStart;
  
  // Start equation for load and electrical torque
  //port.PStart/port.SNom = udPuStart * idPuStart + uqPuStart * iqPuStart;
  
  // Equations to determine the initial state values
  if initOpt == InitializationOption.noInit then
    // No initial equations
  else
    der(slipPu) = 0;
    CePu = CmPuScaled;
  end if;
equation
  // The equations will be written in the stator-aligned synchronous rotating frame
  theta = port.UPhase;
  
  // Mechanical equations
  der(CmPu) = 0;    
  CmPuScaled = CmPu*omegaPu^ExpCm;
  CePu = udpPu * idPu + uqpPu * iqPu;
  slipPu = omegaRefPu - omegaPu;
  der(omegaPu) = 1 / (2 * H) * (CePu - CmPuScaled);
  
  // Stator voltages
  udPu - udpPu = RsPu * idPu - 1 * (LssPu - LmPu^2 / LrrPu) * iqPu;
  uqPu - uqpPu = RsPu * iqPu + 1 * (LssPu - LmPu^2 / LrrPu) * idPu;
  
end AsynchronousMachineBase;
