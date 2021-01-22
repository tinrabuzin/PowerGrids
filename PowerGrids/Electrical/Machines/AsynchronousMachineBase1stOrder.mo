within PowerGrids.Electrical.Machines;

model AsynchronousMachineBase1stOrder
  extends AsynchronousMachineBase;
initial equation
//  0 = -RrPu / LrrPu * (udpPuStart + 1 * LmPu^2/LrrPu*iqPuStart) + slipPuStart * uqpPuStart;
//  0 = -RrPu / LrrPu * (uqpPuStart - 1 * LmPu^2/LrrPu*idPuStart) - slipPuStart * udpPuStart;
equation
  // Ignoring rotor dynamics
  0 = -RrPu / LrrPu * (udpPu + 1 * LmPu^2/LrrPu*iqPu) + slipPu * uqpPu;
  0 = -RrPu / LrrPu * (uqpPu - 1 * LmPu^2/LrrPu*idPu) - slipPu * udpPu;
end AsynchronousMachineBase1stOrder;
