# Microgrid Load Flow Analysis

A comprehensive MATLAB-based simulation tool for analyzing load flow in DC and AC microgrids operating in islanded mode. This project implements advanced numerical methods to solve power flow equations and evaluate steady-state behavior of microgrid systems.

## Overview

This project provides an in-depth analysis of load flow in DC and AC microgrids using MATLAB simulations. It examines the steady-state behavior of power systems, focusing on voltage magnitude, phase angles, current distribution, and power flows across different bus configurations.

### Key Features

- **DC Microgrid Analysis**: Newton-Raphson method for fast convergence
- **AC Microgrid Analysis**: Modified Gauss-Seidel iterative method
- **Multiple Test Systems**: 6-bus and 38-bus microgrid configurations
- **PV Bus Support**: Handles voltage-controlled buses with reactive power limits
- **Voltage-Dependent Loads**: Supports ZIP load models with frequency dependency
- **Droop Control**: Implements P-f and Q-V droop characteristics for DG units
- **Comprehensive Visualization**: Automatic generation of convergence and performance plots

## Mathematical Methods

### DC Microgrid (Newton-Raphson)

The DC microgrid analysis uses the Newton-Raphson method to solve the following equations:

- **Virtual Resistance Droop**: V₀ᵢ - Vᵢ - Rᵥᵢ(PGᵢ/Vᵢ) = 0
- **Power Balance**: PGᵢ - PLᵢ - Gᵢᵢ·Vᵢ² - Vᵢ·∑(Gᵢⱼ·Vⱼ) = 0

The Jacobian matrix is computed analytically for fast convergence (typically 3-5 iterations).

### AC Microgrid (Modified Gauss-Seidel)

The AC microgrid analysis implements a modified Gauss-Seidel method with:

- **P-f Droop Control**: PGᵢ = (1/mₚᵢ)(ω₀ - ω)
- **Q-V Droop Control**: QGᵢ = (1/nqᵢ)(V₀ᵢ - |Vᵢ|)
- **Frequency Update**: ω = (∑(1/mₚᵢ)·ω₀ - (PL_total + Ploss)) / ∑(1/mₚᵢ)

The method handles both PQ and PV buses with reactive power limit enforcement.

## Project Structure

```
.
├── README.md                               # This file
├── MAIN_PRG.m                             # AC microgrid - constant load
├── MAIN_PRG_LOAD.m                        # AC microgrid - voltage-dependent load
├── MAIN_PRG_DC.m                          # DC microgrid analysis
├── MAIN_MG_AC_LF_Gauss_Seidel_FINAL_PV.m # AC microgrid with PV bus support
├── Ybus_matrix.m                          # Y-bus matrix formation
├── Gbus_matrix.m                          # G-bus matrix formation (DC)
├── Jacobian_matrix.m                      # Jacobian calculation (DC)
├── system_data_THESIS.m                   # AC system data (6/38 bus)
├── system_data_WITH_PV.m                  # AC system data with PV buses
└── system_data_DC.m                       # DC system data
```

## Test Systems

### 6-Bus AC Microgrid

- **Buses**: 3 PQ load buses, 3 DG buses
- **Total Load**: ~4.2 kW active, ~2.1 kVAR reactive
- **DG Capacity**: 2.83 kW per unit (max)
- **Base Values**: Sbase = 1000 VA, Vbase = 220/√3 V
- **Features**: Configurable PV buses with Q-limit enforcement

### 38-Bus AC Microgrid

- **Buses**: 33 load buses, 5 DG buses
- **Total Load**: ~4.0 pu active, ~2.5 pu reactive
- **DG Capacity**: Variable (0.5-3.0 pu per unit)
- **Topology**: IEEE 33-bus derived with DG integration
- **Features**: Complex meshed network with multiple feeders

### 6-Bus DC Microgrid

- **Buses**: 3 load buses, 3 DG buses
- **Control**: Virtual resistance droop
- **Load**: Voltage-dependent (ZIP model)
- **Base Values**: Pbase = 5000 W, Vbase = 380 V

## Getting Started

### Prerequisites

- MATLAB R2016b or later
- No additional toolboxes required

### Running the Simulations

#### AC Microgrid (Constant Load)

```matlab
% Run 6-bus system
MG_system = 1;
MAIN_PRG

% Run 38-bus system
MG_system = 2;
MAIN_PRG
```

#### AC Microgrid (Voltage-Dependent Load)

```matlab
% Run with ZIP load model
MAIN_PRG_LOAD
```

#### AC Microgrid (With PV Buses)

```matlab
% Run with voltage-controlled buses
MAIN_MG_AC_LF_Gauss_Seidel_FINAL_PV
```

#### DC Microgrid

```matlab
% Run DC analysis
MAIN_PRG_DC
```

### Configuration

Modify system parameters in the respective `system_data_*.m` files:

```matlab
% In system_data_THESIS.m
MG_system = 1;  % 1 for 6-bus, 2 for 38-bus

% Adjust droop coefficients
mp = [10.4e-5; 9.4e-5; 8.4e-5];  % P-f droop (rad/s/W)
nq = [0.0013; 0.0013; 0.0013];    % Q-V droop (V/VAR)

% Modify load parameters
alpha = 2;  % Voltage exponent for active power
beta = 2;   % Voltage exponent for reactive power
```

## Output

The simulation generates six figures for each run:

1. **Convergence Plot**: Log₁₀(error) vs iterations
2. **Voltage Profile**: Bus voltages vs iterations
3. **Active Power**: Generator outputs vs iterations
4. **Reactive Power**: Generator outputs vs iterations
5. **Frequency**: System frequency vs iterations
6. **System Losses**: P and Q losses vs iterations

### Sample Results (6-Bus AC System)

- **Convergence**: Achieved in ~15-25 iterations (ε < 10⁻⁸)
- **Voltage Range**: 0.95-1.01 pu
- **Frequency**: ~0.9998-1.0 pu
- **Total Losses**: ~2-3% of total generation

## Key Algorithms

### Y-Bus Matrix Formation

```matlab
function [Ybus] = Ybus_matrix(LF,LT,y_line,nbus,nline)
    Ybus = zeros(nbus);
    for i = 1:nline
        Ybus(LF(i),LF(i)) = Ybus(LF(i),LF(i)) + y_line(i);
        Ybus(LT(i),LT(i)) = Ybus(LT(i),LT(i)) + y_line(i);
        Ybus(LF(i),LT(i)) = Ybus(LF(i),LT(i)) - y_line(i);
        Ybus(LT(i),LF(i)) = Ybus(LF(i),LT(i));
    end
end
```

### Voltage Update (Gauss-Seidel)

```matlab
V(i) = (1/Ybus(i,i)) * (conj(Sinj(i)/V(i)) - sum_j(Ybus(i,j)*V(j)))
```

### PV Bus Handling

For PV buses, the algorithm:
1. Computes reactive power from current voltage
2. Checks against Qmax and Qmin limits
3. Enforces voltage magnitude if within limits
4. Treats as PQ bus if limits violated

## Validation

Results have been validated against:
- IEEE test case benchmarks
- Published research papers (Singh et al., EEEIC 2015)
- Commercial power system analysis tools

Typical accuracy: Voltage < 0.01%, Power < 0.1%

## Applications

- Microgrid planning and design
- DG integration studies
- Voltage stability analysis
- Droop control parameter tuning
- Load flow education and research

## Limitations

- Assumes balanced three-phase operation
- Islanded mode only (no grid connection)
- Static analysis (no dynamics)
- No optimal power flow features

## Future Enhancements

- [ ] Unbalanced three-phase analysis
- [ ] Dynamic simulation capability
- [ ] Optimal power flow (OPF)
- [ ] Renewable energy source modeling
- [ ] Battery storage integration
- [ ] Grid-connected mode operation

## References

1. "A Simple and Accurate Approach to Solve the Power Flow for Balanced Islanded Microgrids", IEEE EEEIC 2015
2. Singh et al., "Novel Approach for Load Flow Analysis of Islanded Microgrids"
3. IEEE 33-bus test system documentation

## Author

Developed during pre-final year under the guidance of:
- **Professor Vinay Pant**
- Electrical Engineering Department
- Indian Institute of Technology Roorkee

## License

This project is intended for academic and research purposes.

## Acknowledgments

Special thanks to Prof. Vinay Pant for guidance and the Electrical Engineering Department at IIT Roorkee for providing the research environment.

---

For questions or collaboration, please open an issue in this repository.
