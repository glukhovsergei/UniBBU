

                                            UniBBU

UniBBU is a set of Python scripts for simulation of the regenerative multibunch beam breakup (BBU) instability in recirculating linacs and energy recovery linacs. It consists of the following configuration files and scripts:

machine.ini		Configuration file determining the machine to simulate
MESA.ini, S-DALINAC.ini	Configuration file containing the description of the particular machine
GenerateMachines.py	Generation of cavity mode configurations and lattices according to configuration files
CCP0.py			Stability analysis for monopole modes (approximated in case of multiple recirculations)
CCP1.py			Stability analysis for dipole modes (approximated in case of multiple recirculations)
CCP2.py			Stability analysis for quadrupole modes (approximated in case of multiple recirculations)
BunchTracking.py	Bunch tracking (exact)


Format of a configuration file containing the description of the particular machine:

[General]
CavityParameters	string		Path to the file with the parameters of the modes which can be excited in the cavity
FundamentalFrequency	float		Fundamental frequency, Hz
NumberOfModes		integer		Total number of modes which can be excited in the cavity
NumberOfModesByOrder	list of int	Number of monopole, dipole and quadrupole modes which can be excited in the cavity
LatticeFile		string		Path to the lattice file in .lte format (elegant)
Beamline		string		Name of the main beamline in the lattice file
CavityElement		string		Name of the beamline element representing the cavity in the lattice file
NumberOfCavities	int		Total number of cavities in the lattice
NumberOfCavityPasses	int		Total number of cavities passed by a bunch from its injection to ejection
p_central		float		Initial relativistic beta*gamma
beta_x			float		Initial horizontal Twiss beta-function, m
alpha_x			float		Initial vertical Twiss beta-function, m
beta_y			float		Initial horizontal Twiss alpha-function
alpha_y			float		Initial vertical Twiss alpha-function
emit_x			float		Initial horizontal emittance, m
emit_y			float		Initial vertical emittance, m
sigma_s			float		Initial bunch length, m
sigma_dp		float		Initial bunch relative energy spread

[Run]
UsedCavityPasses	int		Consider only UsedCavityPasses cavity passes in all simulations
Machine			int		Index of the machine (-1 for all machines)
Modes			list of int	Indices of active modes (-1 for all modes of the chosen order, not used in tracking)
UseSweep		bool		Use parameter sweep mode
NameX			string		Name of the first varied parameter in the lattice file
NameY			string		Name of the second varied parameter in the lattice file
StartX			float		Initial value of the first varied parameter
StartY			float		Initial value of the second varied parameter
StepX			float		Step of the first varied parameter
StepY			float		Step of the second varied parameter
NstepsX			int		Number of steps of the first varied parameter
NstepsY			int		Number of steps of the second varied parameter

[Spread]
NumberOfMachines	int		Number of generated machines
PreserveCavities	bool		Use previously generated cavity parameters
SameCavities		bool		Generate cavity imperfections only once
SameLattice		bool		Generate lattice imperfections only once
FrequencySpread		float		Mode frequency spread, Hz (-1 to use values from CavityParameters file)
DistrCutoffSigma	float		Number of sigmas where all imperfection distributions are cut off
RandomPolarization	bool		Random polarization of the first dipole and quadrupole mode in a pair
QuadGradient		float		Quadrupole focusing strength spread, 1/m
QuadRoll		float		Quadrupole roll spread, rad
IdealLatticeFirst	bool		Use lattice without imperfections in the first machine

[CCP]
RelativeRange		float		Relative frequency scan range
Step			float		Frequency scan step, in units of fundamental frequency
MaxCurrent		float		Maximal absolute current value to display in complex current plot, A
OutputFrequencyUnits	float		Frequency units in complex current plot (dimensionless factor)
OutputCurrentUnits	float		Current units in output files (dimensionless factor)

[Tracking]
Current			float		Beam current, A
NumberOfBunches		int		Number of bunches to track
InitialModeAmplitudes	list of float	Initial values of monopole, dipole and quadrupole modes used in tracking
MaxVoltLines		int		Maximal number of lines in the output file with mode voltage values
MaxCoordLines		int		Maximal number of lines in the output file with final bunch parameters
OutputTimeUnits		float		Time units in output files (dimensionless factor)
LongitudinalMotion	bool		Enable longitudinal motion (mandatory for tracking with monopole modes)


Generated machines data and simulation results are saved in the 'machines' folder with the following structure:

machines
  cavities
    xxxx.txt	Mode parameters for all cavities in the machine No. xxxx (including imperfections)
  lattices
    xxxx.txt	Lattice file in .lte format for the machine No. xxxx (including imperfections)
  matrices
    xxxx.txt	Relativistic beta*gamma values, lengths and transport matrices of all lattice regions between subsequent cavities for the machine No. xxxx (including imperfections)
  results
    Simulation results
  transition
    xxxx.txt	Transition coefficients for the machine No. xxxx (including imperfections)


Sergei Glukhov

2019-2022