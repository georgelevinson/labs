from asyncio import constants
from cmath import pi
import math
import matplotlib.pyplot as plt
from scipy import constants

class Experiment:
    exp_no = 0
    data_points = []
    temperature_coef = 0
    base_resistivity_OhmPerM = 0
    thermal_conductivity_coef = 0

    def __init__(self, data):
        self.data_points = data
        self.exp_no = self.data_points[0].test_no
        self.base_resistivity = self.data_points[0].resistivity_microOhmPerMeter * pow(10, -6)
        self.calc_TKp()
        self.calc_th_cond_coef()

    def calc_TKp(self):
        first = self.data_points[0]
        last = self.data_points[-1]
        self.temperature_coef = (last.resistivity_microOhmPerMeter / first.resistivity_microOhmPerMeter - 1)/(last.temperature_C - first.temperature_C)

    def calc_th_cond_coef(self):
        self.thermal_conductivity_coef = pow(self.base_resistivity, -1) * (pow(math.pi, 2) / 3) * pow(constants.k / constants.e, 2) * (20 + 273.15)

    def plot_to_file(self, plt):
        resistivities = [dp.resistivity_microOhmPerMeter for dp in self.data_points]
        temperatures = [dp.temperature_C for dp in self.data_points]
        plt.plot(resistivities, temperatures)
        plt.xlabel("Питомий опір (мкОм*м)")
        plt.ylabel("Температура ($^\circ$C)")
        plt.title("Зразок №{}".format(self.exp_no))
        plt.savefig("./exp{}.png".format(self.exp_no))
        plt.clf()

    def __str__(self) -> str:
        result = "\n==============\nEXPERIMENT #{}\n==============\nTKp = {}\nThermal conductivity coefficient = {}\n\n".format(self.exp_no, self.temperature_coef, self.thermal_conductivity_coef)

        for dp in self.data_points:
            result += str(dp)
        
        return result
            
class DataPoint:
    test_no = 0
    resistivity_microOhmPerMeter = 0
    resistance_ohm = 0
    temperature_C = 0
    length_m = 0
    diameter_m = 0
    cross_section_mSq = 0
    
    def __init__ (self, R, T, D, L, test_no):
        self.test_no = test_no
        self.resistance_ohm = R
        self.temperature_C = T
        self.length_m = L
        self.diameter_m = D * pow(10, -3)
        self.calc_cross_section()
        self.calc_resistivity()

    def calc_cross_section(self):
        self.cross_section_mSq = math.pi * pow(self.diameter_m, 2) / 4
    
    def calc_resistivity(self):
        self.resistivity_microOhmPerMeter = (self.resistance_ohm * self.cross_section_mSq / self.length_m) * pow(10,6)

    def __str__(self) -> str:
        return   """\nresistivity_microOhmPerMeter = {resistivity} temperature_C = {temperature} resistance_ohm = {resistance}\nlength_m = {length} diameter_m = {diameter} cross_section_mSq = {xsection}\n\n""".format(
                        resistivity = self.resistivity_microOhmPerMeter,
                        resistance = self.resistance_ohm, temperature = self.temperature_C,
                        length = self.length_m, diameter = self.diameter_m, xsection = self.cross_section_mSq)
    
# geom.parameters for all experiment subjects, (diameter in mm, length in meters)
geom_params = [(0.2, 7.0), (0.24, 3.5), (0.22, 2.7), (1.0, 15.0), (0.24, 3.0)]
# arrays of resistance values for every conducted exp, temperature is given in degrees Celsius by (20 + index * 10)
resistances = [
                [3.745, 3.907, 4.070, 4.232, 4.394, 4.556, 4.718, 4.880, 5.043, 5.205, 5.367],
                [77.406, 77.522, 77.638, 77.755, 77.871, 77.987, 78.103, 78.219, 78.335, 78.451, 78.567],
                [31.979, 31.984, 31.990, 31.995, 32.000, 32.006, 32.011, 32.017, 32.022, 32.028, 32.033],
                [0.506, 0.527, 0.548, 0.569, 0.589, 0.610, 0.631, 0.652, 0.672, 0.693, 0.714],
                [6.436, 6.838, 7.241, 7.643, 8.045, 8.447, 8.850, 9.252, 9.654, 10.056, 10.459]
              ]

data_points = []

for test_index, (dim, res) in enumerate(zip(geom_params, resistances)):
    for i, val in enumerate(res):
        dp = DataPoint(val, 20 + 10 * i, dim[0], dim[1], test_index + 1)
        data_points.append(dp)

experiments = []

for index in range(0,5):
    e = Experiment(data_points[index * 11 : (index + 1) * 11])
    experiments.append(e)

# region Print Output To Console

for e in experiments:
    print(e)

#endregion

# region Write Output To File

data_output_file = open('mrea_lab1.txt', 'w+')

for e in experiments:
    data_output_file.write(str(e))

#endregion

# region Write Plot Output to Files

for e in experiments:
    e.plot_to_file(plt)

#endregion

print("EXECUTION COMPLETED")
