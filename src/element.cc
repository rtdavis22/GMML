// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Kyle Forrester

#include "gmml/internal/element.h"

#include <stdexcept>

#include "gmml/internal/stubs/common.h"

using std::string;

namespace gmml {

// --- Element Data -------------------------------------------------------- //
// This data is auto-generated from the elements.xml file in the Blue Obelisk
// Element Repository using the bodr.py script in the scripts directory.

namespace {

struct ElementData {
    const std::string symbol;
    const std::string name;
    double mass;  // in g/mol
    double exact_mass;  // in g/mol
    double ionization_energy;  // in eV
    double electron_affinity;  // in eV
    double electron_negativity;  // pauling scale
    double covalent_radius;  // in angstroms
    double van_der_waals_radius; // in angstroms
    double boiling_point;  // in kelvin
    double melting_point;  // in kelvin
};

const ElementData kElementData[] = {
    {"Xx", "Dummy", 0.0000, 0.00000, 0, 0, 0, 0.0, 0, 0, 0},
    {"H", "Hydrogen", 1.00794, 1.007825032, 13.5984, 0.75420375, 2.20, 0.37, 1.2, 20.28, 14.01},
    {"He", "Helium", 4.002602, 4.002603254, 24.5874, 0, 0, 0.32, 1.4, 4.216, 0.95},
    {"Li", "Lithium", 6.941, 7.01600455, 5.3917, 0.618049, 0.98, 1.34, 2.2, 1615, 453.7},
    {"Be", "Beryllium", 9.012182, 9.0121822, 9.3227, 0, 1.57, 0.90, 1.9, 3243, 1560},
    {"B", "Boron", 10.811, 11.0093054, 8.2980, 0.279723, 2.04, 0.82, 1.8, 4275, 2365},
    {"C", "Carbon", 12.0107, 12, 11.2603, 1.262118, 2.55, 0.77, 1.7, 5100, 3825},
    {"N", "Nitrogen", 14.0067, 14.003074, 14.5341, -0.07, 3.04, 0.75, 1.6, 77.344, 63.15},
    {"O", "Oxygen", 15.9994, 15.99491462, 13.6181, 1.4611120, 3.44, 0.73, 1.55, 90.188, 54.8},
    {"F", "Fluorine", 18.9984032, 18.99840322, 17.4228, 3.4011887, 3.98, 0.71, 1.5, 85, 53.55},
    {"Ne", "Neon", 20.1797, 19.99244018, 21.5645, 0, 0, 0.69, 1.54, 27.1, 24.55},
    {"Na", "Sodium", 22.98976928, 22.98976928, 5.1391, 0.547926, 0.93, 1.54, 2.4, 1156, 371},
    {"Mg", "Magnesium", 24.3050, 23.9850417, 7.6462, 0, 1.31, 1.30, 2.2, 1380, 922},
    {"Al", "Aluminium", 26.9815386, 26.98153863, 5.9858, 0.43283, 1.61, 1.18, 2.1, 2740, 933.5},
    {"Si", "Silicon", 28.0855, 27.97692653, 8.1517, 1.389521, 1.90, 1.11, 2.1, 2630, 1683},
    {"P", "Phosphorus", 30.973762, 30.97376163, 10.4867, 0.7465, 2.19, 1.06, 1.95, 553, 317.3},
    {"S", "Sulfur", 32.065, 31.972071, 10.3600, 2.0771029, 2.58, 1.02, 1.8, 717.82, 392.2},
    {"Cl", "Chlorine", 35.453, 34.96885268, 12.9676, 3.612724, 3.16, 0.99, 1.8, 239.18, 172.17},
    {"Ar", "Argon", 39.948, 39.96238312, 15.7596, 0, 0, 0.97, 1.88, 87.45, 83.95},
    {"K", "Potassium", 39.0983, 38.96370668, 4.3407, 0.501459, 0.82, 1.96, 2.8, 1033, 336.8},
    {"Ca", "Calcium", 40.078, 39.96259098, 6.1132, 0.02455, 1.00, 1.74, 2.4, 1757, 1112},
    {"Sc", "Scandium", 44.955912, 44.9559119, 6.5615, 0.188, 1.36, 1.44, 2.3, 3109, 1814},
    {"Ti", "Titanium", 47.867, 47.9479463, 6.8281, 0.084, 1.54, 1.36, 2.15, 3560, 1935},
    {"V", "Vanadium", 50.9415, 50.9439595, 6.7462, 0.525, 1.63, 1.25, 2.05, 3650, 2163},
    {"Cr", "Chromium", 51.9961, 51.9405075, 6.7665, 0.67584, 1.66, 1.27, 2.05, 2945, 2130},
    {"Mn", "Manganese", 54.938045, 54.9380451, 7.4340, 0, 1.55, 1.39, 2.05, 2235, 1518},
    {"Fe", "Iron", 55.845, 55.9349375, 7.9024, 0.151, 1.83, 1.25, 2.05, 3023, 1808},
    {"Co", "Cobalt", 58.933195, 58.933195, 7.8810, 0.6633, 1.88, 1.26, 2, 3143, 1768},
    {"Ni", "Nickel", 58.6934, 57.9353429, 7.6398, 1.15716, 1.91, 1.21, 2, 3005, 1726},
    {"Cu", "Copper", 63.546, 62.9295975, 7.7264, 1.23578, 1.90, 1.38, 2, 2840, 1356.6},
    {"Zn", "Zinc", 65.38, 63.9291422, 9.3942, 0, 1.65, 1.31, 2.1, 1180, 692.73},
    {"Ga", "Gallium", 69.723, 68.9255736, 5.9993, 0.41, 1.81, 1.26, 2.1, 2478, 302.92},
    {"Ge", "Germanium", 72.64, 73.9211778, 7.8994, 1.232712, 2.01, 1.22, 2.1, 3107, 1211.5},
    {"As", "Arsenic", 74.92160, 74.9215965, 9.7886, 0.814, 2.18, 1.19, 2.05, 876, 1090},
    {"Se", "Selenium", 78.96, 79.9165213, 9.7524, 2.02067, 2.55, 1.16, 1.9, 958, 494},
    {"Br", "Bromine", 79.904, 78.9183371, 11.8138, 3.3635880, 2.96, 1.14, 1.9, 331.85, 265.95},
    {"Kr", "Krypton", 83.798, 83.911507, 13.9996, 0, 3.00, 1.10, 2.02, 120.85, 116},
    {"Rb", "Rubidium", 85.4678, 84.91178974, 4.1771, 0.485916, 0.82, 2.11, 2.9, 961, 312.63},
    {"Sr", "Strontium", 87.62, 87.9056121, 5.6949, 0.05206, 0.95, 1.92, 2.55, 1655, 1042},
    {"Y", "Yttrium", 88.90585, 88.9058483, 6.2173, 0.307, 1.22, 1.62, 2.4, 3611, 1795},
    {"Zr", "Zirconium", 91.224, 89.9047044, 6.6339, 0.426, 1.33, 1.48, 2.3, 4682, 2128},
    {"Nb", "Niobium", 92.90638, 92.9063781, 6.7589, 0.893, 1.6, 1.37, 2.15, 5015, 2742},
    {"Mo", "Molybdenum", 95.96, 97.9054082, 7.0924, 0.7472, 2.16, 1.45, 2.1, 4912, 2896},
    {"Tc", "Technetium", 98, 97.907216, 7.28, 0.55, 1.9, 1.56, 2.05, 4538, 2477},
    {"Ru", "Ruthenium", 101.07, 101.9043493, 7.3605, 1.04638, 2.2, 1.26, 2.05, 4425, 2610},
    {"Rh", "Rhodium", 102.90550, 102.905504, 7.4589, 1.14289, 2.28, 1.35, 2, 3970, 2236},
    {"Pd", "Palladium", 106.42, 105.903486, 8.3369, 0.56214, 2.20, 1.31, 2.05, 3240, 1825},
    {"Ag", "Silver", 107.8682, 106.905097, 7.5762, 1.30447, 1.93, 1.53, 2.1, 2436, 1235.1},
    {"Cd", "Cadmium", 112.411, 113.9033585, 8.9938, 0, 1.69, 1.48, 2.2, 1040, 594.26},
    {"In", "Indium", 114.818, 114.903878, 5.7864, 0.404, 1.78, 1.44, 2.2, 2350, 429.78},
    {"Sn", "Tin", 118.710, 119.9021947, 7.3439, 1.112066, 1.96, 1.41, 2.25, 2876, 505.12},
    {"Sb", "Antimony", 121.760, 120.9038157, 8.6084, 1.047401, 2.05, 1.38, 2.2, 1860, 903.91},
    {"Te", "Tellurium", 127.60, 129.9062244, 9.0096, 1.970875, 2.1, 1.35, 2.1, 1261, 722.72},
    {"I", "Iodine", 126.90447, 126.904473, 10.4513, 3.059038, 2.66, 1.33, 2.1, 457.5, 386.7},
    {"Xe", "Xenon", 131.293, 131.9041535, 12.1298, 0, 2.6, 1.30, 2.16, 165.1, 161.39},
    {"Cs", "Caesium", 132.9054519, 132.9054519, 3.8939, 0.471626, 0.79, 2.25, 3, 944, 301.54},
    {"Ba", "Barium", 137.327, 137.9052472, 5.2117, 0.14462, 0.89, 1.98, 2.7, 2078, 1002},
    {"La", "Lanthanum", 138.90547, 138.9063533, 5.5769, 0.47, 1.10, 1.69, 2.5, 3737, 1191},
    {"Ce", "Cerium", 140.116, 139.9054387, 5.5387, 0.5, 1.12, 0, 2.48, 3715, 1071},
    {"Pr", "Praseodymium", 140.90765, 140.9076528, 5.473, 0.5, 1.13, 0, 2.47, 3785, 1204},
    {"Nd", "Neodymium", 144.242, 141.9077233, 5.5250, 0.5, 1.14, 0, 2.45, 3347, 1294},
    {"Pm", "Promethium", 145, 144.912749, 5.582, 0.5, 0, 0, 2.43, 3273, 1315},
    {"Sm", "Samarium", 150.36, 151.9197324, 5.6437, 0.5, 1.17, 0, 2.42, 2067, 1347},
    {"Eu", "Europium", 151.964, 152.9212303, 5.6704, 0.5, 0, 0, 2.4, 1800, 1095},
    {"Gd", "Gadolinium", 157.25, 157.9241039, 6.1498, 0.5, 1.20, 0, 2.38, 3545, 1585},
    {"Tb", "Terbium", 158.92535, 158.9253468, 5.8638, 0.5, 0, 0, 2.37, 3500, 1629},
    {"Dy", "Dysprosium", 162.500, 163.9291748, 5.9389, 0.5, 1.22, 0, 2.35, 2840, 1685},
    {"Ho", "Holmium", 164.93032, 164.9303221, 6.0215, 0.5, 1.23, 0, 2.33, 2968, 1747},
    {"Er", "Erbium", 167.259, 165.9302931, 6.1077, 0.5, 1.24, 0, 2.32, 3140, 1802},
    {"Tm", "Thulium", 168.93421, 168.9342133, 6.1843, 0.5, 1.25, 0, 2.3, 2223, 1818},
    {"Yb", "Ytterbium", 173.054, 173.9388621, 6.2542, 0.5, 0, 0, 2.28, 1469, 1092},
    {"Lu", "Lutetium", 174.9668, 174.9407718, 5.4259, 0.5, 1.27, 1.60, 2.27, 3668, 1936},
    {"Hf", "Hafnium", 178.49, 179.94655, 6.8251, 0, 1.3, 1.50, 2.25, 4875, 2504},
    {"Ta", "Tantalum", 180.94788, 180.9479958, 7.5496, 0.322, 1.5, 1.38, 2.2, 5730, 3293},
    {"W", "Tungsten", 183.84, 183.9509312, 7.8640, 0.815, 2.36, 1.46, 2.1, 5825, 3695},
    {"Re", "Rhenium", 186.207, 186.9557531, 7.8335, 0.15, 1.9, 1.59, 2.05, 5870, 3455},
    {"Os", "Osmium", 190.23, 191.9614807, 8.4382, 1.07780, 2.2, 1.28, 2, 5300, 3300},
    {"Ir", "Iridium", 192.217, 192.9629264, 8.9670, 1.56436, 2.20, 1.37, 2, 4700, 2720},
    {"Pt", "Platinum", 195.084, 194.9647911, 8.9588, 2.12510, 2.28, 1.28, 2.05, 4100, 2042.1},
    {"Au", "Gold", 196.966569, 196.9665687, 9.2255, 2.30861, 2.54, 1.44, 2.1, 3130, 1337.58},
    {"Hg", "Mercury", 200.59, 201.970643, 10.4375, 0, 2.00, 1.49, 2.05, 629.88, 234.31},
    {"Tl", "Thallium", 204.3833, 204.9744275, 6.1082, 0.377, 1.62, 1.48, 2.2, 1746, 577},
    {"Pb", "Lead", 207.2, 207.9766521, 7.4167, 0.364, 2.33, 1.47, 2.3, 2023, 600.65},
    {"Bi", "Bismuth", 208.98040, 208.9803987, 7.2855, 0.942363, 2.02, 1.46, 2.3, 1837, 544.59},
    {"Po", "Polonium", 209, 208.9824304, 8.414, 1.9, 2.0, 0, 2, 0, 527},
    {"At", "Astatine", 210, 209.987148, 0, 2.8, 2.2, 0, 2, 610, 575},
    {"Rn", "Radon", 222, 222.0175777, 10.7485, 0, 0, 1.45, 2, 211.4, 202},
    {"Fr", "Francium", 223, 223.0197359, 4.0727, 0, 0.7, 0, 2, 950, 300},
    {"Ra", "Radium", 226, 226.0254098, 5.2784, 0, 0.9, 0, 2, 1413, 973},
    {"Ac", "Actinium", 227, 227.0277521, 5.17, 0, 1.1, 0, 2, 3470, 1324},
    {"Th", "Thorium", 232.03806, 232.0380553, 6.3067, 0, 1.3, 0, 2.4, 5060, 2028},
    {"Pa", "Protactinium", 231.03588, 231.035884, 5.89, 0, 1.5, 0, 2, 4300, 1845},
    {"U", "Uranium", 238.02891, 238.0507882, 6.1941, 0, 1.38, 0, 2.3, 4407, 1408},
    {"Np", "Neptunium", 237, 237.0481734, 6.2657, 0, 1.36, 0, 2, 4175, 912},
    {"Pu", "Plutonium", 244, 244.064204, 6.0260, 0, 1.28, 0, 2, 3505, 913},
    {"Am", "Americium", 243, 243.0613811, 5.9738, 0, 1.3, 0, 2, 2880, 1449},
    {"Cm", "Curium", 247, 247.070354, 5.9914, 0, 1.3, 0, 2, 3383, 1620},
    {"Bk", "Berkelium", 247, 247.070307, 6.1979, 0, 1.3, 0, 2, 983, 1258},
    {"Cf", "Californium", 251, 251.079587, 6.2817, 0, 1.3, 0, 2, 1173, 1172},
    {"Es", "Einsteinium", 252, 252.08298, 6.42, 0, 1.3, 0, 2, 0, 1130},
    {"Fm", "Fermium", 257, 257.095105, 6.50, 0, 1.3, 0, 2, 0, 1800},
    {"Md", "Mendelevium", 258, 258.098431, 6.58, 0, 1.3, 0, 2, 0, 1100},
    {"No", "Nobelium", 259, 259.10103, 6.65, 0, 1.3, 0, 2, 0, 1100},
    {"Lr", "Lawrencium", 262, 262.10963, 4.9, 0, 0, 0, 2, 0, 1900},
    {"Rf", "Rutherfordium", 265, 261.10877, 6.0, 0, 0, 0, 2, 0, 0},
    {"Db", "Dubnium", 268, 262.11408, 0, 0, 0, 0, 2, 0, 0},
    {"Sg", "Seaborgium", 271, 263.11832, 0, 0, 0, 0, 2, 0, 0},
    {"Bh", "Bohrium", 272, 264.1246, 0, 0, 0, 0, 2, 0, 0},
    {"Hs", "Hassium", 270, 265.13009, 0, 0, 0, 0, 2, 0, 0},
    {"Mt", "Meitnerium", 276, 268.13873, 0, 0, 0, 0, 2, 0, 0},
    {"Ds", "Darmstadtium", 281, 271.14606, 0, 0, 0, 0, 0, 0, 0},
    {"Rg", "Roentgenium", 280, 272.15362, 0, 0, 0, 0, 0, 0, 0},
    {"Cn", "Copernicium", 285, 285.17411, 0, 0, 0, 0, 0, 0, 0},
    {"Uut", "Ununtrium", 284, 284.17808, 0, 0, 0, 0, 0, 0, 0},
    {"Uuq", "Ununquadium", 289, 289.18728, 0, 0, 0, 0, 0, 0, 0},
    {"Uup", "Ununpentium", 288, 288.19249, 0, 0, 0, 0, 0, 0, 0},
    {"Uuh", "Ununhexium", 293, 292.19979, 0, 0, 0, 0, 0, 0, 0},
    {"Uus", "Ununseptium", 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {"Uuo", "Ununoctium", 294, 0, 0, 0, 0, 0, 0, 0, 0},
};

const int kElementDataSize = ARRAY_SIZE(kElementData);

}  // namespace

Element::Element(const std::string& symbol) {
    atomic_number_ = 0;
    for (int i = 1; i < kElementDataSize; i++) {
        if (symbol == kElementData[i].symbol) {
            atomic_number_ = i;
            break;
        }
    }
    if (atomic_number_ == 0) {
        throw std::invalid_argument("Element symbol not recognized.");
    }
}

Element::Element(int atomic_number) : atomic_number_(atomic_number) {
    if (atomic_number_ < 1 || atomic_number_ >= kElementDataSize) {
       throw std::invalid_argument("Atomic Number not recognized.");
    }
}

Element Element::from_name(const std::string &name) {
    for (int i = 1; i < kElementDataSize; i++) {
        if (name == kElementData[i].name)
            return Element(i);
    }
    throw std::invalid_argument("Element Name not recognized.");
}

std::string Element::symbol() const {
    return kElementData[atomic_number_].symbol;
}

std::string Element::name() const {
    return kElementData[atomic_number_].name;
}

int Element::period() const {
    if (atomic_number_ < 3)
        return 1;
    else if (atomic_number_ < 11)
        return 2;
    else if (atomic_number_ < 19)
        return 3;
    else if (atomic_number_ < 37)
        return 4;
    else if (atomic_number_ < 55)
        return 5;
    else if (atomic_number_ < 87)
        return 6;
    else
        return 7;
}

double Element::mass() const {
    return kElementData[atomic_number_].mass;
}

double Element::exact_mass() const {
    return kElementData[atomic_number_].exact_mass;
}

double Element::ionization_energy() const {
    return kElementData[atomic_number_].ionization_energy;
}

double Element::electron_affinity() const {
    return kElementData[atomic_number_].electron_affinity;
}

double Element::electron_negativity() const {
    return kElementData[atomic_number_].electron_negativity;
}

double Element::covalent_radius() const {
    return kElementData[atomic_number_].covalent_radius;
}

double Element::van_der_waals_radius() const {
    return kElementData[atomic_number_].van_der_waals_radius;
}

double Element::boiling_point() const {
    return kElementData[atomic_number_].boiling_point;
}

double Element::melting_point() const {
    return kElementData[atomic_number_].melting_point;
}

}  // namespace gmml
