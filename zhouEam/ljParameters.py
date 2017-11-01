# http://www.reflectometry.org/danse/elements.html
import periodictable as pt
import mendeleev as mdlv

from ase.units import _Nav, kB, kJ
JOULE_TO_EV = kJ / 1000
AVOGADRO = _Nav
EV_AT = JOULE_TO_EV * 1000 / AVOGADRO

class LjParameters():

    def thermal_props(self, element, EV=True):

        element = mdlv.element(element)
        mp = element.melting_point
        bp = element.boiling_point

        if EV == True:
            ev_h = element.evaporation_heat * EV_AT
            h_o_f = element.heat_of_formation * EV_AT
            f_h = element.fusion_heat * EV_AT
        else:
            ev_h = element.evaporation_heat
            h_o_f = element.heat_of_formation
            f_h = element.fusion_heat

        data = {'evaporation_heat': ev_h,
            'heat_of_formation': h_o_f,
            'fusion_heat': f_h,
            'melting_point': mp,
            'boiling_point':bp
           }
        return data

    '''
       REFS:
       Calculation of the parameters of the Lennard-Jones potential for pairs of
       identical atoms based on the properties of solid substances

       V. P. Filippova, S. A. Kunavin, M. S. Pugachev


    '''
    def calc_lj_01 (self, element):

        form = pt.formula(element)
        e_ = form.structure[0][1]

        data = self.thermal_props(element)

        if  e_.crystal_structure['symmetry'] == 'fcc' or \
            e_.crystal_structure['symmetry'] == 'hcp':

            sigma = e_.crystal_structure['a'] * 0.635
            epsilon = data['evaporation_heat'] * 0.172

        elif e_.crystal_structure['symmetry'] == 'BCC':
            sigma = e_.crystal_structure['a'] *0.788
            epsilon = data['evaporation_heat'] * 0.1819

        else:
            sigma = e_.crystal_structure['a'] * 0.7  # 0.648
            epsilon = data['evaporation_heat'] * 0.18
            #epsilon = data['fusion_heat']  * 0.843

        return ({'epsilon': epsilon, 'sigma':sigma})

    def lammpsInteraction(self, atoms):
        self.atoms = atoms

        str_ = 'pair_style lj/cut 7\n'
     
        for i,e in enumerate(self.atoms):
            el = e['ase'].symbol
            paramsLJ = self.calc_lj_01(el)
            # eps sigma cut
            coefPair = '  ' + str(paramsLJ['epsilon']) + ' ' +\
                    str(paramsLJ['sigma']) + '  ' + \
                    str(paramsLJ['sigma'] * 1.5) +  '\n'
            j = i + 1
            str_ += 'pair_coeff ' + str(j) + ' ' + str(j) + coefPair 

        return str_



def test_lj():
    ljp = LjParameters()

    elements = ['U', 'Fe', 'Mo','V', 'Cr','Cu', 'Ni', 'Al',\
                'Zr', 'Si', 'Mg', 'Co', 'Ag', 'Pt', 'W', \
                'Pb','Au', 'Pd', 'Ti', 'Nb', 'Ta' ]

    for e in elements:
        print (e, ljp.calc_lj_01(e))

if __name__ == '__main__':
    test_lj()
