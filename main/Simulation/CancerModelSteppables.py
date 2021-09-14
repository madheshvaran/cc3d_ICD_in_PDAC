
from cc3d.core.PySteppables import *
import numpy as np
from pathlib import Path

max_cancer_cell_volume = 100
start_cancer_cell_volume = 30
cancer_cell_growth_rate = 0.01
cancer_cell_mitosis_prob = 0.0045
apoptosis_prob = 0.02
apoptosis_shrink_rate = -0.5
necrosis_prob = 0.01
necrosis_shrink_rate = -0.8
cancer_min_oxygen_growth = 30
apoptosis_min_oxgen_level = 10
il2_recruiting_threshold = 1
tcell_born_prob = 0.0005
tcell_max_volume = 25
tcell_lambda_volume = 10.0
tcell_start_width = 5
il2_secretion = 1.
apoptosis_by_tcell_prob = 0.005
treg_born_prob = 0.0005
tgfb_secretion = 1.
no_coffee_prob = 0.005
no_coffee_shrink_rate = -0.5
treg_min_oxigen_level = 20
tcell_min_oxigen_level = 20

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):



        self.plot_win = self.add_new_plot_window(title='Tumor Volume',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)

        self.plot_win.add_plot("TumorVol", style='Dots', color='red', size=2)


       
    def step(self,mcs):
        tumorVol = np.sum([cell.volume for cell in self.cell_list_by_type(self.CANCER,self.CANCERAPOPTOSIS,self.CANCERNECROSIS)]) 
        self.plot_win.add_data_point("TumorVol", mcs, tumorVol)
            
        #if mcs == 0: self.import_chem_field("oxygen.field",self.field.Oxygen)
            
        
    def on_stop(self):
        # this gets called each time user stops simulation
        
            
        
        return    
        
    def import_chem_field(self, file,field):
            
        try:
            f = open(file)
            y = 0
            for line in f.readlines():
                values = line.split(" ")[:-1]
                for x in range(len(values)):
                    print(values[x])
                    field[x,y,0] = float(values[x])
                y += 1
        except IOError:
            print ("Could not open file for reading.")
            return
        
        
class CancerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)



    def start (self):
         for cell in self.cell_list_by_type(self.CANCER):
            cell.targetVolume = start_cancer_cell_volume
            cell.lambdaVolume = 2.0
            
            

    def step(self, mcs):
    
        # for cell in self.cell_list_by_type(self.CANCER):
            # cell.dict['nNeighbors'] = len(self.get_cell_neighbor_data_list(cell))
        
        field = self.field.Oxygen
    
        #Growth

        for cell in self.cell_list_by_type(self.CANCER):
            oxygenValue = field[int(cell.xCOM),int(cell.yCOM), int(cell.zCOM)]

            if oxygenValue > cancer_min_oxygen_growth and cell.volume < max_cancer_cell_volume:
                cell.targetVolume += cancer_cell_growth_rate*oxygenValue
        
            if oxygenValue < apoptosis_min_oxgen_level:
                if np.random.rand() < apoptosis_prob:
                    cell.type = self.CANCERAPOPTOSIS

        
        #Shrink
        for cell in self.cell_list_by_type(self.CANCERAPOPTOSIS):
            cell.targetVolume += apoptosis_shrink_rate

        for cell in self.cell_list_by_type(self.CANCERNECROSIS):
            cell.targetVolume += necrosis_shrink_rate

        
        self.perform_intervention(700,35,mcs)
        #self.perform_intervention(2500,35,mcs)
        
        #Dying by TCells
        for cell in self.cell_list_by_type(self.CANCER):
            touching_tcell = False
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.ACTTCELL:
                        touching_tcell = True
                        break
            if touching_tcell:
                if np.random.rand() < apoptosis_by_tcell_prob:
                    cell.type = self.CANCERAPOPTOSIS
        
        
    def perform_intervention(self, start, duration, mcs):
        if mcs > start and mcs < start+duration:
            for cell in self.cell_list_by_type(self.CANCER):
                if np.random.rand() < necrosis_prob:
                    cell.type = self.CANCERNECROSIS
        


        
class TcellSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        
        

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        
        
        #Tcell born
        damp_field = self.field.DAMP
        il2_field = self.field.IL2
        tgfb_field = self.field.TGFB
        oxygen_field = self.field.Oxygen
        
        
        offset = tcell_start_width+50
        tcel_pop_points  =  [(x,0+offset)for x in range(0+offset,self.dim.x-offset,tcell_start_width)]
        tcel_pop_points += [(x,self.dim.y-1-offset)for x in range(0+offset,self.dim.x-offset,tcell_start_width)]
        tcel_pop_points += [(0+offset,y)for y in range(0+offset,self.dim.y-offset,tcell_start_width)]
        tcel_pop_points += [(self.dim.x-offset,y)for y in range(0+offset,self.dim.y-offset,tcell_start_width)]
        
        
        for x,y in tcel_pop_points:
            if damp_field[x,y,0] + il2_field[x,y,0] - tgfb_field[x,y,0]  > il2_recruiting_threshold:
                if np.random.rand() < tcell_born_prob:
                    cell = self.potts.createCell()
                    cell.type = self.ACTTCELL
                    pt=CompuCell.Point3D()
                    for i in range(0,tcell_start_width):
                        for j in range(0,tcell_start_width):
                            pt.x = x+i 
                            pt.y = y+j
                            self.cellField.set(pt,cell)
                    cell.targetVolume = tcell_max_volume
                    cell.lambdaVolume = tcell_lambda_volume
                    
                    
        #Tcell IL2 secretion
        
        secretor = self.get_field_secretor("IL2")
        for cell in self.cell_list_by_type(self.ACTTCELL):
            release_IL2 = False
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.CANCER:
                        release_IL2 = True
                        break
                        
            if release_IL2: secretor.secreteOutsideCellAtBoundary(cell, il2_secretion)   
            
            
        #Tcell missing coffee
        for cell in self.cell_list_by_type(self.ACTTCELL):
            if tgfb_field[cell.xCOM,cell.yCOM,0]>il2_field[cell.xCOM,cell.yCOM,0]:
                if np.random.rand() < no_coffee_prob:
                    cell.type = self.NOCOFFEETCELL
                                  
         
        #Tcell Hipoxia
        for cell in self.cell_list_by_type(self.TREG):     
            oxygenValue = oxygen_field[int(cell.xCOM),int(cell.yCOM), int(cell.zCOM)]       
            if oxygenValue < tcell_min_oxigen_level:
                if np.random.rand() < apoptosis_prob:
                    cell.type = self.NOCOFFEETCELL    
   
   
        #Tcell shrinking
        for cell in self.cell_list_by_type(self.NOCOFFEETCELL):
            cell.targetVolume += no_coffee_shrink_rate
        

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return



class TregSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        
        

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        
        
        #Tcell born
        
        il2_field = self.field.IL2
        damp_field = self.field.DAMP
        tgfb_field = self.field.TGFB
        oxygen_field = self.field.Oxygen
        
        offset = tcell_start_width+50
        tcel_pop_points  =  [(x,0+offset)for x in range(0+offset,self.dim.x-offset,tcell_start_width)]
        tcel_pop_points += [(x,self.dim.y-1-offset)for x in range(0+offset,self.dim.x-offset,tcell_start_width)]
        tcel_pop_points += [(0+offset,y)for y in range(0+offset,self.dim.y-offset,tcell_start_width)]
        tcel_pop_points += [(self.dim.x-offset,y)for y in range(0+offset,self.dim.y-offset,tcell_start_width)]
        
        
        for x,y in tcel_pop_points:
            if damp_field[x,y,0] + il2_field[x,y,0] - tgfb_field[x,y,0] > il2_recruiting_threshold:
                if np.random.rand() < treg_born_prob:
                    cell = self.potts.createCell()
                    cell.type = self.TREG
                    pt=CompuCell.Point3D()
                    for i in range(0,tcell_start_width):
                        for j in range(0,tcell_start_width):
                            pt.x = x+i 
                            pt.y = y+j
                            self.cellField.set(pt,cell)
                    cell.targetVolume = tcell_max_volume
                    cell.lambdaVolume = tcell_lambda_volume
                    
                    
        #Tcell IL2 secretion
        
        secretor = self.get_field_secretor("TGFB")
        for cell in self.cell_list_by_type(self.TREG):
            release_TGFB = False
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.CANCER:
                        release_TGFB = True
                        break
                        
            if release_TGFB: secretor.secreteOutsideCellAtBoundary(cell, tgfb_secretion)   
            
        
        #Treg Hipoxia
        for cell in self.cell_list_by_type(self.TREG):
            oxygenValue = oxygen_field[int(cell.xCOM),int(cell.yCOM), int(cell.zCOM)]            
            if oxygenValue < treg_min_oxigen_level:
                if np.random.rand() < apoptosis_prob:
                    cell.type = self.NOCOFFEETCELL        
            
        if mcs%100 == 0:
            for cell in self.cell_list_by_type(self.TREG):
                if np.random.rand() < apoptosis_prob:
                    cell.type = self.NOCOFFEETCELL      
        

    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return




class CancerMitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        #max_neighbors = np.max([cell.dict['nNeighbors'] for cell in self.cell_list])
        #print(max_neighbors)
        
        field = self.field.Oxygen

        cells_to_divide=[]
        for cell in self.cell_list_by_type(self.CANCER):
            oxygenValue = field[int(cell.xCOM),int(cell.yCOM), int(cell.zCOM)]
            if oxygenValue >= 50 and cell.volume >= 100:
                if np.random.rand() < cancer_cell_mitosis_prob:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            #Other valid options
            #self.divide_cell_orientation_vector_based(cell,1,1,0)
            #self.divide_cell_along_major_axis(cell)
            #self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0                  

        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        # if self.parent_cell.type==1:
            # self.child_cell.type=2
        # else:
            # self.child_cell.type=1
        
        # for cell in self.cell_list:
            # cell.targetVolume += 1
