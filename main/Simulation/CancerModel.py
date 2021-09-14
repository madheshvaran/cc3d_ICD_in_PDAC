
from cc3d import CompuCellSetup
        


from CancerModelSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from CancerModelSteppables import CancerSteppable

CompuCellSetup.register_steppable(steppable=CancerSteppable(frequency=1))


from CancerModelSteppables import TcellSteppable

CompuCellSetup.register_steppable(steppable=TcellSteppable(frequency=1))


from CancerModelSteppables import TregSteppable

CompuCellSetup.register_steppable(steppable=TregSteppable(frequency=1))


from CancerModelSteppables import CancerMitosisSteppable

CompuCellSetup.register_steppable(steppable=CancerMitosisSteppable(frequency=1))

CompuCellSetup.run()
