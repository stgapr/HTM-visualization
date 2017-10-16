# coding=utf-8
from visual import *
from HTM import *
OPEN_VISUAL=True #visualization interface won't show when put to False


class visualizableSynapse(Synapse):
   
    def __init__(self,parent,inputpos=(0, 0, 0), outputpos=(0,0,0),perm=1,connect_to_cell=None):
        Synapse.__init__(self,parent,inputpos=(0, 0, 0), outputpos=(0,0,0),perm=1,connect_to_cell=None)  
  
    def visual_drawsynapse(self,):
        '''
        draw the graph of synapse.
        '''
        self.visual_synapse = curve(ends=[self.inputpos, self.outputpos], radius=0.01,color=(1,1,1))
        self.visual_synapse.pos = self.visual_synapse.ends
        self.parentDendrite.visual_syns.append(self.visual_synapse)
    
    def is_valid(self):
        '''
        check if the synapse is valid. It's a binarization  function.
        '''
        if super(Synapse,self).is_valid()==True:
            self.visual_synapse.color=(0.3,1,0.4) #turn red
        else:
            self.visual_synapse.color=(1,1,1) #turn white

class visualizableDendrite(Dendrite):

    def __init__(self,parent,pos=(0, 0, 0)):
        Dendrite.__init__(self,parent,pos=(0, 0, 0))  
        self.visual_ddrs = []  # dendrites
        self.visual_syns = []
   

    def visual_drawdendrite(self):
        '''
        draw the graph of dendrite
        '''
        if not (hasattr(self, 'visual_ddrs') and hasattr(self,'visual_syns')):
            self.visual_ddrs=[]
            self.visual_syns=[]

        if OPEN_VISUAL:
            self.visual_dendrite = cylinder(pos=self.pos,radius=0.12 )
            self.visual_dendrite.color = self.visual_dendrite.icolor = (0.2, 0.2, 0.2) #dark grey

            self.visual_dendrite.axis = (0, 0, 4)
            self.visual_ddrs.append(self.visual_dendrite)

class visualizableCell(Cell):
 
    def __init__(self,proximaldendrite,pos):
        Dendrite.__init__(self,proximaldendrite,pos)
        if OPEN_VISUAL:
            self.visual_drawcell()
        

    def visual_drawcell(self):
        '''
        draw graph of this cell.
        '''
        self.visual_cell=sphere(pos=self.pos,radius=0.3)
        self.visual_cell.color=self.visual_cell.icolor=(0.3,0.3,0.3) #grey
        self.visual_cell.axis=(0,0,4)
        #self.visual_cells.append(self.visual_cell)
        self._flag=True

class visualizableColumn(Column):
 
    def __init__(self,parent,pos=(0, 0, 0)):
       Column.__init__(self,parent,pos=(0, 0, 0))  
        #self.visual_cells=[] #cells for visualization
 

    @property
    def state(self):
        if super(Column,self).property()==Column.ACTIVE:
            self.proximaldendrite.visual_dendrite.color=(1,0,0) #red
        else:
            self.proximaldendrite.visual_dendrite.color=(1,1,1) #white
   
    @state.setter
    def state(self, state):
        """Update the state of each Cell in this Column."""
        if state == Column.ACTIVE:
            predictive = {
                c for c in self.cells if c.state == Column.PREDICTIVE}

            if len(predictive) != 0:
                for cell in predictive:
                    cell.state = Column.ACTIVE                    
                    cell.visual_cell.color=(1,0,0) #red
                for cell in set(self.cells) - predictive:
                    cell.state = Column.INACTIVE
                    cell.visual_cell.color=(1,1,1) #white

            else:  # 没有predictive的
                for cell in self.cells:
                    cell.state = Column.ACTIVE
                    cell.visual_cell.color=(1,0,0) #red

        elif state == Column.INACTIVE:
            for cell in self.cells:
                cell.state = Column.INACTIVE
                cell.visual_cell.color=(1,1,1) #white

class visualizableRegion(Region): 


    def __init__(self,name=None):
        Region.__init__(self,name=None)

    def makeSynapses(self):
        for column in self.columns:
            if OPEN_VISUAL: column.proximaldendrite.visual_drawdendrite()    
            column.proximaldendrite.createSynapses()

            for syn in column.proximaldendrite.synapses:
                syn.visual_drawsynapse()
    
    def predict(self,activesyn):
            # those cells connected to active synapses turn to be predicitve             
            for syn in activesyn:
                syn.connect_to_cell.visual_cell.color=(0.6,0.8,0) #yellow
                syn.connect_to_cell.state = Column.PREDICTIVE
def main():
    
    """用随机生成的数据跑HTM"""
    random.seed()
    encoder = Encoder()
    region = visualizableRegion()
    
    step=1
    while True:
        step=loop(step=step,encoder=encoder,region=region)  
        processEvent()

if __name__ == "__main__":
    scene.forward = (0 ,1 ,-1)
    drag = None
    rate(100)
    main()

