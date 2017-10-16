# coding=utf-8

import random
from random import uniform, randint,choice
import math
import time
import os
import sys
import weakref

# upper-lower limit  of input data 

class Encoder:
    minval = -30
    maxval = 120

    W = 11     # Length of Sparse representation 
    N = 400    # The amount of bits. Better to be a squared number (must be > w).

    # 提取参数
    inputrange = float(maxval - minval)
    halfwidth = (W - 1)/2
    padding = halfwidth
    resolution = inputrange / (N - 2 * halfwidth)
    def encode(self, integer):
            """Encode a number and turn it to be a sparse representation,

            @param integer: Data (Integer) to be encode
            @returns:  Sparsed bit array which is encoded.
            """
            assert isinstance(integer, int), \
                'There should be an integer，but got{}'.format(type(integer))
            assert integer >= self.minval, \
                'Input{} lower than minimum {}'.format(integer, self.minval)
            assert integer < self.maxval, \
                'Input {} higher than maximum {}'.format(integer, self.maxval)

            output = [0 for _ in range(self.N)]

            # Calculate how data should be interpret to a sparse bit array.
            centerbin = integer - self.minval + self.resolution//2
            centerbin = centerbin // self.resolution + self.padding 
            minbin = round(centerbin - self.halfwidth)
            maxbin = round(minbin + 2 * self.halfwidth)

            assert minbin >= 0
            assert maxbin < self.N
            for i in range(minbin, maxbin + 1):
                output[i] = 1
            return output


def get_bit_at(sparse, x, y, z=0):
        """Get a specify bit from a sparse bit.

        @param sparse: Data (Integer) to be encode
        @param x,y: x,y indicate  the place where the bit should be in 2-dimensions. 
        @param z: z is a preserved parameter, useless now.
        @returns:  the bit value where you appointed .
        """
        side = int(math.sqrt(len(sparse)))
        return sparse[side * x + y]

##############
# HTM classes
##############
#
class Synapse(object):

    THRESHOLD = 0.6
    """If self.permanence > threshold,  the synapse is valid"""

    PERMDELTA = 0.05
    """learning rate."""

    def __init__(self, parent,inputpos=(0, 0, 0), outputpos=(0,0,0),perm=1,connect_to_cell=None):
        """Get a specify bit from a sparse bit.

        @param parent: the dendrite it belongs to / connect from
        @param inputpos: x,y indicate  the place where the bit should be in 2-dimensions. 
        @param outputpos: z is a preserved parameter, useless now.
        @param perm:
        @returns:  the bit value where you appointed .
        """
        self.permanence = perm
        self.inputpos = inputpos
        self.outputpos=outputpos
        self.parentDendrite=parent
        self.connect_to_cell=connect_to_cell


    def is_valid(self):
        '''
        check if the synapse is valid. It's a binarization  function.
        '''
        if self.permanence > Synapse.THRESHOLD:
            return True
        else:
            return False

    def increase_perm(self):
        '''
        increase the permanence value of synapse, by a learning rate (PERMDELTA)
        '''
        self.permanence = min(self.permanence + Synapse.PERMDELTA, 1.0)

    def decrease_perm(self):
        '''
        decrease the permanence value of synapse, by a learning rate (PERMDELTA)
        '''
        self.permanence = max(0.0, self.permanence - Synapse.PERMDELTA)


class Dendrite(object):
    """
    Stems from cells (sometimes shared by a whole column).  Keeps track of
    current Synapses and of cells that might form synapses while learning.
    """

    NPOTENTIAL = 10
    """The number of potential synapses for this Dendrite."""


    def __init__(self, parent,pos=(0, 0, 0)):

        self.pos=pos         
        self.parentColumn=parent


    def createSynapses(self):
        '''
        create synapses which start from ths dendrite, to a random-chosen cell.
        '''
        connect_to_cells=[]
        for i in range(Dendrite.NPOTENTIAL):
            cell=Region.instances[0].choice_a_cell()
            connect_to_cells.append(cell)
            self.synapses = [Synapse(self,inputpos=cell.pos,outputpos=self.pos,connect_to_cell=cell) ]

    

class Cell(object):
    """
    A single computational unit.  Keeps track of its dendrite segments,
    which determine where a Cell gets its input from, and its state.
    """

    def __init__(self, proximaldendrite,pos):
        """
        A Column has one shared Dendrite to each of its cells.
        @param proximal: a Dendrite, shared among Cells in the same Column.
        @param pos:  coordinates of this cell

        """
        self.distal = []
        self.proximaldendrite = proximaldendrite
        self.pos=pos
        self._flag=False
        self.state = Column.INACTIVE
        

class Column(object):
    """
    An array of Cells.  All the cells in one column share a single (proximal
    ) dendrite segment.
    """

    NCELLS = 4
    """Number of cells in this Column."""

    INACTIVE = 0
    ACTIVE = 1
    PREDICTIVE = 2

    def __init__(self, parent,pos=(0, 0, 0)):
        """
        @param parent: the parent region it belongs to.
        @param pos: coordinates of this column (the third parameter always =0)
        """
        self.parentRegion=parent
        self.proximaldendrite = Dendrite(self,pos)
        self.pos=pos
        self.cells=[]

        '''add cell to this column'''
        for i in range(self.NCELLS):
            self.cells.append(Cell(self.proximaldendrite,pos=(pos[0],pos[1],i)))
        self.boost = 1.0  # 这个数应该适当降低，假如这个column太活跃

    def num_active_cells(self):
        """Return the current number of active cells."""
        return sum([cell.state == Column.ACTIVE for cell in self.cells])

    def num_predictive_cells(self):    	
        """Return the current number of predictive cells"""
        return sum([cell.state == Column.PREDICTIVE for cell in self.cells])

    @property
    def state(self):
        """
        A Column is active whenever at least one of its Cells is active.
        Therefore. So we must set the status of cells, then calculate the status
        next step will be in each step, by this function.
        """

        if Column.ACTIVE in [c.state for c in self.cells]:
             return Column.ACTIVE
        else:
            return Column.INACTIVE

    @state.setter
    def state(self, state):
        """Update the state of each Cell in this Column."""
        if state == Column.ACTIVE:
            predictive = {
                c for c in self.cells if c.state == Column.PREDICTIVE}

            if len(predictive) != 0:
                for cell in predictive:
                    cell.state = Column.ACTIVE                    

                for cell in set(self.cells) - predictive:
                    cell.state = Column.INACTIVE


            else:  # 没有predictive的
                for cell in self.cells:
                    cell.state = Column.ACTIVE


        elif state == Column.INACTIVE:
            for cell in self.cells:
                cell.state = Column.INACTIVE
 



class Region(object):
    """A Region is a plain full of Columns."""

    NCOLUMNS = 400  
    """Number of columns in a Region."""
    side_len=math.sqrt(NCOLUMNS/2) *2

    density = 2. / 100
    """The sparsity of the active columns."""
    instances=[]

    def __init__(self,name=None):
        self.__class__.instances.append(weakref.proxy(self))
        #some hack for getting instance of Region in other class when program's running.
        self.name=name
        side_len=math.sqrt(Region.NCOLUMNS)/2
        self.columns = [Column(self,(index // 20, index % 20, 0))
                        for index, _ in enumerate(range(self.NCOLUMNS))]
        #store the array of active/inactive columns, and update them when next step start,
        #in order to avoid 
        self.last_inactivecols=[]
        self.last_activecols=[]

        #when everything done, we could draw the graph of synapses now 
        self.makeSynapses()

    def makeSynapses(self):
        for column in self.columns:
            column.proximaldendrite.createSynapses()

    def choice_a_cell(self):
        return choice(choice(self.columns).cells)

    @property
    def column_matrix(self):
        """Return the Columns as a 2D array."""
        side = int(math.sqrt(Region.NCOLUMNS))

        return [[self.columns[row * side + col] for col in range(side)]
                for row in range(side)]

    def process(self, sparse):
        """
        Process some input data (runs spatial and temporal pooler),
        and update the column status  of last step, then start a new step:
        check the activity and predictivity, increase/decrease permenance of synapses,
        """
        numvalid = {}

        for col in self.columns: 
           #columns compete to be active by comparing the amount of valid synapses provoked them
            numvalid[col] = len([syn for syn in col.proximaldendrite.synapses
                                 if syn.is_valid() and
                                 get_bit_at(sparse, *syn.inputpos) == 1])
            numvalid[col] *= col.boost

        srtd = sorted(self.columns, key=numvalid.get)  
        # according to white paper of HTM, numvalid which >15 is ideal.
        index = Encoder.N - round(Encoder.N * self.density)  # 

        for col in self.last_activecols:
            col.state = Column.ACTIVE   
         
        for col in self.last_inactivecols: 
            col.state = Column.INACTIVE

        self.last_inactivecols = srtd[:index]
        self.last_activecols = srtd[index:]

        for winner in self.last_activecols:
            print("active cols",winner.pos)

        for col in self.last_activecols:
            activesyn = {syn for syn in col.proximaldendrite.synapses
                         if get_bit_at(sparse, *syn.inputpos) == 1}
            for synapse in activesyn:
                synapse.increase_perm()
            inactivesyn = set(col.proximaldendrite.synapses) - activesyn
            for synapse in inactivesyn:
                synapse.decrease_perm()
            self.predict(activesyn)

    def predict(self,activesyn):
        # those cells connected to active synapses turn to be predicitve             
        for syn in activesyn:
            #syn.connect_to_cell.visual_cell.color=(0.6,0.8,0) #yellow
            syn.connect_to_cell.state = Column.PREDICTIVE


        return

    def __str__(self):
        rules = {4: "■", 3: "X", 2: "x", 1: "o", 0: "·"}
        strings = [''.join([rules[col.num_predictive_cells()]
                            for col in row]) for row in self.column_matrix]
        return '\n'.join(strings)




def loop(step=1,encoder=Encoder(),region=Region()):
        value = int(abs(math.sin(0.1 * step) * 0.5 *(encoder.maxval - encoder.minval)))  
        '''abs(sin)'''
        #value=random.randint(encoder.minval+1,encoder.maxval-1)
        '''random data (meaningless)'''
        data = encoder.encode(value)
        region.process(data) 
        print(region) 
        
        time.sleep(0.2)
        os.system("cls") #clear the output in command window.
        step=step+1
        return step



def processEvent():
        pass
        #process of mouse dragging. some bug now''
        '''
        if scene.mouse.events:
            visual_synapse = scene.mouse.getevent()

            if drag and (visual_synapse.drop or visual_synapse.click):   # drop the selected object
                newpos = visual_synapse.project(
                    normal=scene.up, d=yoffset) + offset
                if abs(newpos.x) <= Region.side_len and abs(newpos.z) <= Region.side_len:
                    drag.pos = newpos
                drag.color = drag.icolor
                drag = None
            # pick up the object
            elif visual_synapse.pick and hasattr(visual_synapse.pick, "icolor"):
                drag = visual_synapse.pick
                drag.color = color.white
                yoffset = visual_synapse.pickpos.y
                offset = drag.pos - \
                    visual_synapse.project(normal=scene.up, d=yoffset)

        if drag:
            newpos = scene.mouse.project(normal=scene.up, d=yoffset) + offset
            if abs(newpos.x) <= Region.side_len and abs(newpos.z) <= Region.side_len:
                drag.pos = newpos

  
        for column in region.columns:
            for visual_syn in column.proximaldendrite.visual_syns:
                 visual_syn.pos = visual_syn.ends 
        '''
def main():
    
    """用随机生成的数据跑HTM"""
    random.seed()
    encoder = Encoder()
    region = Region()
    
    step=1
    while True:
        step=loop(step=step,encoder=encoder,region=region)  
        processEvent()

if __name__ == "__main__":
    main()
