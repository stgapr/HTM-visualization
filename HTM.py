# coding=utf-8
from visual import *
import random
from random import uniform, randint,choice
import math
import time
import os
import sys
import weakref
OPEN_VISUAL=True

scene.forward = (-0.25, -0.25, -1)
class Encoder:

    # 记录NYC高低温
    minval = -30
    maxval = 120

    # 输出长度参数
    W = 11     # Number of 1's (必须是奇数）
    N = 400    # bit数量 (must be > w). 最好是平方数

    # 提取参数
    inputrange = float(maxval - minval)
    halfwidth = (W - 1)/2
    padding = halfwidth
    resolution = inputrange / (N - 2 * halfwidth)


    def encode(self, integer):
        """在 SDR中编码一个整数.

        @param integer: 要编码的整数
        @returns:  Encoder.W 编码后的bit列.
        """
        assert isinstance(integer, int), \
            '应该输入integer，但输入了{}'.format(type(integer))
        assert integer >= self.minval, \
            '输入{} 低于最小值 {}'.format(integer, self.minval)
        assert integer < self.maxval, \
            '输入 {} 高于最大值 {}'.format(integer, self.maxval)

        output = [0 for _ in range(self.N)]

        # 计算起终索引
        centerbin = integer - self.minval + self.resolution//2
        centerbin = centerbin // self.resolution + self.padding
        minbin = round(centerbin - self.halfwidth)
        maxbin = round(minbin + 2 * self.halfwidth)
        print("MINBIN,MAXBIN",minbin,maxbin)

        assert minbin >= 0
        assert maxbin < self.N
        for i in range(minbin, maxbin + 1):
            output[i] = 1
        #print(output)
        return output


def get_bit_at(sparse, x, y, z=0):

    side = int(math.sqrt(len(sparse)))
    return sparse[side * x + y]


class Synapse(object):

    THRESHOLD = 0.6
    """假如self.permanence > threshold, 则树突可连通。"""

    PERMDELTA = 0.05
    """学习率"""

    def __init__(self, parent,inputpos=(0, 0, 0), outputpos=(0,0,0),perm=1,connect_to_cell=None):
        self.permanence = perm
        self.inputpos = inputpos
        self.outputpos=outputpos
        self.parentDendrite=parent
        self.connect_to_cell=connect_to_cell

    def visual_drawsynapse(self):
        self.visual_synapse = curve(ends=[self.inputpos, self.outputpos], radius=0.01,color=(1,1,1))
        self.visual_synapse.pos = self.visual_synapse.ends
        self.parentDendrite.visual_syns.append(self.visual_synapse)

    def is_valid(self):
        if self.permanence > Synapse.THRESHOLD:
            self.visual_synapse.color=(0.3,1,0.4) #turn red
            return True
        else:
            self.visual_synapse.color=(1,1,1) #turn white
            return False
        return self.permanence > Synapse.THRESHOLD

    def increase_perm(self):
        self.permanence = min(self.permanence + Synapse.PERMDELTA, 1.0)

    def decrease_perm(self):
        self.permanence = max(0.0, self.permanence - Synapse.PERMDELTA)


class Dendrite(object):
    """
    细胞树突 (有时候被整个column共享)。  追踪现有的突触和学习时可能生长出突触的神经元
    """

    NPOTENTIAL = 10
    """该树突的潜在突触数量"""


    def __init__(self, parent,pos=(0, 0, 0)):
        self.visual_ddrs = []  # dendrites
        self.visual_syns = []
        self.pos=pos
         
        self.parentColumn=parent

    def visual_drawdendrite(self):
        if OPEN_VISUAL:
            self.visual_dendrite = cylinder(pos=self.pos,radius=0.12 )
            self.visual_dendrite.color = self.visual_dendrite.icolor = (0.2, 0.2, 0.2) #dark grey

            self.visual_dendrite.axis = (0, 0, 4)
            self.visual_ddrs.append(self.visual_dendrite)

            #实例化Synapses（绘制在最后完成）
    def makeSynapses(self):
        connect_to_cells=[]
        for i in range(Dendrite.NPOTENTIAL):
            cell=Region.instances[0].choice_a_cell()
            connect_to_cells.append(cell)
            self.synapses = [Synapse(self,inputpos=cell.pos,outputpos=self.pos,connect_to_cell=cell) ]

    

class Cell(object):
    """
    单个计算单元  ，持续追踪其树突部分
    """

    def __init__(self, proximaldendrite,pos):

        self.distal = []
        self.proximaldendrite = proximaldendrite
        self.pos=pos
        self._flag=False


        self.state = Column.INACTIVE
        if OPEN_VISUAL:self.visual_drawcell()

    def visual_drawcell(self):
        self.visual_cell=sphere(pos=self.pos,radius=0.3)
        self.visual_cell.color=self.visual_cell.icolor=(0.3,0.3,0.3) #grey
        self.visual_cell.axis=(0,0,4)
        #self.visual_cells.append(self.visual_cell)
        self._flag=True
    '''
    def color(self,color):
        print("setclolor")
        if  self._flag==True:
            self.visual_cell.color=self.visual_cell.icolor=color
    '''

class Column(object):
    """
   一列Cells. 一个column包含的所有cells 共享一个单 一的proximaldendrite段。
    """

    NCELLS = 4
    """Column里的cells数量"""

    INACTIVE = 0
    ACTIVE = 1
    PREDICTIVE = 2

    def __init__(self, parent,pos=(0, 0, 0)):
        self.parentRegion=parent
        self.proximaldendrite = Dendrite(self,pos)
        self.pos=pos
        self.cells=[]
        self.visual_cells=[] #cells for visualization
        for i in range(self.NCELLS):
            self.cells.append(Cell(self.proximaldendrite,pos=(pos[0],pos[1],i)))
        self.boost = 1.0  # 这个数应该适当降低，假如这个column太活跃

    def num_active_cells(self):
        """返回活跃细胞的当前数量"""
        return sum([cell.state == Column.ACTIVE for cell in self.cells])

    def num_predictive_cells(self):    	
        """返回预测细胞的当前数量"""
        return sum([cell.state == Column.PREDICTIVE for cell in self.cells])

    @property
    def state(self):
        """只要它含有的cells中的一个被激活，那么这个Column就是激活的 """
        if Column.ACTIVE in [c.state for c in self.cells]:
             self.proximaldendrite.visual_dendrite.color=(1,0,0) #red
             return Column.ACTIVE
        else:
            self.proximaldendrite.visual_dendrite.color=(1,1,1) #white
            return Column.INACTIVE

    @state.setter
    def state(self, state):
        """设置该Column包含的所有细胞的状态."""
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



class Region(object):

    NCOLUMNS = 400  # needs to be a perfect square for now!
    side_len=math.sqrt(NCOLUMNS/2) *2

    density = 2. / 100
    """活跃column的稀疏度"""
    instances=[]

    def __init__(self,name=None):
        self.__class__.instances.append(weakref.proxy(self))
        self.name=name

        side_len=math.sqrt(Region.NCOLUMNS)/2
        self.columns = [Column(self,(index // 20, index % 20, 0))
                        for index, _ in enumerate(range(self.NCOLUMNS))]
        self.last_inactivecols=[]
        self.last_activecols=[]

        #现在绘制syn
     
        for column in self.columns:
            if OPEN_VISUAL: column.proximaldendrite.visual_drawdendrite()    
            column.proximaldendrite.makeSynapses()

            for syn in column.proximaldendrite.synapses:
                syn.visual_drawsynapse()
      

    def choice_a_cell(self):
        return choice(choice(self.columns).cells)



    @property
    def column_matrix(self):
        side = int(math.sqrt(Region.NCOLUMNS))

        return [[self.columns[row * side + col] for col in range(side)]
                for row in range(side)]

    def process(self, sparse):

        numvalid = {}

        for col in self.columns: #算一下有效刺激本column的synapses数量
            #SIMPLE METHOD
            numvalid[col] = len([syn for syn in col.proximaldendrite.synapses
                                 if syn.is_valid() and
                                 get_bit_at(sparse, *syn.inputpos) == 1])
            numvalid[col] *= col.boost

        srtd = sorted(self.columns, key=numvalid.get)  # 一般来说这个numvalid>15
        index = Encoder.N - round(Encoder.N * self.density)  # ???

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

            # 判断那些active的syn连接到的细胞是否足够进入predictive状态
            
            for syn in activesyn:
                print("go predicitve!",syn.connect_to_cell.pos)
                syn.connect_to_cell.visual_cell.color=(0.6,0.8,0) #yello
                syn.connect_to_cell.state = Column.PREDICTIVE

            #现在预测位已经准备好，开始看哪些细胞能够被column带活。           
            

       
           





        return

    def __str__(self):
        rules = {4: "■", 3: "X", 2: "x", 1: "o", 0: "·"}
        strings = [''.join([rules[col.num_predictive_cells()]
                            for col in row]) for row in self.column_matrix]
        return '\n'.join(strings)


def main():
    scene.forward = (0 ,1 ,-1)
    """用随机生成的数据跑HTM"""
    random.seed()

    encoder = Encoder()
    region = Region()
    drag = None
    step=1
    while True:
        value = int(abs(math.sin(0.1 * step) * 0.5 *(encoder.maxval - encoder.minval)))  # abs(sin)
        #value=random.randint(encoder.minval+1,encoder.maxval-1)
        data = encoder.encode(value)
        region.process(data) #等下
        print(region)
       
        rate(100)
        #time.sleep(1)
        os.system("cls")
        step=step+1
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




if __name__ == "__main__":
    main()
