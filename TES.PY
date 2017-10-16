class BaseClass():
    def __init__(self):
        self.variable="123"

class newClass(BaseClass):
    def __init__(self):
        super(BaseClass,self).__init__()
        self.b="bbb"
        print(self.variable)
  


newinstance=newClass()
print("1111111111111",newinstance.variable)