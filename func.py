#!/usr/bin/env python
# coding: utf-8

# In[7]:

class func:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def Hough_Trans(self):
        # Generate a dataframe for storage
        df = pd.DataFrame({'x':[]})
        a = lambda i: self.x * math.cos(i * math.pi/180) + self.y * math.sin(i * math.pi/180)  # the distance to the line
        distance = [a(i) for i in range(181) ]
        df = df.append(pd.DataFrame({'%s'%self.x: distance}),sort = True)
        df.dropna(axis = 'columns')
        return df




