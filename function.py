#!/usr/bin/env python
# coding: utf-8

from math import pi
from math import sin
from math import cos

class func:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def Hough_Trans(self):
        # the distance to the line
        a = lambda i: self.x * cos(i * pi/180) + self.y * sin(i * pi/180)  
        distance = [a(i) for i in range(181)]
        return distance
