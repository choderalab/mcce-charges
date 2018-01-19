#!/usr/bin/python
import math


class VECTOR:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z


def ddvv(v1, v2):
    dx = v2.x - v1.x
    dy = v2.y - v1.y
    dz = v2.z - v1.z
    return dx*dx+dy*dy+dz*dz


def vector_normalize(v):
    d = math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
    vn = VECTOR(0, 0, 0)
    if (d > 1.0e-20):
        vn.x = v.x/d
        vn.y = v.y/d
        vn.z = v.z/d

    return vn


def vector_vminusv(v1, v2):
    z = VECTOR()
    z.x = v1.x - v2.x
    z.y = v1.y - v2.y
    z.z = v1.z - v2.z
    return z


def avv(v1, v2):
    v1 = vector_normalize(v1)
    v2 = vector_normalize(v2)
    t = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
    if t>1.0:
        t =1.0
    elif t< -1.0:
        t = -1.0

    return math.acos(t)
